# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the Isovar phase-resolver integration (#269).

These don't import Isovar. A tiny in-memory ``StubAssemblyProvider``
implements :class:`IsovarAssemblyProvider`'s three methods and stands
in for the real thing — verifying the plumbing without tying the test
suite to Isovar's release cadence.
"""

import pytest
from pyensembl import cached_release

from varcode import (
    IsovarAssemblyProvider,
    IsovarPhaseResolver,
    MutantTranscript,
    Variant,
    VariantCollection,
    apply_phase_resolver_to_effects,
)


ensembl_grch38 = cached_release(81)
CFTR_ID = "ENST00000003084"
BRCA1_ID = "ENST00000357654"


class StubAssemblyProvider:
    """In-memory provider for tests. Maps
    ``(variant, transcript.id) -> (variants_on_contig, mutant_transcript)``.
    """

    def __init__(self, contigs):
        self._contigs = {
            (v, tid): (partners, mt)
            for (v, tid), (partners, mt) in contigs.items()
        }

    def has_contig(self, variant, transcript):
        return (variant, transcript.id) in self._contigs

    def variants_in_contig(self, variant, transcript):
        entry = self._contigs.get((variant, transcript.id))
        return entry[0] if entry is not None else ()

    def mutant_transcript(self, variant, transcript):
        entry = self._contigs.get((variant, transcript.id))
        return entry[1] if entry is not None else None


# --------------------------------------------------------------------
# Protocol conformance
# --------------------------------------------------------------------


def test_stub_provider_matches_protocol():
    """Sanity check: the stub used throughout these tests satisfies
    the :class:`IsovarAssemblyProvider` runtime-checkable Protocol."""
    provider = StubAssemblyProvider({})
    assert isinstance(provider, IsovarAssemblyProvider)


# --------------------------------------------------------------------
# IsovarPhaseResolver basics
# --------------------------------------------------------------------


def _cftr():
    return ensembl_grch38.transcript_by_id(CFTR_ID)


def test_resolver_source_is_isovar():
    resolver = IsovarPhaseResolver(StubAssemblyProvider({}))
    assert resolver.source == "isovar"


def test_in_cis_true_when_both_on_same_contig():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    contigs = {
        (v1, transcript.id): ((v1, v2), None),
        (v2, transcript.id): ((v1, v2), None),
    }
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))
    assert resolver.in_cis(v1, v2, transcript=transcript) is True


def test_in_cis_false_when_on_different_contigs():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    # v1 is on a contig, v2 is on a *different* contig. Each contig
    # only contains its own variant.
    contigs = {
        (v1, transcript.id): ((v1,), None),
        (v2, transcript.id): ((v2,), None),
    }
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))
    assert resolver.in_cis(v1, v2, transcript=transcript) is False


def test_in_cis_none_when_no_evidence():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    resolver = IsovarPhaseResolver(StubAssemblyProvider({}))
    assert resolver.in_cis(v1, v2, transcript=transcript) is None


def test_in_cis_handles_one_sided_evidence():
    """v1 has no contig but appears in v2's contig → cis by virtue
    of being on the same physical molecule as v2."""
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    contigs = {
        (v2, transcript.id): ((v1, v2), None),
        # Deliberately no entry keyed on v1 — a real Isovar result
        # might only build one assembly per somatic variant, and
        # germline variants only show up inside it.
    }
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))
    assert resolver.in_cis(v1, v2, transcript=transcript) is True


def test_phased_partners_returns_cis_set():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    contigs = {(v1, transcript.id): ((v1, v2), None)}
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))
    partners = resolver.phased_partners(v1, transcript)
    assert v2 in partners


# --------------------------------------------------------------------
# Effect-attachment flow — the headline integration
# --------------------------------------------------------------------


def _stub_mutant_transcript(transcript):
    """Build a MutantTranscript that's clearly distinguishable from
    anything the DNA-only annotator would construct: tag the
    annotator_name so tests can tell the Isovar payload landed."""
    return MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="ACGT" * 50,
        mutant_protein_sequence="M" + "A" * 99,
        annotator_name="isovar",
    )


def test_effect_gets_isovar_mutant_transcript_attached():
    """End-to-end: pass a ``phase_resolver`` into
    ``Variant.effects()`` and verify the returned effect's
    ``mutant_transcript`` is the Isovar-provided one, not None."""
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    stub_mt = _stub_mutant_transcript(transcript)
    contigs = {(variant, transcript.id): ((variant,), stub_mt)}
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))

    effects = variant.effects(phase_resolver=resolver)
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is stub_mt
    # The assembled protein is what's attached — not varcode's
    # reference-diff inference.
    assert target.mutant_transcript.annotator_name == "isovar"


def test_effects_without_resolver_have_no_mutant_transcript_attached():
    """Back-compat: without a resolver, point-variant effects carry
    the class-default ``mutant_transcript = None`` — nothing has
    been populated."""
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    effects = variant.effects()
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is None


def test_variant_collection_effects_accepts_phase_resolver():
    """Same plumbing on the collection. Returns an EffectCollection
    whose covered effects carry the Isovar MutantTranscript."""
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    mt1 = _stub_mutant_transcript(transcript)
    contigs = {(v1, transcript.id): ((v1, v2), mt1)}
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))

    collection = VariantCollection([v1, v2])
    effects = collection.effects(phase_resolver=resolver)

    v1_effects = [e for e in effects if e.variant == v1
                  and getattr(e, "transcript", None) is transcript]
    assert v1_effects
    assert v1_effects[0].mutant_transcript is mt1
    # v2 has no contig of its own — mutant_transcript stays None.
    v2_effects = [e for e in effects if e.variant == v2
                  and getattr(e, "transcript", None) is transcript]
    for e in v2_effects:
        assert e.mutant_transcript is None


def test_resolver_without_contig_leaves_effects_alone():
    """Resolver present but no contig for this (variant, transcript) →
    effect's ``mutant_transcript`` stays None."""
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    resolver = IsovarPhaseResolver(StubAssemblyProvider({}))
    effects = variant.effects(phase_resolver=resolver)
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is None


def test_apply_phase_resolver_helper_is_idempotent():
    """Calling apply_phase_resolver_to_effects twice with the same
    resolver leaves the same objects in place."""
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    mt = _stub_mutant_transcript(transcript)
    contigs = {(variant, transcript.id): ((variant,), mt)}
    resolver = IsovarPhaseResolver(StubAssemblyProvider(contigs))
    effects = variant.effects()
    apply_phase_resolver_to_effects(effects, resolver)
    apply_phase_resolver_to_effects(effects, resolver)
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is mt


def test_user_example_shape_works():
    """The canonical example from #269 works as-written: build a
    resolver, pass it to effects(), iterate and read
    effect.mutant_transcript.mutant_protein_sequence."""
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    stub_mt = _stub_mutant_transcript(transcript)
    contigs = {(variant, transcript.id): ((variant,), stub_mt)}

    isovar_results = StubAssemblyProvider(contigs)
    resolver = IsovarPhaseResolver(isovar_results)

    variants = VariantCollection([variant])
    effects = variants.effects(phase_resolver=resolver)

    found_protein = None
    for e in effects:
        t = getattr(e, "transcript", None)
        if t is transcript and isovar_results.has_contig(variant, t):
            found_protein = e.mutant_transcript.mutant_protein_sequence
            break
    assert found_protein == stub_mt.mutant_protein_sequence
