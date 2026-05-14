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

"""Tests for ReadPhaseResolver and the read-phasing protocols (#269).

These don't import any upstream tool. A tiny in-memory
``StubReadPhasingSource`` implements the two narrow Protocols and
stands in for the real thing (e.g. an Isovar adapter shipped by
openvax/isovar) — verifying the plumbing without tying the test suite
to any one upstream's release cadence.
"""

from pyensembl import cached_release

from varcode import (
    MutantTranscript,
    MutantTranscriptSource,
    ReadPhaseResolver,
    ReadPhasingSource,
    Variant,
    VariantCollection,
    apply_phase_resolver_to_effects,
)


ensembl_grch38 = cached_release(81)
CFTR_ID = "ENST00000003084"


class StubReadPhasingSource:
    """In-memory source for tests. Variant-keyed phasing channel +
    optional ``(variant, transcript)``-keyed mutant-transcript channel.
    """

    def __init__(self, phasing=None, mutant_transcripts=None):
        self._phasing = dict(phasing) if phasing else {}
        self._mts = dict(mutant_transcripts) if mutant_transcripts else {}

    def has_evidence(self, variant):
        return variant in self._phasing

    def partners_in_cis(self, variant):
        return self._phasing.get(variant, ())

    def mutant_transcript(self, variant, transcript):
        return self._mts.get((variant, transcript.id))


# --------------------------------------------------------------------
# Protocol conformance
# --------------------------------------------------------------------


def test_stub_satisfies_read_phasing_protocol():
    source = StubReadPhasingSource()
    assert isinstance(source, ReadPhasingSource)


def test_stub_satisfies_mutant_transcript_protocol():
    source = StubReadPhasingSource()
    assert isinstance(source, MutantTranscriptSource)


# --------------------------------------------------------------------
# ReadPhaseResolver basics
# --------------------------------------------------------------------


def _cftr():
    return ensembl_grch38.transcript_by_id(CFTR_ID)


def test_resolver_default_source_tag():
    resolver = ReadPhaseResolver(StubReadPhasingSource())
    assert resolver.source == "read_phasing"


def test_in_cis_true_when_both_observed_together():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    source = StubReadPhasingSource(phasing={v1: (v1, v2), v2: (v1, v2)})
    resolver = ReadPhaseResolver(source)
    assert resolver.in_cis(v1, v2, transcript=transcript) is True


def test_in_cis_false_when_observed_on_different_molecules():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    source = StubReadPhasingSource(phasing={v1: (v1,), v2: (v2,)})
    resolver = ReadPhaseResolver(source)
    assert resolver.in_cis(v1, v2, transcript=transcript) is False


def test_in_cis_none_when_no_evidence():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    resolver = ReadPhaseResolver(StubReadPhasingSource())
    assert resolver.in_cis(v1, v2, transcript=transcript) is None


def test_in_cis_handles_one_sided_evidence():
    """v1 has no direct evidence but appears in v2's partner set →
    cis. Matches the case where an upstream tool builds one entry per
    somatic variant and germline partners appear only inside it."""
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    source = StubReadPhasingSource(phasing={v2: (v1, v2)})
    resolver = ReadPhaseResolver(source)
    assert resolver.in_cis(v1, v2, transcript=transcript) is True


def test_phased_partners_returns_cis_set():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    source = StubReadPhasingSource(phasing={v1: (v1, v2)})
    resolver = ReadPhaseResolver(source)
    partners = resolver.phased_partners(v1, transcript)
    assert v2 in partners


def test_phased_partners_empty_without_evidence():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    resolver = ReadPhaseResolver(StubReadPhasingSource())
    assert resolver.phased_partners(v1, transcript) == ()


# --------------------------------------------------------------------
# MutantTranscript channel — optional; routed through wrapped source
# --------------------------------------------------------------------


def _stub_mutant_transcript(transcript):
    return MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="ACGT" * 50,
        mutant_protein_sequence="M" + "A" * 99,
        annotator_name="upstream_rna_tool",
    )


def test_mutant_transcript_returns_none_when_source_lacks_method():
    """A ReadPhasingSource that doesn't implement MutantTranscriptSource
    causes the resolver's mutant_transcript() to return None silently
    rather than raise."""

    class PhasingOnly:
        def has_evidence(self, variant):
            return False

        def partners_in_cis(self, variant):
            return ()

    transcript = _cftr()
    v = Variant("7", 117531100, "T", "A", ensembl_grch38)
    resolver = ReadPhaseResolver(PhasingOnly())
    assert resolver.mutant_transcript(v, transcript) is None


def test_effect_gets_observed_mutant_transcript_attached():
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    stub_mt = _stub_mutant_transcript(transcript)
    source = StubReadPhasingSource(
        phasing={variant: (variant,)},
        mutant_transcripts={(variant, transcript.id): stub_mt},
    )
    resolver = ReadPhaseResolver(source)

    effects = variant.effects(phase_resolver=resolver)
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is stub_mt
    assert target.mutant_transcript.annotator_name == "upstream_rna_tool"


def test_effects_without_resolver_have_no_mutant_transcript_attached():
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    effects = variant.effects()
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is None


def test_variant_collection_effects_accepts_phase_resolver():
    transcript = _cftr()
    v1 = Variant("7", 117531100, "T", "A", ensembl_grch38)
    v2 = Variant("7", 117531114, "G", "T", ensembl_grch38)
    mt1 = _stub_mutant_transcript(transcript)
    source = StubReadPhasingSource(
        phasing={v1: (v1, v2)},
        mutant_transcripts={(v1, transcript.id): mt1},
    )
    resolver = ReadPhaseResolver(source)

    collection = VariantCollection([v1, v2])
    effects = collection.effects(phase_resolver=resolver)

    v1_effects = [e for e in effects if e.variant == v1
                  and getattr(e, "transcript", None) is transcript]
    assert v1_effects
    assert v1_effects[0].mutant_transcript is mt1
    v2_effects = [e for e in effects if e.variant == v2
                  and getattr(e, "transcript", None) is transcript]
    for e in v2_effects:
        assert e.mutant_transcript is None


def test_resolver_without_evidence_leaves_effects_alone():
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    resolver = ReadPhaseResolver(StubReadPhasingSource())
    effects = variant.effects(phase_resolver=resolver)
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is None


def test_apply_phase_resolver_helper_is_idempotent():
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    mt = _stub_mutant_transcript(transcript)
    source = StubReadPhasingSource(
        phasing={variant: (variant,)},
        mutant_transcripts={(variant, transcript.id): mt},
    )
    resolver = ReadPhaseResolver(source)
    effects = variant.effects()
    apply_phase_resolver_to_effects(effects, resolver)
    apply_phase_resolver_to_effects(effects, resolver)
    target = next(e for e in effects if getattr(e, "transcript", None) is
                  transcript)
    assert target.mutant_transcript is mt


def test_user_example_shape_works():
    transcript = _cftr()
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    stub_mt = _stub_mutant_transcript(transcript)
    source = StubReadPhasingSource(
        phasing={variant: (variant,)},
        mutant_transcripts={(variant, transcript.id): stub_mt},
    )
    resolver = ReadPhaseResolver(source)

    variants = VariantCollection([variant])
    effects = variants.effects(phase_resolver=resolver)

    found_protein = None
    for e in effects:
        t = getattr(e, "transcript", None)
        if t is transcript and source.has_evidence(variant):
            found_protein = e.mutant_transcript.mutant_protein_sequence
            break
    assert found_protein == stub_mt.mutant_protein_sequence
