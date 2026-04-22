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

"""Phase-resolver interfaces for cis/trans-aware effect prediction
(openvax/varcode#269).

A :class:`PhaseResolver` answers "are these two variants on the same
haplotype?" from some evidence source — DNA ``PS`` tags, RNA
read co-occurrence, or an assembled-contig tool like
`Isovar <https://github.com/openvax/isovar>`_.

Varcode defines the protocol here; the resolvers plug into
:meth:`Variant.effects` / :meth:`VariantCollection.effects` via the
``phase_resolver`` kwarg. Concrete DNA-based and read-based resolvers
can live here too or be shipped by downstream tools.

This module deliberately imports nothing from Isovar. An
:class:`IsovarAssemblyProvider` is a duck-typed interface — Isovar
(or any contig assembly tool) implements it; varcode consumes it.
"""

from typing import Optional, Protocol, Sequence, runtime_checkable


@runtime_checkable
class IsovarAssemblyProvider(Protocol):
    """Minimal interface for a source of per-locus RNA assemblies
    with phased variants (#269, #259).

    An "assembly" is a consensus contig built from RNA reads spanning
    a somatic variant. Variants co-occurring on the same contig are
    in cis — they came from the same physical molecule, so
    co-occurrence is direct evidence, not inference. The contig's
    translated protein IS the observed mutant product, not something
    inferred from reference + edits.

    Isovar is the reference implementation, but varcode intentionally
    does not import it. Any tool exposing this shape plugs in — a
    long-read haplotype caller, a custom assembler, or a mock in
    tests.

    Implementers need to be consistent about ``transcript`` — if the
    tool returns different contigs per isoform, the methods must key
    on ``(variant, transcript)``. If the tool is isoform-agnostic,
    all ``transcript`` arguments may be ignored (the provider just
    returns the same contig for any).
    """

    def has_contig(self, variant, transcript) -> bool:
        """True if this provider has an assembled contig covering
        ``variant`` on ``transcript``."""
        ...

    def variants_in_contig(self, variant, transcript) -> Sequence:
        """Variants observed on the same contig as ``variant``.

        May include germline SNPs, nearby somatic variants, or
        ``variant`` itself (implementations pick their convention).
        Empty sequence when no contig covers ``variant``.
        """
        ...

    def mutant_transcript(self, variant, transcript):
        """The :class:`~varcode.MutantTranscript` assembled from the
        contig. ``None`` when no contig covers
        ``(variant, transcript)``.
        """
        ...


class IsovarPhaseResolver:
    """Phase resolver backed by an :class:`IsovarAssemblyProvider`
    (#269, #259).

    Two variants are cis if they appear on the same assembled contig.
    That's direct molecular evidence — not a probabilistic call.

    Usage::

        isovar_results = run_isovar(bam, vcf)
        resolver = IsovarPhaseResolver(isovar_results)
        effects = variants.effects(phase_resolver=resolver)

    Any effect whose ``(variant, transcript)`` is covered by an
    assembled contig gets its :attr:`~MutationEffect.mutant_transcript`
    populated with the contig-derived :class:`MutantTranscript` —
    the protein attached to the effect is the protein actually
    observed in RNA, not one inferred from the reference.
    """

    #: Provenance string for downstream consumers that filter by
    #: phase source. Matches the pattern established by the
    #: :class:`Outcome.source` field.
    source = "isovar"

    def __init__(self, provider: IsovarAssemblyProvider):
        self.provider = provider

    def has_contig(self, variant, transcript) -> bool:
        """Convenience passthrough to the provider."""
        return self.provider.has_contig(variant, transcript)

    def mutant_transcript(self, variant, transcript):
        """Return the assembled :class:`MutantTranscript`, or ``None``
        when this provider has no contig for ``(variant, transcript)``.
        """
        if not self.provider.has_contig(variant, transcript):
            return None
        return self.provider.mutant_transcript(variant, transcript)

    def in_cis(self, v1, v2, transcript=None) -> Optional[bool]:
        """Return ``True`` if ``v1`` and ``v2`` appear on the same
        Isovar contig, ``False`` if they're each on a different
        contig (distinct physical molecules — trans), ``None`` when
        neither variant has a contig (no evidence).

        ``transcript`` is required because assemblies may be
        isoform-specific; pass ``None`` only if the provider is
        isoform-agnostic (the protocol allows this).
        """
        if transcript is None:
            return None
        v1_has = self.provider.has_contig(v1, transcript)
        v2_has = self.provider.has_contig(v2, transcript)
        if not v1_has and not v2_has:
            return None
        if v1_has:
            return v2 in self.provider.variants_in_contig(v1, transcript)
        return v1 in self.provider.variants_in_contig(v2, transcript)

    def phased_partners(self, variant, transcript) -> Sequence:
        """Variants observed on the same contig as ``variant`` on
        ``transcript`` — i.e. the cis set. Empty if no contig."""
        if not self.provider.has_contig(variant, transcript):
            return ()
        return tuple(self.provider.variants_in_contig(variant, transcript))


def apply_phase_resolver_to_effects(effects, phase_resolver):
    """Post-process an :class:`EffectCollection` (or any iterable of
    :class:`MutationEffect`) to attach contig-derived
    :class:`MutantTranscript` objects when the resolver has evidence.

    Mutates each effect in place by setting
    ``effect.mutant_transcript``. Effects whose transcript isn't
    resolvable or whose ``(variant, transcript)`` has no contig are
    left untouched — so this is safe to call on a mixed collection
    where only some variants have RNA evidence.
    """
    if phase_resolver is None:
        return effects
    if not hasattr(phase_resolver, "mutant_transcript"):
        return effects
    for e in effects:
        transcript = getattr(e, "transcript", None)
        variant = getattr(e, "variant", None)
        if variant is None or transcript is None:
            continue
        mt = phase_resolver.mutant_transcript(variant, transcript)
        if mt is not None:
            # Intentional mutation: the effect's mutant_transcript
            # slot was either None (point variants, cryptic stubs,
            # etc.) or populated from DNA-only inference. Isovar's
            # assembly is higher-confidence evidence, so it wins.
            e.mutant_transcript = mt
    return effects
