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

A *phase resolver* answers "are these two variants on the same
haplotype?" from some evidence source — DNA ``PS`` tags, RNA
read co-occurrence, an assembled contig from a long-read pipeline,
or anything else that can express variant co-occurrence.

Varcode defines two narrow Protocols here and plugs the resolvers into
:meth:`Variant.effects` / :meth:`VariantCollection.effects` via the
``phase_resolver`` kwarg. Varcode does **not** import or name any
specific upstream tool — implementations live in whatever package
produces the evidence (Isovar, long-read callers, custom assemblers,
test stubs).
"""

from typing import Optional, Protocol, Sequence, runtime_checkable


@runtime_checkable
class ReadPhasingSource(Protocol):
    """Reports per-variant read-level evidence and co-occurring partners.

    The minimum interface needed to answer ``in_cis(v1, v2)`` from
    read-level data — co-observation on the same supporting reads,
    same long-read fragment, same assembled contig, etc. Implementations
    decide how strong the evidence is; consumers only see a boolean
    membership question.
    """

    def has_evidence(self, variant) -> bool:
        """True if this source has any alt-supporting evidence for
        ``variant``."""
        ...

    def partners_in_cis(self, variant) -> Sequence:
        """Variants observed in cis with ``variant`` — i.e. on the
        same supporting reads / fragment / contig. May include
        germline SNPs, nearby somatic variants, or ``variant`` itself
        (implementations pick their convention). Empty sequence when
        no evidence covers ``variant``."""
        ...


@runtime_checkable
class MutantTranscriptSource(Protocol):
    """Reports an observed mutant transcript for ``(variant, transcript)``.

    Independent of phasing. A source can implement just
    :class:`ReadPhasingSource`, just :class:`MutantTranscriptSource`,
    or both. Consumers iterating effects use this channel to substitute
    the observed mutant protein (from RNA assembly, long-read calling,
    etc.) for the reference-inferred one.
    """

    def mutant_transcript(self, variant, transcript):
        """The :class:`~varcode.MutantTranscript` for
        ``(variant, transcript)``, or ``None`` when this source has
        no observed transcript for that pair."""
        ...


class ReadPhaseResolver:
    """Phase resolver backed by a :class:`ReadPhasingSource` (#269, #259).

    Two variants are cis if the source reports them as co-observed.
    That's direct molecular evidence — not a probabilistic call.

    Usage::

        # Any object satisfying ReadPhasingSource works. Common
        # implementation: an Isovar adapter shipped by openvax/isovar.
        resolver = ReadPhaseResolver(source)
        effects = variants.effects(phase_resolver=resolver)

    If ``source`` also satisfies :class:`MutantTranscriptSource`, the
    resolver routes :meth:`mutant_transcript` calls to it — so any
    effect whose ``(variant, transcript)`` is covered by the source
    gets its :attr:`~MutationEffect.mutant_transcript` populated with
    the observed mutant transcript.
    """

    #: Provenance string for downstream consumers that filter by
    #: phase source. Matches the pattern established by the
    #: :class:`Outcome.source` field.
    source = "read_phasing"

    def __init__(self, source: ReadPhasingSource):
        self._source = source

    def has_evidence(self, variant) -> bool:
        """Convenience passthrough to the wrapped source."""
        return self._source.has_evidence(variant)

    def mutant_transcript(self, variant, transcript):
        """Return the observed :class:`MutantTranscript` for
        ``(variant, transcript)``, or ``None`` when the wrapped source
        doesn't implement :class:`MutantTranscriptSource` or has no
        transcript for that pair."""
        get = getattr(self._source, "mutant_transcript", None)
        if get is None:
            return None
        return get(variant, transcript)

    def in_cis(self, v1, v2, transcript=None) -> Optional[bool]:
        """Return ``True`` if ``v1`` and ``v2`` are co-observed by the
        wrapped source, ``False`` if exactly one has evidence (so they
        are on distinct physical molecules — trans), ``None`` when
        neither has evidence.

        ``transcript`` is accepted for interface symmetry with
        :class:`VCFPhaseResolver.in_cis` but isn't consulted at the
        Protocol layer — isoform-specific sources are not yet a
        first-class concern. Reintroducible later as an optional
        Protocol extension.
        """
        v1_has = self._source.has_evidence(v1)
        v2_has = self._source.has_evidence(v2)
        if not v1_has and not v2_has:
            return None
        if v1_has:
            return v2 in self._source.partners_in_cis(v1)
        return v1 in self._source.partners_in_cis(v2)

    def phased_partners(self, variant, transcript=None) -> Sequence:
        """Variants co-observed with ``variant`` — i.e. the cis set.
        Empty when the source has no evidence for ``variant``."""
        if not self._source.has_evidence(variant):
            return ()
        return tuple(self._source.partners_in_cis(variant))


class VCFPhaseResolver:
    """Phase resolver backed by VCF ``GT`` + ``PS`` FORMAT fields.

    Reads the phase data that varcode's VCF loader already parses
    into :class:`~varcode.Genotype` (via #267): whether the
    ``GT`` delimiter was ``|`` (phased) or ``/`` (unphased), the
    ``PS`` phase-set identifier, and the per-haplotype allele indices
    in :attr:`Genotype.alleles`.

    Two variants are **cis** when they sit in the same phase set on
    the same haplotype slot, **trans** when they sit in the same
    phase set on different slots, and the resolver returns ``None``
    ("no evidence") for variants that aren't both phased, don't share
    a phase set, or lack called alleles.

    Compatible with any tool that writes standard-shaped VCF:
    WhatsHap, HapCUT2, DeepVariant, GATK HaplotypeCaller, long-read
    callers (PEPPER-DeepVariant, Clair3), population phasers
    (SHAPEIT5, Eagle2). varcode doesn't care which one wrote the
    file — it only reads ``GT`` and ``PS``.

    Multi-allelic sites are handled: varcode splits those rows into
    one :class:`~varcode.Variant` per ALT, each with an
    ``alt_allele_index`` preserved on the
    :class:`~varcode.VariantCollection` metadata. The resolver maps
    each variant to its GT-encoded index and asks "which haplotype
    slot carries this specific alt?".

    Single-sample by construction. Phase is per-sample; multi-sample
    VCFs need one resolver per sample.

    Currently supplies the cis/trans query but does **not** attach a
    :class:`~varcode.MutantTranscript` — DNA phasing alone doesn't
    produce an assembled contig. The natural next step is a
    ``HaplotypeEffect`` / multi-variant ``apply_variants_to_transcript``
    helper that, when two or more cis variants overlap the same
    transcript, builds a single joint :class:`MutantTranscript`
    applying all edits at once. That's a separate PR — this
    resolver already has the inputs it needs (``in_cis``) to drive
    the grouping.
    """

    #: Provenance tag, matching :attr:`ReadPhaseResolver.source`.
    source = "vcf_ps"

    def __init__(self, variant_collection, sample):
        self._collection = variant_collection
        self._sample = sample

    def _genotype(self, variant):
        try:
            return self._collection.genotype(variant, self._sample)
        except Exception:
            return None

    def _haplotype_slot(self, genotype, variant) -> Optional[int]:
        """Index into ``genotype.alleles`` where this variant's alt
        allele sits, or ``None`` when the sample doesn't carry the
        alt on any called haplotype.

        For a multi-allelic VCF row split into multiple Variants,
        each Variant has its own ``alt_allele_index``; this function
        honors that so ``GT=1|2`` correctly reports the first alt
        on slot 0 and the second alt on slot 1.
        """
        if genotype is None:
            return None
        alt_idx = self._collection._alt_index_for(variant)
        for i, allele in enumerate(genotype.alleles):
            if allele == alt_idx:
                return i
        return None

    @staticmethod
    def _is_homozygous_alt(genotype, alt_idx: int) -> bool:
        """True when every called haplotype carries ``alt_idx``. In
        that case phase is deterministic regardless of the ``|``/``/``
        delimiter — both copies carry the alt, so co-occurring
        homozygous-alt variants are trivially cis on both haplotypes.
        """
        if genotype is None:
            return False
        called = [a for a in genotype.alleles if a is not None]
        if not called:
            return False
        return all(a == alt_idx for a in called)

    def in_cis(self, v1, v2, transcript=None) -> Optional[bool]:
        """Return ``True`` if ``v1`` and ``v2`` are on the same
        haplotype in the same phase set, ``False`` if they're on
        different haplotypes in the same phase set, ``None`` when
        the phase relationship can't be determined (unphased GT,
        different phase sets, uncalled alleles).

        ``transcript`` is accepted for interface symmetry with
        :class:`ReadPhaseResolver.in_cis` but isn't consulted —
        DNA-level phase is isoform-agnostic.
        """
        g1 = self._genotype(v1)
        g2 = self._genotype(v2)
        if g1 is None or g2 is None:
            return None
        alt1 = self._collection._alt_index_for(v1)
        alt2 = self._collection._alt_index_for(v2)
        # Homozygous-alt on either side makes phase deterministic:
        # every haplotype carries this alt, so it's cis with anything
        # the OTHER variant sits on. Shortcut before requiring
        # phased=True.
        hom1 = self._is_homozygous_alt(g1, alt1)
        hom2 = self._is_homozygous_alt(g2, alt2)
        if hom1 and hom2:
            return True
        if hom1:
            # v1 is on every haplotype → cis with v2 iff v2 is called.
            return self._haplotype_slot(g2, v2) is not None
        if hom2:
            return self._haplotype_slot(g1, v1) is not None
        # Otherwise both must be phased and in the same phase set.
        if not g1.phased or not g2.phased:
            return None
        if g1.phase_set != g2.phase_set or g1.phase_set is None:
            return None
        s1 = self._haplotype_slot(g1, v1)
        s2 = self._haplotype_slot(g2, v2)
        if s1 is None or s2 is None:
            return None
        return s1 == s2

    def phased_partners(self, variant, transcript=None):
        """Variants in the collection that are cis with ``variant``
        under this resolver — i.e. sit in the same phase set on the
        same haplotype slot. Empty when ``variant`` isn't phased or
        has no called alt in the sample.
        """
        partners = []
        for other in self._collection:
            if other == variant:
                continue
            if self.in_cis(variant, other) is True:
                partners.append(other)
        return tuple(partners)


def apply_phase_resolver_to_effects(effects, phase_resolver):
    """Post-process an :class:`EffectCollection` (or any iterable of
    :class:`MutationEffect`) to attach observed
    :class:`MutantTranscript` objects when the resolver has evidence.

    Mutates each effect in place by setting
    ``effect.mutant_transcript``. Effects whose transcript isn't
    resolvable or whose ``(variant, transcript)`` has no observed
    transcript are left untouched — so this is safe to call on a mixed
    collection where only some variants have RNA evidence.
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
            # etc.) or populated from DNA-only inference. An observed
            # mutant transcript is higher-confidence evidence, so it
            # wins.
            e.mutant_transcript = mt
    return effects


def build_haplotype_effects(variant_collection, effects, phase_resolver):
    """Enrich ``effects`` with :class:`HaplotypeEffect` entries when
    ``phase_resolver`` groups two or more cis variants on the same
    transcript (#269).

    Additive — per-variant effects stay on the collection. Returns a
    flat list of new :class:`HaplotypeEffect` objects the caller can
    concatenate onto the existing :class:`EffectCollection`.

    Groups variants by transcript first, then uses
    ``phase_resolver.in_cis(v_i, v_j, transcript)`` to partition each
    transcript's variants into cis sets. Any cis set with ≥ 2
    resolvable edits becomes one HaplotypeEffect via
    :func:`apply_variants_to_transcript`.

    Edit conflicts (overlapping cDNA ranges) cause the joint build
    to return ``None`` — the group is skipped silently (the
    per-variant effects still describe each variant individually).
    """
    if phase_resolver is None or not hasattr(phase_resolver, "in_cis"):
        return []
    from .effects.effect_classes import HaplotypeEffect
    from .mutant_transcript import apply_variants_to_transcript

    # Group variants by transcript via existing per-variant effects —
    # each effect already knows its transcript, which avoids
    # re-running variant.transcripts for every variant.
    by_transcript = {}
    for e in effects:
        t = getattr(e, "transcript", None)
        v = getattr(e, "variant", None)
        if t is None or v is None:
            continue
        by_transcript.setdefault(t.id, (t, []))[1].append(v)
    # Dedup per transcript (a variant can appear in multiple effect
    # rows if it has multiple outcomes, e.g. SpliceOutcomeSet).
    haplotype_effects = []
    for transcript_id, (transcript, variants) in by_transcript.items():
        unique = []
        seen = set()
        for v in variants:
            key = (v.contig, v.start, v.end, v.ref, v.alt)
            if key not in seen:
                seen.add(key)
                unique.append(v)
        if len(unique) < 2:
            continue
        # Partition variants into cis groups using the resolver.
        # Greedy union-find: for each pair (i, j), if in_cis, merge
        # their components.
        parent = list(range(len(unique)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj

        for i in range(len(unique)):
            for j in range(i + 1, len(unique)):
                try:
                    cis = phase_resolver.in_cis(
                        unique[i], unique[j], transcript=transcript)
                except Exception:
                    cis = None
                if cis is True:
                    union(i, j)
        groups = {}
        for i in range(len(unique)):
            groups.setdefault(find(i), []).append(unique[i])
        for members in groups.values():
            if len(members) < 2:
                continue
            # Prefer a resolver-provided MutantTranscript (RNA assembly,
            # long-read, etc.) when available — it's the actual
            # observed molecule, not inferred from reference + edits.
            # Members of a cis group share a contig by definition, so
            # any member's contig covers the whole group.
            mt = None
            if hasattr(phase_resolver, "mutant_transcript"):
                for v in members:
                    mt = phase_resolver.mutant_transcript(v, transcript)
                    if mt is not None:
                        break
            if mt is None:
                mt = apply_variants_to_transcript(members, transcript)
            if mt is None:
                continue
            haplotype_effects.append(HaplotypeEffect(
                variants=members,
                transcript=transcript,
                mutant_transcript=mt,
                phase_source=getattr(phase_resolver, "source", None),
            ))
    return haplotype_effects
