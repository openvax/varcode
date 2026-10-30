"""Reference open reading frames using Varcode's transcript-model engine.

The optional orfdata catalogue supplies reference observations, never patient
translation evidence. Canonical ``effects()`` results are unchanged.
"""

from dataclasses import asdict, dataclass, field, fields, is_dataclass
from hashlib import sha256
from typing import TYPE_CHECKING, Mapping, Optional, Tuple

from .effect_hypotheses import EffectHypothesis, RealizedTranscriptProduct
from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import codon_table_for_transcript, translate_sequence
from .effects.effect_classes import MutationEffect, StartLoss, Unresolved
from .genomic_layout import SequenceUnavailable, UnsupportedLayoutEdit
from .mutant_transcript import MutantTranscript
from .nucleotides import reverse_complement
from .transcript_model import predict_transcript_model_effect

if TYPE_CHECKING:
    from orfdata import InitiationContext, ORFRecord, ReferenceORF


def load_orf_catalog(*args, **kwargs):
    """Load the shared orfdata reference; install ``varcode[orf]`` to use it."""
    try:
        from orfdata import load_catalog as load
    except ImportError as error:
        raise ImportError("Install varcode[orf] for the shared reference catalogue") from error
    return load(*args, **kwargs)


def _contig(value, reference_name=None):
    from .contigs import normalize_contig_name
    return normalize_contig_name(str(value), reference_name=reference_name)


def _positions(transcript):
    exons = sorted(transcript.exons, key=lambda e: e.start,
                   reverse=transcript.on_backward_strand)
    return tuple(p for e in exons for p in (
        range(e.end, e.start - 1, -1) if transcript.on_backward_strand
        else range(e.start, e.end + 1)))


def _model_positions(record):
    positions = tuple(p for a, b in record.coding_blocks for p in range(a, b + 1))
    return positions[::-1] if record.strand == "-" else positions


def _orf_key(reference):
    transcript = reference.transcript
    return (reference.context.key, transcript.id,
            getattr(transcript, "versioned_id", transcript.id), reference.source,
            reference.orf_id, tuple(reference.coding_blocks), reference.translation_table,
            reference.protein_sequence, reference.start_codon, reference.kind,
            reference.start_complete, reference.stop_complete, reference.coding_origin,
            reference.cdna_interval)


class _ORFTranscript:
    """Private coding view used only by Varcode's transcript-model engine.

    Public callers receive a ReferenceORF, whose parent transcript retains its
    annotation identity and biotype. This adapter isolates the legacy engine's
    requirement that the selected coding frame look like a coding transcript.
    """

    def __init__(self, reference):
        self.reference_orf = reference
        self.record = reference.record
        transcript = reference.transcript
        self.reference_transcript = transcript
        self.positions = reference.genomic_positions
        self.sequence = reference.transcript_sequence
        start, end = reference.cdna_interval
        wanted = reference.coding_positions
        self.protein_sequence = reference.protein_sequence
        self.start_codon_spliced_offsets = tuple(range(start, start + 3)) if reference.start_complete else ()
        self.stop_codon_spliced_offsets = tuple(range(end - 3, end)) if reference.stop_complete else ()
        self.start_codon_positions = wanted[:3] if reference.start_complete else ()
        self.stop_codon_positions = wanted[-3:] if reference.stop_complete else ()
        self.coding_sequence = reference.coding_sequence
        self.coding_sequence_position_ranges = reference.coding_blocks
        self.five_prime_utr_sequence = self.sequence[:start]
        self.three_prime_utr_sequence = self.sequence[end:]
        self.is_protein_coding = True
        self.complete = reference.start_complete is True and reference.stop_complete is True
        self.contains_start_codon = reference.start_complete
        self.contains_stop_codon = reference.stop_complete
        self.translation_table = reference.translation_table
        self.start_codon = reference.start_codon
        # Distinct coding models on one RNA must survive candidate deduplication.
        self.id = "%s:orf:%s" % (transcript.id, sha256(repr(_orf_key(reference)).encode()).hexdigest())
        self.name = "%s:%s" % (transcript.id, reference.orf_id)
        if transcript.genome is None:
            from .orf_serialization import _Feature, _ReferenceGenome
            context = reference.context
            self.genome = _ReferenceGenome(dict(
                reference_name=context.genome_build, species=context.species, taxon_id=context.taxon_id,
                annotation_name=context.annotation_name, annotation_version=context.annotation_version,
                assembly_accession=context.assembly_accession, strain=context.strain,
                strain_taxon_id=context.strain_taxon_id))
            self.exons = tuple(_Feature(a, b) for a, b in transcript.exons)
            self.start = min(self.positions) if self.positions else 0
            self.end = max(self.positions) if self.positions else 0
            self.on_backward_strand = transcript.strand == "-"
            self.on_forward_strand = transcript.strand == "+"

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "reference_transcript"), name)

    def spliced_offset(self, position):
        return self.positions.index(position)

    def to_dict(self):
        raise ValueError("A private ORF coding view cannot serialize as its parent transcript; "
                         "use ORFAnnotation.to_dict()")

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return isinstance(other, _ORFTranscript) and self.id == other.id


def _bind(record, transcript, genome_build, *, reference_context=None):
    if genome_build != record.genome_build:
        raise ValueError("ORF assembly mismatch")
    from orfdata import ReferenceORF
    return ReferenceORF(record, transcript, context=reference_context)


def _reference_provider(view, fallback=None):
    offsets = {p: i for i, p in enumerate(view.positions)}
    build = view.reference_orf.genome_build

    def sequence(contig, start, end):
        if _contig(contig, build) == _contig(view.contig, build) and all(p in offsets for p in range(start, end + 1)):
            genomic = "".join(view.sequence[offsets[p]] for p in range(start, end + 1))
            return genomic.translate(str.maketrans("ACGT", "TGCA")) if view.on_backward_strand else genomic
        if fallback is not None:
            return fallback(contig, start, end)
        raise SequenceUnavailable("Intronic/genomic reference sequence unavailable")
    return sequence


def _plain(value):
    """Lossless JSON values for hypotheses without effect/transcript cycles."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if hasattr(value, "contig") and hasattr(value, "ref") and hasattr(value, "alt"):
        return {k: _plain(getattr(value, k, None)) for k in
                ("contig", "start", "end", "ref", "alt", "reference_name",
                 "original_start", "original_ref", "original_alt", "sv_type",
                 "affected_start", "affected_end", "mate_contig", "mate_start",
                 "mate_orientation", "ci_start", "ci_end", "alt_assembly", "info")}
    if is_dataclass(value):
        return {f.name: _plain(getattr(value, f.name)) for f in fields(value)}
    raise TypeError("Unsupported annotation evidence value: %s" % type(value).__name__)


def _context(sequence, start):
    from orfdata import InitiationContext, annotate_initiation_context

    if start is None or sequence is None or start + 3 > len(sequence):
        return InitiationContext("", "", "", unavailable_reason="reference_start_not_retained")
    return annotate_initiation_context(sequence, start)


def _peptide_overlay(reference, protein, origins, start):
    """Track native residues by codon ancestry, not coincidental peptide matches."""
    native = reference.coding_positions
    retained = []
    if protein is not None and start is not None:
        origin = reference.coding_origin
        by_codon = {} if origin is None else {
            tuple(native[origin + 3 * i:origin + 3 * i + 3]): i + 1
            for i in range(len(reference.protein_sequence))
            if 0 <= origin + 3 * i and origin + 3 * i + 3 <= len(native)}
        for index, aa in enumerate(protein):
            codon = tuple(origins[start + index * 3:start + index * 3 + 3])
            ref = by_codon.get(codon)
            if ref is not None and reference.protein_sequence[ref - 1] == aa:
                retained.append((ref, index + 1))
    mapping = {}
    for ref, alt in retained:
        mapping.setdefault(ref, []).append(alt)
    observations = []
    for peptide in reference.evidence.peptides:
        intervals = []
        for a, b in peptide.protein_intervals:
            for start in mapping.get(a, ()):
                if all(start + i - a in mapping.get(i, ()) for i in range(a, b + 1)):
                    intervals.append((start, start + b - a))
        observations.append(dict(reference_observation=asdict(peptide),
                                 retained_mutant_intervals=intervals,
                                 reference_peptide_retained=bool(intervals) if protein is not None else None))
    return dict(retained_native_residues=retained, peptide_observations=observations,
                reference_ms_support_on_retained_sequence=any(
                    o["reference_peptide_retained"] and p.is_supporting
                    for o, p in zip(observations, reference.evidence.peptides)),
                mutant_translation_evidence="not_assessed",
                mutant_mhc_presentation_evidence="not_assessed",
                mutant_hla_presentation_evidence="not_assessed" if reference.context.taxon_id == 9606 else "not_applicable")


def _product(view, product):
    origins = tuple(product.evidence.get("genomic_positions", ()))
    target = view.start_codon_positions[0]
    indices = [i for i, p in enumerate(origins) if p == target]
    start = indices[0] if len(indices) == 1 else None
    sequence = product.cdna_sequence
    protein = product.protein_sequence
    context = _context(sequence, start)
    table = codon_table_for_transcript(view)
    stop = None
    if start is not None and sequence is not None:
        stop = next((i for i in range(start, len(sequence) - 2, 3)
                     if sequence[i:i + 3] in table.stop_codons), None)
    return ORFProduct(cdna_sequence=sequence, protein_sequence=protein,
                start_codon_present=product.start_codon_present,
                termination="observed_stop" if stop is not None else "sequence_boundary" if sequence is not None else "unavailable",
                stop_cdna_offset=stop, initiation_context=context,
                junction_signature=product.junction_signature,
                evidence=product.evidence,
                reference_evidence=_peptide_overlay(view.reference_orf, protein, origins, start))


def _effect(effect):
    return {"type": type(effect).__name__, **{
        k: _plain(getattr(effect, k)) for k in (
            "aa_mutation_start_offset", "aa_mutation_end_offset", "aa_ref", "aa_alt",
            "reason", "mechanism", "max_hypotheses") if hasattr(effect, k)}}


@dataclass(frozen=True)
class ORFProduct:
    """A realized ORF's sequences, initiation context and reference evidence.

    Unavailable sequences remain None. Reference peptide retention is not
    evidence that this product was translated or presented by MHC.
    """

    cdna_sequence: Optional[str]
    protein_sequence: Optional[str]
    start_codon_present: Optional[bool]
    termination: str
    stop_cdna_offset: Optional[int]
    initiation_context: "InitiationContext"
    junction_signature: tuple = ()
    evidence: Mapping = field(default_factory=dict)
    reference_evidence: Mapping = field(default_factory=dict)

    def to_dict(self):
        """Export the established product schema without changing native fields."""
        context = self.initiation_context
        if context.unavailable_reason:
            initiation = dict(status="unavailable", reason=context.unavailable_reason)
        else:
            initiation = dict(
                upstream=context.upstream, start_codon=context.start_codon,
                downstream=context.downstream, start_offset=context.start_offset,
                five_prime_complete=context.five_prime_complete,
                interpretation="sequence_prior_only; initiation_not_measured")
        return dict(
            cdna_sequence=self.cdna_sequence, protein_sequence=self.protein_sequence,
            start_codon_present=self.start_codon_present, termination=self.termination,
            stop_cdna_offset=self.stop_cdna_offset, initiation_context=initiation,
            junction_signature=_plain(self.junction_signature),
            evidence=_plain(self.evidence), reference_evidence=_plain(self.reference_evidence))


@dataclass(frozen=True)
class ORFOutcome:
    """One consequence and its baseline/mutant products, with uncertainty intact."""

    effect: MutationEffect
    baseline: ORFProduct
    mutant: Optional[ORFProduct]
    hypothesis: Optional[EffectHypothesis] = None
    probability: Optional[float] = None
    phase_mass: Optional[float] = None
    source: Optional[str] = None
    unresolved_reason: Optional[str] = None

    def to_dict(self):
        result = dict(effect=_effect(self.effect), baseline=self.baseline.to_dict(),
                      mutant=self.mutant.to_dict() if self.mutant is not None else None)
        if self.hypothesis is not None:
            result.update(hypothesis=_plain(self.hypothesis),
                          probability=self.probability, phase_mass=self.phase_mass)
        else:
            result["unresolved_reason"] = self.unresolved_reason
            if self.mutant is not None:
                result["source"] = self.source
        return result


@dataclass(frozen=True)
class ORFAnnotation:
    """Per-ORF consequences available through ``outcomes`` as Python objects.

    A failed reference binding leaves ``reference`` and ``effect`` as None,
    with an explicit ``unresolved_reason``. ``to_dict()`` is for export only.
    """

    record: Optional["ORFRecord"]
    effect: Optional[MutationEffect]
    reference: Optional["ReferenceORF"] = None
    variants: tuple = ()
    outcomes: Tuple[ORFOutcome, ...] = ()
    catalogue_metadata: Optional[Mapping] = None
    variant_locations: tuple = ()
    unresolved_reason: Optional[str] = None
    unavailable_operations: Mapping = field(default_factory=dict)

    @property
    def orf_id(self):
        return self.reference.orf_id if self.reference is not None else self.record.orf_id

    @property
    def data(self):
        """Compatibility export; use attributes for in-process access."""
        return self.to_dict()

    def to_dict(self):
        result = dict(
            schema="varcode.orf_annotation.v1", model=None if self.record is None else self.record.to_dict(),
            catalogue_metadata=self.catalogue_metadata, variants=_plain(self.variants),
            effect=_effect(self.effect) if self.effect is not None else dict(
                type="Unresolved", mechanism="uorf_model_unavailable", reason=self.unresolved_reason),
            outcomes=[outcome.to_dict() for outcome in self.outcomes],
            discovery=dict(status="reference_models_only", new_starts="not_assessed"))
        if self.unavailable_operations:
            result["unavailable_operations"] = dict(self.unavailable_operations)
        if self.reference is not None:
            transcript = self.reference.transcript
            result.update(
                reference=self.reference.to_dict(),
                resolved_transcript=dict(id=transcript.id,
                    version=getattr(transcript, "version", None),
                    release=getattr(transcript.genome, "release", None),
                    genome_build=self.reference.context.genome_build),
                variant_locations=_plain(self.variant_locations),
                downstream_cds_regulation="not_assessed")
        return _plain(result)


@dataclass(frozen=True)
class ORFProductAnnotation:
    """Consequence of a mapped RNA product relative to a reference ORF.

    Attribute access retains native objects. Use ``to_dict()`` for JSON export.
    """

    reference: "ReferenceORF"
    status: str
    effect: Optional[MutationEffect] = None
    product: Optional[ORFProduct] = None
    reason: Optional[str] = None
    cdna_sequence: Optional[str] = None
    evidence: Mapping = field(default_factory=dict)

    def to_dict(self):
        result = dict(schema="varcode.orf_product.v1", model=(None if self.reference.record is None
                                                            else self.reference.record.to_dict()),
                      reference=self.reference.to_dict(), status=self.status)
        if self.product is None:
            result.update(reason=self.reason, cdna_sequence=self.cdna_sequence,
                          evidence=_plain(self.evidence))
        else:
            result.update(effect=_effect(self.effect), product=self.product.to_dict())
        return result


def _finish_annotation(annotation):
    from .version import __version__

    for effect in (annotation.effect, *(outcome.effect for outcome in annotation.outcomes)):
        if effect is not None:
            effect.annotator = "transcript_model"
            effect.annotator_version = __version__
    if annotation.effect is not None:
        annotation.effect.orf_annotation = annotation
    return annotation


def _boundary_operations(reference):
    unavailable = {}
    for boundary, operation in (("start", "start_loss"), ("stop", "stop_loss")):
        state = getattr(reference, boundary + "_complete")
        if state is not True:
            unavailable[operation] = "reference_%s_%s" % (boundary, "unknown" if state is None else "incomplete")
    if unavailable:
        unavailable["complete_mutant_protein"] = "reference_coding_fragment"
    return unavailable


def _baseline_product(view):
    reference = view.reference_orf
    origin = reference.coding_origin
    return ORFProduct(
        cdna_sequence=view.sequence, protein_sequence=reference.protein_sequence,
        start_codon_present=reference.start_complete,
        termination=("observed_stop" if reference.stop_complete else
                     "sequence_boundary" if reference.stop_complete is False else "unknown_boundary"),
        stop_cdna_offset=reference.end - 3 if reference.stop_complete else None,
        initiation_context=_context(view.sequence, reference.start if reference.start_complete else None),
        evidence=dict(genomic_positions=view.positions, coding_interval=reference.cdna_interval,
                      coding_phase=reference.coding_phase),
        reference_evidence=_peptide_overlay(reference, reference.protein_sequence, view.positions,
                                            None if origin is None else reference.start + origin))


def _unavailable_annotation(view, variants, reason, metadata=None):
    effect = Unresolved(variants[0], view, mechanism="orf_consequence_unavailable", reason=reason)
    return _finish_annotation(ORFAnnotation(
        record=view.record, reference=view.reference_orf, effect=effect, variants=variants,
        catalogue_metadata=metadata, unresolved_reason=reason,
        unavailable_operations={**_boundary_operations(view.reference_orf), "variant_consequence": reason},
        variant_locations=_variant_locations(variants, view.reference_orf),
        outcomes=(ORFOutcome(effect, _baseline_product(view), None, unresolved_reason=reason),)))


def _unsupported_reference(view):
    """Check the ordinary coding engine's assumptions against the actual source model."""
    reference = view.reference_orf
    if reference.wraps_junction:
        return "circular_orf_consequence_unavailable"
    if not reference.coding_positions:
        return "orf_genomic_mapping_unavailable"
    if reference.coding_origin is None:
        return "orf_reading_frame_unknown"
    if reference.subcodon_introns:
        return "orf_programmed_frame_consequence_unavailable"
    table = codon_table_for_transcript(view)
    phase = reference.coding_phase
    coding = reference.coding_sequence[phase:]
    coding = coding[:len(coding) // 3 * 3]
    if set(coding) - set("ACGT"):
        return "orf_reference_codon_ambiguous"
    translated = translate_sequence(coding, codon_table=table, to_stop=True)
    if reference.start_complete and translated:
        translated = "M" + translated[1:]
    skip = 1 if reference.coding_origin < 0 else 0
    if translated != reference.protein_sequence[skip:]:
        return "orf_reference_translation_requires_recoding"
    return None


def _annotate_partial(view, variants, *, germline_variants=(), metadata=None):
    """Calculate local products within known coding bounds, without inventing missing ends."""
    from .mutant_transcript import _resolve_variant_edit

    reference = view.reference_orf
    if germline_variants:
        return _unavailable_annotation(view, variants, "partial_orf_germline_phasing_unavailable", metadata)
    first = reference.start + reference.coding_phase
    end = reference.end - (reference.end - first) % 3
    if reference.stop_complete:
        end -= 3
    edits = []
    for variant in variants:
        if getattr(variant, "is_structural", False):
            return _unavailable_annotation(view, variants, "partial_orf_structural_product_unavailable", metadata)
        resolved = _resolve_variant_edit(variant, view, view.sequence)
        if resolved is None:
            return _unavailable_annotation(view, variants, "partial_orf_splice_edit_unavailable", metadata)
        edit, _ = resolved
        if not first <= edit.cdna_start < end or not edit.cdna_start <= edit.cdna_end <= end:
            return _unavailable_annotation(view, variants, "partial_orf_boundary_or_partial_codon_edit", metadata)
        edits.append(edit)
    edits.sort(key=lambda e: e.cdna_start)
    if any(a.cdna_end > b.cdna_start or a.cdna_start == b.cdna_start for a, b in zip(edits, edits[1:])):
        return _unavailable_annotation(view, variants, "overlapping_partial_orf_edits", metadata)
    sequence, positions = view.sequence, list(view.positions)
    delta = 0
    for edit in reversed(edits):
        a, b = edit.cdna_start, edit.cdna_end
        sequence = sequence[:a] + edit.alt_bases + sequence[b:]
        if len(edit.alt_bases) != b - a:
            positions[a:b] = [None] * len(edit.alt_bases)
        delta += len(edit.alt_bases) - (b - a)
    bound_end = reference.end + delta
    table = codon_table_for_transcript(view)
    coding = sequence[first:bound_end]
    coding = coding[:len(coding) // 3 * 3]
    if set(coding) - set("ACGT"):
        return _unavailable_annotation(view, variants, "partial_orf_mutant_codon_ambiguous", metadata)
    protein = translate_sequence(coding, codon_table=table, to_stop=True)
    start_present = reference.start_complete
    if reference.start_complete:
        start_present = coding[:3] in table.start_codons
        protein = ("M" + protein[1:]) if start_present and protein else ""
    skip = 1 if reference.coding_origin < 0 else 0
    protein = reference.protein_sequence[:skip] + protein
    stop = next((first + i for i in range(0, len(coding), 3) if coding[i:i + 3] in table.stop_codons), None)
    mt = MutantTranscript(view, edits=tuple(edits), cdna_sequence=sequence,
                          mutant_protein_sequence=protein, annotator_name="transcript_model",
                          evidence=dict(coding_interval=(reference.start, bound_end),
                                        coding_phase=reference.coding_phase, product_extent="coding_fragment"))
    effect = classify_from_protein_diff(
        variants[0], view, reference.protein_sequence, protein, length_delta=delta, mutant_transcript=mt,
        frame_change_aa_offset=max(0, (edits[0].cdna_start - reference.start - reference.coding_origin) // 3),
        start_complete=reference.start_complete is True,
        stop_complete=reference.stop_complete is True,
        coding_start=reference.start + reference.coding_origin, stop_observed=stop is not None)
    baseline = _baseline_product(view)
    product = ORFProduct(
        cdna_sequence=sequence, protein_sequence=protein, start_codon_present=start_present,
        termination=("observed_stop" if stop is not None else "sequence_boundary"
                     if reference.stop_complete is not None else "unknown_boundary"),
        stop_cdna_offset=stop, initiation_context=_context(sequence, reference.start if start_present else None),
        evidence=dict(genomic_positions=tuple(positions), coding_interval=(reference.start, bound_end),
                      coding_phase=reference.coding_phase, product_extent="coding_fragment"),
        reference_evidence=_peptide_overlay(reference, protein, positions,
                                            reference.start + reference.coding_origin))
    return _finish_annotation(ORFAnnotation(
        record=reference.record, reference=reference, effect=effect, variants=variants,
        outcomes=(ORFOutcome(effect, baseline, product),), catalogue_metadata=metadata,
        variant_locations=_variant_locations(variants, reference),
        unavailable_operations=_boundary_operations(reference)))


def _annotate_orf_variants(variants, record, transcript, *, genome_build,
                           germline_variants=(), phase_resolver=None,
                           sequence_provider=None, max_hypotheses=64,
                           homozygous_germline=(), metadata=None, reference=None,
                           reference_context=None):
    """Annotate one known-cis group in a validated upstream frame.

    Literal REF mismatches raise. Unavailable intronic sequence and unsupported
    events remain unresolved. A caller-supplied transcript must have the exact
    model's spliced path and protein; canonical labels are not replaced.
    """
    variants = tuple(variants)
    if not variants:
        raise ValueError("At least one variant is required")
    if reference is not None and reference_context is not None and reference.context.key != reference_context.key:
        raise ValueError("Bound ORF and explicit reference context disagree")
    view = _ORFTranscript(_bind(record, transcript, genome_build,
                              reference_context=reference_context) if reference is None else reference)
    germline_variants = tuple(germline_variants)
    for variant in (*variants, *germline_variants):
        if variant.reference_name != genome_build:
            raise ValueError("Variant assembly mismatch")
        if view.contig is not None and _contig(variant.contig, genome_build) != _contig(view.contig, genome_build):
            raise ValueError("Variant contig mismatch")
        view.reference_orf.context.validate_genome(variant.genome, include_annotation=False)
    if not view.positions:
        return _unavailable_annotation(view, variants, "orf_genomic_mapping_unavailable", metadata)
    if sequence_provider is None:
        from .transcript_model import _provider_for_genome
        sequence_provider = _provider_for_genome(variants[0].genome, transcript)
    provider = _reference_provider(view, sequence_provider)
    try:
        # Validate original genomic REF, including VCF padding anchors. Layout
        # tier 0 deliberately permits unknown intronic REF; this API does not.
        from .errors import ReferenceMismatchError
        for variant in (*variants, *germline_variants):
            if not getattr(variant, "is_structural", False) and variant.original_ref:
                expected = provider(variant.contig, variant.original_start, variant.original_start + len(variant.original_ref) - 1)
                if expected != variant.original_ref:
                    observed = variant.original_ref
                    if view.on_backward_strand:
                        expected = reverse_complement(expected)
                        observed = reverse_complement(observed)
                    raise ReferenceMismatchError(
                        variant=variant, transcript=view.reference_transcript,
                        expected_ref=expected, observed_ref=observed,
                        genome_start=variant.original_start,
                        genome_end=variant.original_start + len(variant.original_ref) - 1)
        reason = _unsupported_reference(view)
        if reason is not None:
            return _unavailable_annotation(view, variants, reason, metadata)
        if not view.complete:
            return _annotate_partial(view, variants, germline_variants=germline_variants, metadata=metadata)
        effect = predict_transcript_model_effect(
            variants, view, germline_variants=germline_variants,
            phase_resolver=phase_resolver,
            sequence_provider=provider,
            max_hypotheses=max_hypotheses, homozygous_germline=homozygous_germline)
    except (SequenceUnavailable, UnsupportedLayoutEdit) as error:
        effect = Unresolved(variants[0], view, mechanism="uorf_sequence_unavailable", reason=str(error))
    outcomes = []
    variant_locations = _variant_locations(variants, view.reference_orf)
    for candidate in getattr(effect, "realized_candidates", ()):
        for outcome in candidate.outcomes:
            outcomes.append(ORFOutcome(
                effect=outcome.effect, hypothesis=outcome.hypothesis,
                baseline=_product(view, outcome.baseline), mutant=_product(view, outcome.mutant),
                probability=candidate.probability, phase_mass=candidate.phase_mass))
    if not outcomes and getattr(variants[0], "is_structural", False):
        baseline = RealizedTranscriptProduct(
            view, view.sequence, view.protein_sequence, start_codon_present=True,
            evidence={"genomic_positions": view.positions})
        for candidate in getattr(effect, "candidates", ()):
            alternative = candidate.effect
            mt = getattr(alternative, "mutant_transcript", None)
            if mt is None or mt.cdna_sequence is None:
                outcomes.append(ORFOutcome(effect=alternative,
                    baseline=_product(view, baseline), mutant=None,
                    unresolved_reason="complete_structural_product_unavailable"))
                continue
            origins = []
            for segment in mt.reference_segments or ():
                source = segment.source
                if not hasattr(source, "exons"):
                    origins.extend([None] * segment.length)
                    continue
                positions = _positions(source)[segment.start:segment.end]
                strand = "-" if source.on_backward_strand else "+"
                if segment.strand == "-":
                    positions = positions[::-1]
                    strand = "+" if strand == "-" else "-"
                origins.extend((source.contig, position, strand) for position in positions)
            if mt.edits or len(origins) != len(mt.cdna_sequence):
                outcomes.append(ORFOutcome(effect=alternative,
                    baseline=_product(view, baseline), mutant=None,
                    unresolved_reason="structural_base_origins_unavailable"))
                continue
            product = RealizedTranscriptProduct(
                view, mt.cdna_sequence, mt.mutant_protein_sequence,
                evidence={**(mt.evidence or {}), "base_origins": origins})
            observed = annotate_orf_rna_product(view.reference_orf, product)
            outcomes.append(ORFOutcome(effect=alternative,
                baseline=_product(view, baseline), mutant=observed.product,
                source=getattr(candidate, "source", None),
                unresolved_reason=observed.reason))
    return _finish_annotation(ORFAnnotation(
        record=record, effect=effect, reference=view.reference_orf,
        variants=variants, outcomes=tuple(outcomes), catalogue_metadata=metadata,
        variant_locations=tuple(variant_locations)))


def _variant_locations(variants, reference):
    variant_locations = []
    record = reference.record or reference
    genome_build = reference.context.genome_build
    for variant in variants:
        context_snvs = []
        context = record.initiation_context
        positions = set(_model_positions(record))
        if context is not None:
            positions.update(context.genomic_positions)
        locations = []
        for position in sorted(p for p in positions if variant.start <= p <= variant.end):
            coding = record.map_genomic_position(position, genome_build=genome_build)
            initiation = record.map_initiation_position(position, genome_build=genome_build)
            locations.append(dict(position=position, coding=coding, initiation=initiation))
            if initiation and len(variant.ref) == len(variant.alt) == 1 and variant.ref != variant.alt:
                if initiation["genomic_ref"] == variant.ref:
                    context_snvs.append(record.annotate_initiation_snv(
                        position, variant.ref, variant.alt, genome_build=genome_build))
                else:
                    context_snvs.append(dict(position=position,
                        status="source_context_differs_from_resolved_transcript"))
        variant_locations.append(dict(variant=variant, locations=locations,
            initiation_snvs=context_snvs,
            insertion_after_position=variant.start if getattr(variant, "is_insertion", False) else None))
    return tuple(variant_locations)


def predict_orf_effects(variants, *, catalog=None, transcripts=None, **kwargs):
    """Return additive annotations for overlapping catalogue models.

    ``variants`` is one variant or an explicitly known-cis group. Transcript
    lookup uses the variant's genome unless an iterable is supplied. Failed
    model validation is exported as unresolved, never silently discarded.
    An empty tuple means no reference model overlapped this query, not that
    the variants cannot create or affect translation upstream of a CDS.
    """
    variants = (variants,) if hasattr(variants, "contig") else tuple(variants)
    if not variants:
        return ()
    catalog = load_orf_catalog() if catalog is None else catalog
    genome = variants[0].genome
    build = variants[0].reference_name
    if build != catalog.genome_build:
        raise ValueError("Catalogue assembly mismatch")
    by_id = None
    if transcripts is not None:
        by_id = {}
        for transcript in transcripts:
            key = str(transcript.id).split(".")[0]
            if key in by_id:
                raise ValueError("Ambiguous supplied transcript identifier: %s" % key)
            by_id[key] = transcript
    results = []
    for record in catalog.records:
        # Include splice sites inside the ORF envelope and initiation context.
        positions = (record.coding_blocks[0][0], record.coding_blocks[-1][1])
        context = record.initiation_context
        if context is not None:
            positions += tuple(context.genomic_positions)
        if not any(_contig(v.contig, build) == _contig(record.contig, build)
                   and v.start <= max(positions) + 6 and v.end >= min(positions) - 6
                   for v in variants):
            continue
        try:
            key = record.transcript_id.split(".")[0]
            transcript = genome.transcript_by_id(key) if by_id is None else by_id[key]
            annotation = _annotate_orf_variants(
                variants, record, transcript, genome_build=build,
                metadata=catalog.metadata, **kwargs)
        except (KeyError, ValueError) as error:
            from .errors import ReferenceMismatchError
            if isinstance(error, ReferenceMismatchError):
                raise
            # There is no validated transcript to attach a MutationEffect to.
            # Keep the failure at the ORF-annotation boundary instead of
            # manufacturing a transcript or passing None to a transcript effect.
            annotation = ORFAnnotation(
                record=record, effect=None, variants=variants,
                catalogue_metadata=catalog.metadata, unresolved_reason=str(error))
        results.append(annotation)
    return tuple(results)


def annotate_orf_rna_product(reference, product):
    """Annotate an observed/reconstructed RNA product with mapped base ancestry.

    ``product`` is a :class:`RealizedTranscriptProduct`. Its evidence must carry
    ``base_origins``: one ``(contig, position, strand)`` or None per cDNA base.
    These placements must come from the reconstruction, not peptide similarity.
    RNA support in ``product.evidence`` remains separate from reference MS.
    Only a uniquely retained complete reference start anchors translation.

    Returns an :class:`ORFProductAnnotation` with native ``reference``,
    ``effect`` and ``product`` attributes. Unresolved results have no product
    and carry a ``reason``. Use ``to_dict()`` for JSON export.
    """
    view = _ORFTranscript(reference)
    record = reference
    product_genome = getattr(product.transcript, "genome", None)
    if product_genome is None:
        raise ValueError("RNA product reference Genome is required")
    from orfdata import ReferenceContext
    # The bound reference is the explicit context for custom Genome metadata
    # gaps. Existing metadata must agree; aliases use the shared taxonomy map.
    ReferenceContext.from_genome(product_genome, context=reference.context)
    sequence = product.cdna_sequence
    origins = product.evidence.get("base_origins", ())
    if sequence is None or len(sequence) != len(origins):
        raise ValueError("A sequence and one base origin per nucleotide are required")
    reason = _unsupported_reference(view) or (None if view.complete else "partial_orf_rna_product_unavailable")
    if reason is not None:
        return ORFProductAnnotation(reference=reference, status="unresolved", reason=reason,
                                    cdna_sequence=sequence, evidence=product.evidence)
    positions = tuple(o[1] if o is not None and _contig(o[0], reference.genome_build) == _contig(record.contig, reference.genome_build)
                      and o[2] == record.strand else None for o in origins)
    starts = [i for i in range(len(positions) - 2)
              if positions[i:i + 3] == tuple(view.start_codon_positions)]
    if len(starts) != 1:
        return ORFProductAnnotation(reference=reference,
                    status="unresolved", reason="reference_start_missing_or_duplicated",
                    cdna_sequence=sequence, evidence=product.evidence)
    start = starts[0]
    table = codon_table_for_transcript(view)
    start_present = sequence[start:start + 3] in table.start_codons
    coding = sequence[start:]
    protein = translate_sequence(coding[:len(coding) // 3 * 3], codon_table=table, to_stop=True) if start_present else ""
    if protein:
        protein = "M" + protein[1:]
    from copy import copy
    product_view = copy(view)
    product_view.start_codon_spliced_offsets = tuple(range(start, start + 3))
    mt = MutantTranscript(reference_transcript=product_view, cdna_sequence=sequence,
                          mutant_protein_sequence=protein,
                          evidence={**product.evidence, "start_codon_offset": start})
    # The shared classifier compares protein outcomes; DNA linkage and source
    # placements remain in the RNA evidence supplied by the caller.
    effect = (classify_from_protein_diff(None, product_view, record.protein_sequence, protein,
                                       mutant_transcript=mt) if start_present else StartLoss(None, view, mutant_transcript=mt))
    observed = RealizedTranscriptProduct(
        view, sequence, protein, product.junction_signature, start_present,
        {**product.evidence, "genomic_positions": positions})
    return ORFProductAnnotation(reference=reference,
                status="conditional_translation", effect=effect,
                product=_product(view, observed), cdna_sequence=sequence,
                evidence=product.evidence)


def annotate_orf_product(reference, product):
    """Compatibility JSON export; prefer :func:`annotate_orf_rna_product`.

    This original entry point returns a dictionary because downstream RNA
    exporters embed its result directly in JSON. New Python callers should use
    ``annotate_orf_rna_product`` and access the returned object's attributes.
    """
    return annotate_orf_rna_product(reference, product).to_dict()


def annotate_orf_variants(variants, reference, **kwargs):
    """Annotate variants relative to a validated ``orfdata.ReferenceORF``.

    The reference's coding fields and ``reference.evidence`` describe the base
    ORF; catalogue ``record``/``model`` are optional. The returned annotation retains that reference and all alternative
    modified products. The parent RNA can be annotated protein-coding or
    noncoding. Native peptide detection never establishes mutant translation.
    """
    return _annotate_orf_variants(variants, reference.record, reference.transcript,
                                  genome_build=reference.context.genome_build,
                                  reference=reference, **kwargs)
