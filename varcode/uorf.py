"""Alternative upstream frames using Varcode's transcript-model engine.

The optional oncoref catalogue supplies reference observations, never patient
translation evidence. Canonical ``effects()`` results are unchanged.
"""

from dataclasses import asdict, dataclass, fields, is_dataclass
from typing import Mapping

from .effect_hypotheses import RealizedTranscriptProduct
from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import codon_table_for_transcript, translate_sequence
from .effects.effect_classes import StartLoss, Unresolved
from .genomic_layout import SequenceUnavailable, UnsupportedLayoutEdit
from .mutant_transcript import MutantTranscript
from .transcript_model import predict_transcript_model_effect


def load_uorf_catalog(*args, **kwargs):
    """Load the shared oncoref reference; install ``varcode[uorf]`` to use it."""
    try:
        from oncoref import load_uorf_catalog as load
    except ImportError as error:
        raise ImportError("Install varcode[uorf] for the shared reference catalogue") from error
    return load(*args, **kwargs)


def _contig(value):
    value = str(value).removeprefix("chr")
    return "MT" if value == "M" else value


def _positions(transcript):
    exons = sorted(transcript.exons, key=lambda e: e.start,
                   reverse=transcript.on_backward_strand)
    return tuple(p for e in exons for p in (
        range(e.end, e.start - 1, -1) if transcript.on_backward_strand
        else range(e.start, e.end + 1)))


def _model_positions(record):
    positions = tuple(p for a, b in record.coding_blocks for p in range(a, b + 1))
    return positions[::-1] if record.strand == "-" else positions


class ORFTranscript:
    """Validated alternative coding frame on one complete reference RNA path.

    Transcript ID is not sufficient: every spliced coding base, the stop and
    the protein must agree. A newer annotation may be used only after these
    checks; the resolved release/version is retained in the annotation export.
    """

    def __init__(self, record, transcript, *, genome_build):
        if genome_build != record.genome_build:
            raise ValueError("ORF assembly mismatch: %s != %s" % (genome_build, record.genome_build))
        reference_build = getattr(getattr(transcript, "genome", None), "reference_name", genome_build)
        if reference_build != genome_build:
            raise ValueError("Transcript assembly mismatch")
        if str(transcript.id).split(".")[0] != record.transcript_id:
            raise ValueError("ORF transcript identity mismatch")
        if (_contig(transcript.contig) != _contig(record.contig)
                or transcript.on_backward_strand != (record.strand == "-")):
            raise ValueError("ORF transcript locus mismatch")
        self.record = record
        self.reference_transcript = transcript
        self.positions = _positions(transcript)
        self.sequence = str(transcript.sequence).upper()
        if len(self.sequence) != len(self.positions):
            raise ValueError("Transcript sequence length mismatch")
        wanted = _model_positions(record)
        try:
            start = self.positions.index(wanted[0])
        except ValueError as error:
            raise ValueError("ORF splice path mismatch") from error
        if self.positions[start:start + len(wanted)] != wanted:
            raise ValueError("ORF splice path mismatch")
        coding = self.sequence[start:start + len(wanted)]
        table = codon_table_for_transcript(transcript)
        if any(b not in "ACGT" for b in coding):
            raise ValueError("ORF reference sequence is ambiguous")
        protein = translate_sequence(coding, codon_table=table, to_stop=True)
        # Preserve the source's initiation convention, including its accepted
        # near-cognate starts. No new initiation efficiency is inferred.
        if coding[:3] in table.start_codons and protein:
            protein = "M" + protein[1:]
        if (len(wanted) != 3 * (len(record.protein_sequence) + 1)
                or coding[-3:] not in table.stop_codons
                or protein != record.protein_sequence):
            raise ValueError("ORF reference protein/stop mismatch")
        self.protein_sequence = protein
        self.start_codon_spliced_offsets = tuple(range(start, start + 3))
        stop = start + len(wanted) - 3
        self.stop_codon_spliced_offsets = tuple(range(stop, stop + 3))
        self.start_codon_positions = wanted[:3]
        self.stop_codon_positions = wanted[-3:]
        self.coding_sequence = coding
        self.coding_sequence_position_ranges = record.coding_blocks
        self.five_prime_utr_sequence = self.sequence[:start]
        self.three_prime_utr_sequence = self.sequence[stop + 3:]
        self.is_protein_coding = True
        self.complete = True
        # Model identity must survive candidate deduplication.
        self.id = "%s:%s" % (transcript.id, record.orf_id)
        self.name = self.id

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "reference_transcript"), name)


def _reference_provider(view, fallback=None):
    offsets = {p: i for i, p in enumerate(view.positions)}

    def sequence(contig, start, end):
        if _contig(contig) == _contig(view.contig) and all(p in offsets for p in range(start, end + 1)):
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
    if start is None or sequence is None or start + 3 > len(sequence):
        return {"status": "unavailable", "reason": "reference_start_not_retained"}
    # Raw contexts remain usable with any supplied ORF record, without oncoref.
    return dict(upstream=sequence[max(0, start - 9):start],
                start_codon=sequence[start:start + 3], downstream=sequence[start + 3:start + 9],
                start_offset=start, five_prime_complete=None,
                interpretation="sequence_prior_only; initiation_not_measured")


def _peptide_overlay(record, protein, origins, start):
    """Track native residues by codon ancestry, not coincidental peptide matches."""
    native = _model_positions(record)
    retained = []
    if protein is not None and start is not None:
        by_codon = {tuple(native[i:i + 3]): i // 3 + 1
                    for i in range(0, len(native) - 3, 3)}
        for index, aa in enumerate(protein):
            codon = tuple(origins[start + index * 3:start + index * 3 + 3])
            ref = by_codon.get(codon)
            if ref is not None and record.protein_sequence[ref - 1] == aa:
                retained.append((ref, index + 1))
    mapping = dict(retained)
    observations = []
    for peptide in record.peptides:
        intervals = []
        for a, b in peptide.protein_intervals:
            mapped = [mapping.get(i) for i in range(a, b + 1)]
            if (mapped and all(i is not None for i in mapped)
                    and mapped == list(range(mapped[0], mapped[0] + len(mapped)))):
                intervals.append((mapped[0], mapped[-1]))
        observations.append(dict(reference_observation=asdict(peptide),
                                 retained_mutant_intervals=intervals,
                                 reference_peptide_retained=bool(intervals) if protein is not None else None))
    return dict(retained_native_residues=retained, peptide_observations=observations,
                reference_ms_support_on_retained_sequence=any(
                    o["reference_peptide_retained"] and p.is_supporting
                    for o, p in zip(observations, record.peptides)),
                mutant_translation_evidence="not_assessed",
                mutant_hla_presentation_evidence="not_assessed")


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
    return dict(cdna_sequence=sequence, protein_sequence=protein,
                start_codon_present=product.start_codon_present,
                termination="observed_stop" if stop is not None else "sequence_boundary" if sequence is not None else "unavailable",
                stop_cdna_offset=stop, initiation_context=context,
                junction_signature=_plain(product.junction_signature),
                evidence=_plain(product.evidence),
                reference_evidence=_peptide_overlay(view.record, protein, origins, start))


def _effect(effect):
    return {"type": type(effect).__name__, **{
        k: _plain(getattr(effect, k)) for k in (
            "aa_mutation_start_offset", "aa_mutation_end_offset", "aa_ref", "aa_alt",
            "reason", "mechanism", "max_hypotheses") if hasattr(effect, k)}}


@dataclass
class UORFAnnotation:
    """Per-model result, preserving all realized phase/splice alternatives."""

    record: object
    effect: object
    data: dict

    @property
    def orf_id(self):
        return self.record.orf_id

    def to_dict(self):
        return _plain(self.data)


def annotate_uorf_variants(variants, record, transcript, *, genome_build,
                           germline_variants=(), phase_resolver=None,
                           sequence_provider=None, max_hypotheses=64,
                           homozygous_germline=(), metadata=None):
    """Annotate one known-cis group in a validated upstream frame.

    Literal REF mismatches raise. Unavailable intronic sequence and unsupported
    events remain unresolved. A caller-supplied transcript must have the exact
    model's spliced path and protein; canonical labels are not replaced.
    """
    variants = tuple(variants)
    if not variants:
        raise ValueError("At least one variant is required")
    view = ORFTranscript(record, transcript, genome_build=genome_build)
    germline_variants = tuple(germline_variants)
    provider = _reference_provider(view, sequence_provider)
    for variant in (*variants, *germline_variants):
        if variant.reference_name != genome_build:
            raise ValueError("Variant assembly mismatch")
        if _contig(variant.contig) != _contig(record.contig):
            raise ValueError("Variant contig mismatch")
    try:
        # Validate original genomic REF, including VCF padding anchors. Layout
        # tier 0 deliberately permits unknown intronic REF; this API does not.
        from .errors import ReferenceMismatchError
        for variant in (*variants, *germline_variants):
            if not getattr(variant, "is_structural", False) and variant.original_ref:
                expected = provider(variant.contig, variant.original_start, variant.original_start + len(variant.original_ref) - 1)
                if expected != variant.original_ref:
                    raise ReferenceMismatchError(variant, transcript, variant.original_ref, expected)
        effect = predict_transcript_model_effect(
            variants, view, germline_variants=germline_variants,
            phase_resolver=phase_resolver,
            sequence_provider=provider,
            max_hypotheses=max_hypotheses, homozygous_germline=homozygous_germline)
    except (SequenceUnavailable, UnsupportedLayoutEdit) as error:
        effect = Unresolved(variants[0], view, mechanism="uorf_sequence_unavailable", reason=str(error))
    data = dict(schema="varcode.uorf_annotation.v1", model=record.to_dict(),
                catalogue_metadata=metadata, variants=_plain(variants),
                resolved_transcript=dict(id=transcript.id,
                    version=getattr(transcript, "version", None),
                    release=getattr(getattr(transcript, "genome", None), "release", None),
                    genome_build=genome_build), effect=_effect(effect), outcomes=[],
                discovery=dict(status="reference_models_only", new_starts="not_assessed"),
                downstream_cds_regulation="not_assessed")
    data["variant_locations"] = []
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
        data["variant_locations"].append(dict(variant=_plain(variant), locations=locations,
            initiation_snvs=context_snvs,
            insertion_after_position=variant.start if getattr(variant, "is_insertion", False) else None))
    for candidate in getattr(effect, "realized_candidates", ()):
        for outcome in candidate.outcomes:
            data["outcomes"].append(dict(
                effect=_effect(outcome.effect), hypothesis=_plain(outcome.hypothesis),
                baseline=_product(view, outcome.baseline), mutant=_product(view, outcome.mutant),
                probability=candidate.probability, phase_mass=candidate.phase_mass))
    if not data["outcomes"] and getattr(variants[0], "is_structural", False):
        baseline = RealizedTranscriptProduct(
            view, view.sequence, view.protein_sequence, start_codon_present=True,
            evidence={"genomic_positions": view.positions})
        for candidate in getattr(effect, "candidates", ()):
            alternative = candidate.effect
            mt = getattr(alternative, "mutant_transcript", None)
            if mt is None or mt.cdna_sequence is None:
                data["outcomes"].append(dict(effect=_effect(alternative),
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
                data["outcomes"].append(dict(effect=_effect(alternative),
                    baseline=_product(view, baseline), mutant=None,
                    unresolved_reason="structural_base_origins_unavailable"))
                continue
            product = RealizedTranscriptProduct(
                view, mt.cdna_sequence, mt.mutant_protein_sequence,
                evidence={**(mt.evidence or {}), "base_origins": origins})
            observed = annotate_uorf_product(record, transcript, product, genome_build=genome_build)
            data["outcomes"].append(dict(effect=_effect(alternative),
                baseline=_product(view, baseline), mutant=observed.get("product"),
                source=getattr(candidate, "source", None),
                unresolved_reason=observed.get("reason")))
    return UORFAnnotation(record, effect, data)


def predict_uorf_effects(variants, *, catalog=None, transcripts=None, **kwargs):
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
    catalog = load_uorf_catalog() if catalog is None else catalog
    genome = variants[0].genome
    build = variants[0].reference_name
    if build != catalog.genome_build:
        raise ValueError("Catalogue assembly mismatch")
    by_id = None if transcripts is None else {str(t.id).split(".")[0]: t for t in transcripts}
    results = []
    for record in catalog.records:
        # Include splice sites inside the ORF envelope and initiation context.
        positions = (record.coding_blocks[0][0], record.coding_blocks[-1][1])
        context = record.initiation_context
        if context is not None:
            positions += tuple(context.genomic_positions)
        if not any(_contig(v.contig) == _contig(record.contig)
                   and v.start <= max(positions) + 6 and v.end >= min(positions) - 6
                   for v in variants):
            continue
        try:
            transcript = genome.transcript_by_id(record.transcript_id) if by_id is None else by_id[record.transcript_id]
            annotation = annotate_uorf_variants(
                variants, record, transcript, genome_build=build,
                metadata=catalog.metadata, **kwargs)
        except (KeyError, ValueError) as error:
            from .errors import ReferenceMismatchError
            if isinstance(error, ReferenceMismatchError):
                raise
            effect = Unresolved(variants[0], None, mechanism="uorf_model_unavailable", reason=str(error))
            annotation = UORFAnnotation(record, effect, dict(
                schema="varcode.uorf_annotation.v1", model=record.to_dict(),
                catalogue_metadata=catalog.metadata, variants=_plain(variants),
                effect=_effect(effect), outcomes=[],
                discovery=dict(status="reference_models_only", new_starts="not_assessed")))
        results.append(annotation)
    return tuple(results)


def annotate_uorf_product(record, transcript, product, *, genome_build):
    """Annotate an observed/reconstructed RNA product with mapped base ancestry.

    ``product`` is a :class:`RealizedTranscriptProduct`. Its evidence must carry
    ``base_origins``: one ``(contig, position, strand)`` or None per cDNA base.
    These placements must come from the reconstruction, not peptide similarity.
    RNA support in ``product.evidence`` remains separate from reference MS.
    Only a uniquely retained complete reference start anchors translation.
    """
    view = ORFTranscript(record, transcript, genome_build=genome_build)
    sequence = product.cdna_sequence
    origins = product.evidence.get("base_origins", ())
    if sequence is None or len(sequence) != len(origins):
        raise ValueError("A sequence and one base origin per nucleotide are required")
    positions = tuple(o[1] if o is not None and _contig(o[0]) == _contig(record.contig)
                      and o[2] == record.strand else None for o in origins)
    starts = [i for i in range(len(positions) - 2)
              if positions[i:i + 3] == tuple(view.start_codon_positions)]
    if len(starts) != 1:
        return dict(schema="varcode.uorf_product.v1", model=record.to_dict(),
                    status="unresolved", reason="reference_start_missing_or_duplicated",
                    cdna_sequence=sequence, evidence=_plain(product.evidence))
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
                          mutant_protein_sequence=protein, evidence=product.evidence)
    # The shared classifier compares protein outcomes; DNA linkage and source
    # placements remain in the RNA evidence supplied by the caller.
    effect = (classify_from_protein_diff(None, product_view, record.protein_sequence, protein,
                                       mutant_transcript=mt) if start_present else StartLoss(None, view))
    observed = RealizedTranscriptProduct(
        view, sequence, protein, product.junction_signature, start_present,
        {**product.evidence, "genomic_positions": positions})
    return dict(schema="varcode.uorf_product.v1", model=record.to_dict(),
                status="conditional_translation", effect=_effect(effect),
                product=_product(view, observed))
