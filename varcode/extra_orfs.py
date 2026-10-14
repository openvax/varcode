"""Supply additional reference coding frames to the ordinary effects API."""

from .effects.effect_classes import Failure
from .orf import (
    ORFAnnotation, _annotate_orf_variants, _bind, _contig, _orf_key,
)
from .version import __version__


def _overlaps(variant, record):
    # Include splice sites between coding blocks and the source's start context.
    positions = [p for block in record.coding_blocks for p in block]
    if record.initiation_context is not None:
        positions.extend(record.initiation_context.genomic_positions)
    return (_contig(variant.contig) == _contig(record.contig)
            and variant.start <= max(positions) + 6
            and variant.end >= min(positions) - 6)


class ExtraORFInput:
    """Snapshot caller iterables once; preserve collection scope and exclusions."""

    def __init__(self, value):
        from orfdata import ORFCatalog, ReferenceORF

        self.catalog = value if isinstance(value, ORFCatalog) else None
        self.context = getattr(value, "reference", None) if self.catalog is None else None
        self.metadata = getattr(value, "metadata", None)
        self.excluded = ()
        self.references = ()
        if self.catalog is not None:
            return
        self.references = (value,) if isinstance(value, ReferenceORF) else tuple(value)
        if any(not isinstance(ref, ReferenceORF) for ref in self.references):
            raise TypeError("extra_orfs requires an ORFdata catalogue or bound reference ORFs")
        if any(ref.record is None for ref in self.references):
            raise NotImplementedError("Annotation-derived ORFs require source coding models")
        catalog = getattr(value, "catalog", None)
        if catalog is not None:
            self.metadata = catalog.metadata
        # ORFdata collections retain records that could not be placed. Do not
        # erase those failures by coercing the collection to just its iterable.
        self.excluded = tuple((catalog.get(item.orf_id), item.reason)
                              for item in getattr(value, "excluded", ()))

    def validate(self, genome):
        if self.catalog is not None:
            from orfdata.reference import canonical_species, species_taxon

            if self.catalog.genome_build != genome.reference_name:
                raise ValueError("Extra ORF catalogue assembly mismatch")
            species = getattr(genome, "species", None)
            taxon = getattr(genome, "taxon_id", None) or species_taxon(species)
            if taxon is not None and taxon != self.catalog.taxon_id:
                raise ValueError("Extra ORF catalogue taxon mismatch")
            if species is not None and canonical_species(species) != canonical_species(self.catalog.species):
                raise ValueError("Extra ORF catalogue species mismatch")
        if self.context is not None:
            self.context.validate_genome(genome)
        for reference in self.references:
            reference.context.validate_genome(genome)
            reference.context.validate_genome(reference.transcript.genome)

    def effects(self, variant, *, raise_on_error=True, germline=None,
                phase_resolver=None, rna_resolver=None):
        self.validate(variant.genome)
        if self.catalog is not None:
            entries = ((record, None, None) for record in self.catalog.records)
        else:
            entries = [(reference.record, reference, None) for reference in self.references]
            entries.extend((record, None, reason) for record, reason in self.excluded)
        result, seen = [], set()
        for record, reference, reason in entries:
            if not _overlaps(variant, record):
                continue
            key = (_orf_key(reference) if reference is not None else
                   (record.source, record.orf_id, record.transcript_id))
            if key in seen:
                continue
            seen.add(key)
            try:
                if reason is not None:
                    raise ValueError(reason)
                if reference is None:
                    transcript = variant.genome.transcript_by_id(record.transcript_id.split(".")[0])
                    reference = _bind(record, transcript, variant.reference_name)
                annotation = _annotate(
                    (variant,), reference, metadata=self.metadata,
                    germline=germline, phase_resolver=phase_resolver,
                    rna_resolver=rna_resolver)
                effect = annotation.effect
            except Exception as error:
                if raise_on_error:
                    raise
                effect = Failure(variant, error=str(error))
                annotation = ORFAnnotation(
                    record=record, effect=None, reference=reference,
                    variants=(variant,), catalogue_metadata=self.metadata,
                    unresolved_reason=str(error))
            effect.orf_annotation = annotation
            effect.annotator = "transcript_model"
            effect.annotator_version = __version__
            result.append(effect)
        return result


def prepare_extra_orfs(value):
    return value if value is None or isinstance(value, ExtraORFInput) else ExtraORFInput(value)


class ORFPhaseResolver:
    """Phase is about the original RNA, never its private coding-view ID."""

    def __init__(self, resolver):
        self.resolver = resolver
        self.phase_source = getattr(resolver, "phase_source", None)

    def in_cis(self, first, second, transcript=None):
        from .phasing import query_in_cis

        transcript = getattr(transcript, "reference_transcript", transcript)
        return query_in_cis(self.resolver, first, second, transcript)


def _annotate(variants, reference, *, metadata=None, germline=None,
              phase_resolver=None, rna_resolver=None):
    transcript = reference.transcript
    # The ordinary RNA resolver reports products for the annotated CDS. Those
    # cannot be attached to a different ORF without remapping base ancestry.
    for variant in variants:
        observed = rna_resolver.observed_outcomes(variant, transcript) if hasattr(
            rna_resolver, "observed_outcomes") else None
        if observed is not None and tuple(observed):
            raise NotImplementedError(
                "RNA outcomes for extra ORFs require annotate_orf_rna_product() with base origins")
        if phase_resolver is not None and hasattr(phase_resolver, "mutant_transcript"):
            if phase_resolver.mutant_transcript(variant, transcript) is not None:
                raise NotImplementedError(
                    "Observed RNA for extra ORFs requires annotate_orf_rna_product() with base origins")
    germline_variants = () if not germline else germline.variants_in_window(
        transcript.contig, transcript.start, transcript.end)
    annotation = _annotate_orf_variants(
        variants, reference.record, transcript, genome_build=reference.model.genome_build,
        reference=reference, metadata=metadata, germline_variants=germline_variants,
        homozygous_germline=() if not germline else germline.homozygous(germline_variants),
        max_hypotheses=64 if not germline else germline.max_phase_hypotheses,
        phase_resolver=None if phase_resolver is None else ORFPhaseResolver(phase_resolver))
    annotation.effect.orf_annotation = annotation
    return annotation


class ORFHaplotypeAnnotator:
    """Reuse collection phasing with the same ORF annotation engine."""

    name = "transcript_model"
    version = __version__

    def annotate_haplotype(self, variants, transcript, germline_ctx=None, phase_resolver=None):
        annotation = _annotate(
            variants, transcript.reference_orf, germline=germline_ctx,
            phase_resolver=getattr(phase_resolver, "resolver", phase_resolver))
        return annotation.effect
