"""Supply additional reference coding frames to the ordinary effects API."""

from bisect import bisect_left, bisect_right
from collections import OrderedDict

from .effects.effect_classes import Failure
from .orf import (
    ORFAnnotation, _annotate_orf_variants, _bind, _contig, _orf_key,
)
from .version import __version__


def _span(record):
    # Include splice sites between coding blocks and the source's start context.
    positions = [p for block in record.coding_blocks for p in block]
    if record.initiation_context is not None:
        positions.extend(record.initiation_context.genomic_positions)
    return _contig(record.contig), min(positions) - 6, max(positions) + 6


def _overlaps(variant, record):
    contig, start, end = _span(record)
    return _contig(variant.contig) == contig and variant.start <= end and variant.end >= start


class _SpanIndex:
    """Entries by contig, sorted by span start, so a variant never scans them all."""

    def __init__(self, entries):
        by_contig = {}
        for order, entry in enumerate(entries):
            contig, start, end = _span(entry[0])
            by_contig.setdefault(contig, []).append((start, end, order, entry))
        self._contigs = {}
        for contig, items in by_contig.items():
            items.sort(key=lambda item: (item[0], item[2]))
            self._contigs[contig] = ([item[0] for item in items], max(item[1] - item[0] for item in items), items)

    def overlapping(self, variant):
        found = self._contigs.get(_contig(variant.contig))
        if found is None:
            return []
        starts, widest, items = found
        window = items[bisect_left(starts, variant.start - widest):bisect_right(starts, variant.end)]
        return [item[3] for item in sorted(
            (item for item in window if item[1] >= variant.start), key=lambda item: item[2])]


class ExtraORFInput:
    """Snapshot caller iterables once; preserve collection scope and exclusions."""

    def __init__(self, value, *, metadata=None, report_excluded=True):
        from orfdata import ORFCatalog, ReferenceORF

        self.catalog = value if isinstance(value, ORFCatalog) else None
        self.context = getattr(value, "reference", None) if self.catalog is None else None
        self.metadata = getattr(value, "metadata", None)
        self.excluded = ()
        self.references = ()
        self._index = None
        self._validated = []
        if self.catalog is not None:
            return
        self.references = (value,) if isinstance(value, ReferenceORF) else tuple(value)
        if any(not isinstance(ref, ReferenceORF) for ref in self.references):
            raise TypeError("extra_orfs requires an ORFdata catalogue or bound reference ORFs")
        if any(ref.record is None for ref in self.references):
            raise NotImplementedError("Annotation-derived ORFs require source coding models")
        catalog = getattr(value, "catalog", None)
        if catalog is not None:
            self.metadata = catalog.metadata if metadata is None else metadata
        # ORFdata collections retain records that could not be placed. Do not
        # erase those failures by coercing the collection to just its iterable,
        # unless the caller only asked for ORFdata's ORFs on the genome.
        if report_excluded:
            self.excluded = tuple((catalog.get(item.orf_id), item.reason)
                                  for item in getattr(value, "excluded", ()))

    def validate(self, genome):
        if any(genome is seen for seen in self._validated):
            return
        self._validate(genome)
        self._validated = (self._validated + [genome])[-4:]

    def _validate(self, genome):
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
        if self._index is None:
            if self.catalog is not None:
                entries = [(record, None, None) for record in self.catalog.records]
            else:
                entries = [(reference.record, reference, None) for reference in self.references]
                entries.extend((record, None, reason) for record, reason in self.excluded)
            self._index = _SpanIndex(entries)
        result, seen = [], set()
        for record, reference, reason in self._index.overlapping(variant):
            key = (_orf_key(reference) if reference is not None else
                   (record.source, record.orf_id, record.transcript_id))
            if key in seen:
                continue
            seen.add(key)
            try:
                if reason is not None:
                    raise ValueError("Extra ORF %s: %s" % (record.orf_id, reason))
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
                unresolved = reason if reason is not None else str(error)
                effect = Failure(variant, error=unresolved)
                annotation = ORFAnnotation(
                    record=record, effect=None, reference=reference,
                    variants=(variant,), catalogue_metadata=self.metadata,
                    unresolved_reason=unresolved)
            effect.orf_annotation = annotation
            effect.annotator = "transcript_model"
            effect.annotator_version = __version__
            result.append(effect)
        return result


class GenomeExtraORFs:
    """ORFdata's extra ORFs for each variant's own genome, placed once per genome.

    ``evidence`` is what ``orfdata.extra_orfs`` takes: ``None`` for every
    catalogue ORF, an evidence name such as ``"ms"`` or
    ``"immunopeptidomics"``, or a rule such as ``RiboSeq() & MS()``. ORFdata picks the catalogue for the
    genome's species and assembly, lifting one from another assembly when
    needed, and places each ORF on every transcript carrying its exact
    coding path and protein.

    The ORFs ORFdata could not place on a genome (an older annotation may
    lack their transcript) are listed by :meth:`excluded`; a variant in one
    gets the ordinary effects, never an error. Placements for the last few
    genomes are kept; :meth:`clear` drops them.
    """

    max_genomes = 4

    def __init__(self, evidence=None):
        self.evidence = evidence
        self._inputs = OrderedDict()

    def _placed(self, genome):
        if genome in self._inputs:
            self._inputs.move_to_end(genome)
            return self._inputs[genome]
        try:
            from orfdata import extra_orfs
        except ImportError as error:
            raise ImportError("extra_orfs needs orfdata>=0.28.0: pip install 'varcode[orf]'") from error
        placed = extra_orfs(genome, self.evidence)
        found = (placed, ExtraORFInput(placed, metadata=_summary(placed.catalog), report_excluded=False))
        self._inputs[genome] = found
        while len(self._inputs) > self.max_genomes:
            self._inputs.popitem(last=False)
        return found

    def for_genome(self, genome):
        return self._placed(genome)[1]

    def excluded(self, genome):
        """ORFdata's ORFs that could not be placed on ``genome``, each with its reason."""
        placed = self._placed(genome)[0]
        return tuple(placed.excluded) + tuple(getattr(placed, "unlifted", ()))

    def clear(self):
        self._inputs.clear()

    def effects(self, variant, **kwargs):
        return self.for_genome(variant.genome).effects(variant, **kwargs)


def _summary(catalog):
    """The catalogue's identity and provenance, without per-record lists, for each annotation to carry."""
    metadata = catalog.metadata
    summary = {k: metadata[k] for k in ("genome_build", "taxon_id", "species", "coordinates", "dataset_version",
                                         "sources", "ribosome_study_pmids") if k in metadata}
    summary["source"] = catalog.source
    if metadata.get("merged_from"):
        summary["merged_from"] = [{k: m[k] for k in ("source", "records", "sha256") if k in m}
                                  for m in metadata["merged_from"]]
    if metadata.get("liftover"):
        summary["liftover"] = {k: v for k, v in metadata["liftover"].items() if k != "unmapped"}
    return summary


# One resolver per evidence selection, so repeated calls reuse placements.
_RESOLVERS = {}


def _resolver(evidence):
    import json
    try:
        from orfdata.evidence import evidence_rule
    except ImportError as error:
        raise ImportError("extra_orfs needs orfdata>=0.28.0: pip install 'varcode[orf]'") from error
    evidence = evidence_rule(evidence)  # "ms" and MS() share placements; an unknown name raises here.
    key = None if evidence is None else json.dumps(evidence.to_dict(), sort_keys=True)
    if key not in _RESOLVERS:
        _RESOLVERS[key] = GenomeExtraORFs(evidence)
    return _RESOLVERS[key]


def clear_extra_orfs_cache():
    """Drop the ORF placements kept for each genome when ``extra_orfs`` names ORFdata's ORFs."""
    for resolver in _RESOLVERS.values():
        resolver.clear()


def _is_rule(value):
    try:
        from orfdata import EvidenceRule
    except ImportError:
        return False
    return isinstance(value, EvidenceRule)


def _is_names(value):
    """A non-empty list, tuple or set of evidence names or rules: ORFs with any of them."""
    return (isinstance(value, (list, tuple, set, frozenset)) and len(value) > 0
            and all(isinstance(v, str) or _is_rule(v) for v in value))


def prepare_extra_orfs(value):
    """Normalize ``extra_orfs``: None or False for none; ``"all"``, an evidence
    name, a list of names (any of them) or an ORFdata rule for ORFdata's ORFs on
    each variant's genome; else a catalogue or bound ORFs."""
    if type(value).__module__ == "numpy" and type(value).__name__ in ("bool_", "bool"):
        value = bool(value)
    if value is None or value is False:
        return None
    if isinstance(value, (ExtraORFInput, GenomeExtraORFs)):
        return value
    if value is True:
        raise ValueError(
            "extra_orfs=True doesn't say which ORFs: pass 'all' for every ORF ORFdata has for the "
            "variant's genome, an evidence name such as 'ms', 'ribosome', 'proteomics' or "
            "'immunopeptidomics', a list of names for any of them, an ORFdata rule, or a catalogue")
    if value == "all":
        return _resolver(None)
    if isinstance(value, str) or _is_rule(value) or _is_names(value):
        return _resolver(value)
    return ExtraORFInput(value)


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
