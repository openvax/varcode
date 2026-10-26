"""Export native predictions to ORFdata without reannotating or translating.

The ``varcode.proteome.v1`` producer graph has one event per supplied native
result and one comparison per outcome. Baselines live in comparison records,
not in the search collection. The native graph is an inert JSON archive;
ORFdata's loader never reconstructs Python classes from it.
"""

import hashlib
import json
from dataclasses import replace
from typing import TYPE_CHECKING, Iterable

from serializable.helpers import to_json

from .effect_candidates import EffectCandidate
from .effect_hypotheses import ClassifiedOutcome, RealizedEffectCandidate
from .effects import EffectCollection, Failure, HypothesisLimit, MutationEffect, Unresolved
from .effects.structural_serialization import structural_effect_to_dict
from .orf import ORFAnnotation, ORFOutcome
from .orf_serialization import _ReferenceGenome

if TYPE_CHECKING:
    from orfdata import Proteome


_SCHEMA = "varcode.proteome.v1"
_EXPORTER = "varcode"
_COORDINATES = {
    "aa_intervals": "0-based half-open; empty mutant intervals retain deletion junctions",
    "cdna_offsets_and_edits": "0-based; edit intervals are half-open",
    "reference_segments": "0-based half-open in each segment source's coordinate system",
    "variants": "native Varcode coordinates; small variants are 1-based, inclusive",
}
_BUILTIN_ANNOTATORS = {"fast", "protein_diff", "realized", "transcript_model"}


def _json(value):
    """Retain native extension fields, rejecting values without a portable form."""
    def portable(item):
        if isinstance(item, BaseException):
            return dict(exception_type="%s.%s" % (type(item).__module__, type(item).__qualname__),
                        message=str(item), args=portable(item.args))
        if isinstance(item, dict):
            return {key: portable(child) for key, child in item.items()}
        if isinstance(item, (tuple, list)):
            return [portable(child) for child in item]
        return item
    return json.loads(to_json(portable(value)))


def _archive(value, *, comparison=False):
    schema = "orf_comparison_schema" if comparison else "orf_effect_schema"
    return _json(structural_effect_to_dict(value, schema=schema))


def _digest(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(data.encode("utf-8")).hexdigest()[:20]


def _identifier(label, identity):
    # The shared exporter escapes native IDs exactly once for each format.
    return "%s:%s" % (label[:160], _digest(identity))


def _ids(values):
    return tuple(sorted({str(value) for value in values if value is not None and str(value)}))


def _effect(value):
    return value if isinstance(value, MutationEffect) else getattr(value, "effect", None)


def _variants(effect, annotation=None):
    if annotation is not None:
        return annotation.variants
    return getattr(effect, "variants", ()) or (
        (effect.variant,) if effect is not None else ())


def _model(effect, annotation=None):
    bound = getattr(effect, "orf", None)
    return annotation.record if annotation is not None else None if bound is None else bound.record


def _genes(effect, transcripts=()):
    result: list = []
    for source in (effect, *transcripts):
        gene = getattr(source, "gene", None)
        if gene is not None and all(gene is not previous for previous in result):
            result.append(gene)
    return result


def _outcome_contexts(value, contexts):
    return (*contexts, dict(source=value.source)) if isinstance(value, ORFOutcome) and value.source is not None else contexts


def _limits(effect):
    return dict(max_hypotheses=effect.max_hypotheses) if isinstance(effect, HypothesisLimit) else {}


def _label(effect, annotation=None, kind=None):
    variants = _variants(effect, annotation)
    variant_labels = ["%s:%s:%s>%s" % (v.contig, v.start, v.ref, v.alt) for v in variants]
    record = _model(effect, annotation)
    models = ([record.orf_id] if record is not None else
              [t.id for t in _transcripts(effect)] or list(_ids(g.id for g in _genes(effect))))
    return ":".join([kind or type(effect).__name__, *models, *variant_labels])


def _variant_identity(variant):
    # to_dict includes original alleles and structural breakend orientation.
    state = dict(variant.to_dict())
    state["genome"] = _ReferenceGenome.capture(variant.genome).metadata
    return _json(state)


def _transcript_identity(transcript):
    return dict(id=transcript.id, genome=_ReferenceGenome.capture(transcript.genome).metadata)


def _transcripts(effect, annotation=None):
    result = []
    if annotation is not None and annotation.reference is not None:
        result.append(annotation.reference.transcript)
    for name in ("reference_transcript", "partner_transcript",
                 "five_prime_transcript", "three_prime_transcript"):
        transcript = getattr(effect, name, None)
        if transcript is not None and all(_transcript_identity(transcript) != _transcript_identity(t) for t in result):
            result.append(transcript)
    return result


def _mechanism_context(effect):
    """Source coordinates and model evidence, excluding predicted products."""
    context = {name: getattr(effect, name) for name in (
        "cryptic_genomic_position", "exon_length_delta", "in_frame",
        "retained_intron_start", "retained_intron_end", "side",
        "contig", "interval_start", "interval_end", "distance_to_exon")
        if getattr(effect, name, None) is not None}

    def feature_identity(feature):
        return {name: getattr(feature, name, None) for name in ("id", "contig", "start", "end", "strand")}

    for name in ("affected_exon", "nearest_exon", "exon"):
        feature = getattr(effect, name, None)
        if feature is not None:
            context[name] = feature_identity(feature)
    for name in ("affected_exons", "exons"):
        features = getattr(effect, name, None)
        if features is not None:
            context[name] = [feature_identity(feature) for feature in features]
    for name in ("five_prime_transcript", "three_prime_transcript"):
        transcript = getattr(effect, name, None)
        if transcript is not None:
            context[name] = _transcript_identity(transcript)
    model = getattr(effect, "mutant_transcript", None)
    if model is not None:
        context["transcript_model"] = dict(
            annotator=model.annotator_name, evidence=_json(model.evidence),
            reference=None if model.reference_transcript is None else _transcript_identity(model.reference_transcript),
            edits=[dict(start=edit.cdna_start, end=edit.cdna_end,
                        variant=None if edit.source_variant is None else _variant_identity(edit.source_variant))
                   for edit in model.edits],
            segments=None if model.reference_segments is None else [dict(
                start=segment.start, end=segment.end, strand=segment.strand, label=segment.label,
                source=dict(feature_identity(segment.source), genome=_ReferenceGenome.capture(
                    getattr(segment.source, "genome", segment.source)).metadata)) for segment in model.reference_segments])
    return context


def _germline_context(effect):
    context = {name: getattr(effect, name) for name in ("germline_phase_state", "germline_unknown")
               if hasattr(effect, name)}
    variants = getattr(effect, "germline_variants_in_window", None)
    if variants is not None:
        context["germline_variants_in_window"] = sorted((_variant_identity(v) for v in variants), key=_digest)
    return context


def _context(effect, annotation=None):
    """Validate available source context without filling missing reference fields."""
    from orfdata import ReferenceContext
    from orfdata.reference import species_taxon

    transcripts = _transcripts(effect, annotation)
    genomes = [source.genome for source in [*transcripts, *_genes(effect)]
               if getattr(source, "genome", None) is not None]
    variants = _variants(effect, annotation)
    variant_genomes = [v.genome for v in variants]
    bound = annotation.reference if annotation is not None else getattr(effect, "orf", None)
    reference = (getattr(annotation.record, "reference_context", None) if annotation is not None and bound is None
                 else None if bound is None else bound.context)
    all_genomes = genomes + variant_genomes
    taxa = set()
    for genome in all_genomes:
        taxon = (getattr(genome, "taxon_id", None) or species_taxon(getattr(genome, "species", None))
                 or getattr(genome, "species_taxon_id", None))
        if taxon is not None:
            taxa.add(taxon)
        if reference is None and taxon is not None and all(getattr(genome, key, None) is not None for key in
                                                         ("species", "reference_name", "annotation_name")):
            reference = ReferenceContext.from_genome(genome)
    if reference is not None:
        # ORFdata resolves reference vs. species taxa (dog vs. wolf, or a yeast
        # strain vs. its species). Use that resolution for both records and
        # entries, and validate even sources encountered before this context.
        for genome in all_genomes:
            reference.validate_genome(genome, include_annotation=any(genome is g for g in genomes))
        return reference, reference.taxon_id, transcripts
    # Without a qualified reference, retain explicit source taxa; do not replace
    # a reference organism by its species ancestor or invent reference metadata.
    if len(taxa) > 1:
        raise ValueError("A predicted product cannot silently combine different source taxa")
    return reference, next(iter(taxa), None), transcripts


def _source(effect, defaults, contexts, annotation=None):
    annotator = getattr(effect, "annotator", None)
    version = getattr(effect, "annotator_version", None)
    if annotator in (None, defaults.get("annotator")):
        annotator = annotator or defaults.get("annotator")
        version = version or defaults.get("annotator_version")
    sources = [item["source"] for item in contexts if item.get("source")]
    source = next((source for source in reversed(sources) if source != "varcode"),
                  sources[-1] if sources else None)
    # A source string such as "rna" describes an import path, not its generating tool.
    native_source = source in (None, "varcode", "varcode_splice", "varcode_motif", "varcode_germline")
    tool = "varcode" if native_source and (annotator in _BUILTIN_ANNOTATORS or source is not None) else "unknown"
    model = _model(effect, annotation)
    provider = source if not native_source else getattr(model, "source", None) or None
    return dict(method=(annotator if source in (None, "varcode") else source) or "unknown",
                tool=tool, tool_version=version if native_source else None, provider=provider)


def _product(product=None, *, effect=None, baseline=False):
    if product is not None:
        termination = getattr(product, "termination", None)
        stop = None if termination is None else {"observed_stop": True, "sequence_boundary": False}.get(termination)
        return dict(
            protein_sequence=product.protein_sequence,
            cdna_sequence=product.cdna_sequence,
            start_complete=product.start_codon_present, stop_complete=stop,
            native=_archive(product))
    sequence = (getattr(effect, "original_protein_sequence", None) if baseline else
                getattr(effect, "mutant_protein_sequence", None))
    transcript = getattr(effect, "reference_transcript" if baseline else "mutant_transcript", None)
    # Only explicit baseline features describe reference boundaries. They do
    # not prove the extent of a predicted mutant product.
    return dict(
        protein_sequence=sequence,
        cdna_sequence=getattr(transcript, "sequence" if baseline else "cdna_sequence", None),
        start_complete=getattr(transcript, "contains_start_codon", None) if baseline else None,
        stop_complete=getattr(transcript, "contains_stop_codon", None) if baseline else None)


def _status(effect, baseline, mutant, reason=None):
    if isinstance(effect, HypothesisLimit):
        return "hypothesis_limit", effect.reason
    if isinstance(effect, Failure):
        return "annotation_failure", str(effect.error) if effect.error is not None else "annotation failed"
    if isinstance(effect, Unresolved):
        return "unresolved", effect.reason or reason or effect.mechanism
    sequence = mutant["protein_sequence"]
    if sequence is None:
        return "sequence_unavailable", reason or "native result supplies no mutant protein sequence"
    if sequence == "":
        return "no_product", "native result supplies an empty protein sequence"
    if sequence == baseline["protein_sequence"]:
        return "unchanged", None
    if mutant["start_complete"] is False or mutant["stop_complete"] is False:
        return "partial_product", None
    return "product_available", None


def _intervals(effect):
    start = getattr(effect, "aa_mutation_start_offset", None)
    end = getattr(effect, "aa_mutation_end_offset", None)
    aa_ref = getattr(effect, "aa_ref", None)
    return dict(mutant_aa=None if start is None or end is None else [start, end],
                baseline_aa=None if start is None or aa_ref is None else [start, start + len(aa_ref)])


def _identity(value, defaults, contexts=()):
    """Source identity excludes the predicted sequence, permitting conflict checks."""
    annotation = value if isinstance(value, ORFAnnotation) else None
    effect = _effect(value)
    transcripts = _transcripts(effect, annotation)
    model = _model(effect, annotation)
    contexts = _outcome_contexts(value, contexts)
    identity = dict(
        type=type(value).__name__,
        variants=sorted((_variant_identity(v) for v in _variants(effect, annotation)), key=_digest),
        transcripts=[_transcript_identity(t) for t in transcripts],
        genes=[dict(id=gene.id, contig=getattr(gene, "contig", None),
                    start=getattr(gene, "start", None), end=getattr(gene, "end", None),
                    strand=getattr(gene, "strand", None),
                    genome=_ReferenceGenome.capture(getattr(gene, "genome", None)).metadata)
               for gene in _genes(effect) if not transcripts],
        model=None if model is None else dict(
            orf_id=model.orf_id, source=model.source, transcript_id=model.transcript_id,
            reference=None if model.reference_context is None else model.reference_context.key),
        source=_source(effect, defaults, contexts, annotation), contexts=contexts,
        mechanism=getattr(effect, "mechanism", None))
    for name, context in (("mechanism_context", _mechanism_context(effect)),
                          ("germline_context", _germline_context(effect))):
        if context:
            identity[name] = context
    if isinstance(value, (ClassifiedOutcome, ORFOutcome)):
        identity.update(hypothesis=None if value.hypothesis is None else _archive(value.hypothesis),
                        baseline=_product(value.baseline))
    if isinstance(value, RealizedEffectCandidate):
        identity["outcomes"] = sorted((_identity(o, defaults, contexts) for o in value.outcomes), key=_digest)
    return _identifier(_label(effect, annotation, type(value).__name__), identity)


class _Export:
    def __init__(self, sample_id, run_id, defaults):
        self.sample_id = sample_id
        self.run_id = run_id
        self.defaults = defaults
        self.proteins = []
        self.records = []
        self.links = []
        self.event_ids = {}
        self.references = {}
        self.unknown_reference = False

    def record(self, kind, native_id, **kwargs):
        from orfdata import ProducerRecord
        record = ProducerRecord(
            namespace="varcode", schema_version=_SCHEMA, kind=kind,
            native_id=native_id, exporters=(_EXPORTER,), sample_id=self.sample_id,
            run_id=self.run_id, **kwargs)
        self.records.append(record)
        return record

    def scope(self, effect, annotation=None):
        reference, taxon, transcripts = _context(effect, annotation)
        if reference is None:
            self.unknown_reference = True
        else:
            self.references[reference.key] = reference
        return reference, taxon, transcripts

    def link(self, source, target, relation, target_kind="record"):
        from orfdata import ProducerLink
        self.links.append(ProducerLink(source, target, target_kind, relation))

    def comparison(self, effect, event_id, *, baseline=None, mutant=None,
                   outcome=None, annotation=None, contexts=()):
        from orfdata import ProducerInfo, ProteinSequence

        supplied_products = outcome is not None
        baseline_data = _product(baseline, effect=effect, baseline=True)
        mutant_data = _product(mutant, effect=None if supplied_products else effect)
        reason = getattr(outcome, "unresolved_reason", None)
        status, omitted = _status(effect, baseline_data, mutant_data, reason)
        reference, taxon, transcripts = self.scope(effect, annotation)
        hypothesis = getattr(outcome, "hypothesis", None)
        record_model = _model(effect, annotation)
        genes = _genes(effect, transcripts)
        gene_ids = _ids([gene.id for gene in genes] + [getattr(record_model, "gene_id", None)])
        gene_names = _ids([gene.name for gene in genes] + [getattr(record_model, "gene_name", None)])
        payload = dict(
            effect_type=type(effect).__name__, effect_description=effect.short_description,
            native=_archive(effect, comparison=True),
            gene_ids=gene_ids, gene_names=gene_names,
            baseline=baseline_data, mutant=mutant_data,
            intervals=_intervals(effect), coordinate_conventions=_COORDINATES,
            hypothesis=None if hypothesis is None else _archive(hypothesis),
            contexts=contexts, omission_reason=omitted,
            focal_variants=[_variant_identity(v) for v in _variants(effect, annotation)],
            probability=getattr(outcome, "probability", None),
            phase_mass=getattr(outcome, "phase_mass", None),
            change_aa_offset=getattr(outcome, "change_aa_offset", None),
            change_cdna_offset=getattr(outcome, "change_cdna_offset", None),
            source_reference=None if reference is None else reference.to_dict(),
            source_taxon_id=taxon)
        # Full products, baselines and assumptions distinguish comparisons;
        # the parent event independently detects conflicting source results.
        native_id = _identifier(_label(effect, annotation) + ":" + effect.short_description,
                                dict(event=self.event_ids[event_id], comparison=payload))
        limits = _limits(effect)
        record = self.record("comparison", native_id, status=status, payload=payload,
                             complete=False if limits else None, limits=limits)
        self.link(event_id, record.key, "has_outcome")
        sequence = mutant_data["protein_sequence"]
        if not sequence:
            return status
        origin = ("nc_orf" if record_model is not None else "reference_cds" if any(
            getattr(t, "is_protein_coding", False) for t in transcripts) else "unknown")
        producer = ProducerInfo(
            origin=origin, native_id=native_id, **_source(effect, self.defaults, contexts, annotation),
            exporters=(_EXPORTER,), sample_id=self.sample_id, run_id=self.run_id,
            gene_ids=gene_ids, gene_names=gene_names,
            transcript_ids=_ids([t.id for t in transcripts] + [getattr(record_model, "transcript_id", None)]),
            orf_ids=() if record_model is None else (record_model.orf_id,),
            orf_kind=getattr(record_model, "biotype", None), reference=reference,
            start_complete=mutant_data["start_complete"], stop_complete=mutant_data["stop_complete"])
        protein = ProteinSequence.from_producer(sequence, taxon_id=taxon, producer=producer)
        self.proteins.append(protein)
        self.link(record.key, protein.accession, "produces", "protein")
        return status

    def walk(self, value, event_id, *, contexts=(), annotation=None, ancestors=()):
        if isinstance(value, EffectCandidate):
            context = dict(source=value.source, evidence=_json(dict(value.evidence)))
            return self.walk(value.effect, event_id, contexts=(*contexts, context),
                             annotation=annotation, ancestors=ancestors)
        if isinstance(value, RealizedEffectCandidate):
            context = dict(ordinal_key=value.ordinal_key, probability=value.probability,
                           phase_mass=value.phase_mass, first_enumeration_index=value.first_enumeration_index)
            return [status for outcome in value.outcomes for status in self.walk(
                outcome, event_id, contexts=(*contexts, context), annotation=annotation, ancestors=ancestors)]
        if isinstance(value, (ClassifiedOutcome, ORFOutcome)):
            return [self.comparison(value.effect, event_id, baseline=value.baseline,
                                    mutant=value.mutant, outcome=value, annotation=annotation,
                                    contexts=_outcome_contexts(value, contexts))]
        if isinstance(value, ORFAnnotation):
            if value.outcomes:
                return [status for outcome in value.outcomes for status in self.walk(
                    outcome, event_id, annotation=value, contexts=contexts, ancestors=ancestors)]
            if value.reference is None:
                status, reason = "failed_binding", value.unresolved_reason or "reference ORF binding unavailable"
            elif value.effect is None:
                status, reason = "unresolved", value.unresolved_reason or "native annotation supplies no effect or outcomes"
            else:
                status, reason = _status(value.effect, {"protein_sequence": None},
                                         {"protein_sequence": None}, value.unresolved_reason)
            limits = _limits(value.effect)
            record = self.record("comparison", self.event_ids[event_id] + ":" + status, status=status,
                                 limits=limits, complete=False if limits else None,
                                 payload=dict(orf_id=value.orf_id, omission_reason=reason,
                                              contexts=contexts))
            self.link(event_id, record.key, "has_outcome")
            return [status]
        if not isinstance(value, MutationEffect):
            raise TypeError("to_proteome expects native effects, candidates, or ORFAnnotations")
        if value.orf_annotation is not None and annotation is None:
            return self.walk(value.orf_annotation, event_id, contexts=contexts, ancestors=ancestors)
        if id(value) in ancestors:
            # Self-candidates expose the primary outcome of structural effects.
            return [self.comparison(value, event_id, contexts=contexts, annotation=annotation)]
        ancestors = (*ancestors, id(value))
        candidates = getattr(value, "realized_candidates", None)
        if candidates is not None:
            candidates = tuple(candidates) + tuple(getattr(value, "_extra_candidates", ()))
        else:
            candidates = getattr(value, "candidates", None)
        if candidates:
            return [status for candidate in candidates for status in self.walk(
                candidate, event_id, contexts=contexts, annotation=annotation, ancestors=ancestors)]
        return [self.comparison(value, event_id, contexts=contexts, annotation=annotation)]

    def add(self, value):
        contexts: tuple = ()
        if isinstance(value, EffectCandidate):
            contexts = (dict(source=value.source, evidence=_json(dict(value.evidence))),)
            value = value.effect
        if isinstance(value, MutationEffect) and value.orf_annotation is not None:
            value = value.orf_annotation
        if not isinstance(value, (MutationEffect, ORFAnnotation, RealizedEffectCandidate,
                                  ClassifiedOutcome, ORFOutcome)):
            raise TypeError("to_proteome expects native effects, candidates, or ORFAnnotations")
        self.scope(_effect(value), value if isinstance(value, ORFAnnotation) else None)
        native_id = _identity(value, self.defaults, contexts)
        # The temporary record provides the shared source-qualified graph key.
        event_index = len(self.records)
        event = self.record("event", native_id, status="supplied_result")
        self.event_ids[event.key] = native_id
        statuses = self.walk(value, event.key, contexts=contexts)
        outcome_contexts = _outcome_contexts(value, contexts)
        payload = dict(native=_archive(value), contexts=outcome_contexts,
                       source=_source(_effect(value), self.defaults, outcome_contexts,
                                      value if isinstance(value, ORFAnnotation) else None),
                       coordinate_conventions=_COORDINATES,
                       sequence_convention="supplied amino acids, including any supplied terminal *; no trimming",
                       coverage="supplied outcomes only; hypothesis enumeration completeness unknown")
        # Replace the provisional event without disturbing links to its key.
        self.records[event_index] = replace(
            event, status=("no_supplied_outcomes" if not statuses else
                           statuses[0] if len(set(statuses)) == 1 else "alternative_outcomes"),
            complete=False if "hypothesis_limit" in statuses else None,
            limits=_limits(_effect(value)), payload=payload)
        return event.key


def to_proteome(effects: Iterable | MutationEffect | EffectCandidate | RealizedEffectCandidate |
                ClassifiedOutcome | ORFAnnotation | ORFOutcome = (), *, sample_id: str | None = None,
                run_id: str | None = None, orf_annotations: Iterable[ORFAnnotation] = ()) -> "Proteome":
    """Expose supplied native outcomes through ORFdata's shared collection API.

    Parameters
    ----------
    effects : native result, EffectCollection or iterable
        Native effects, EffectCandidates, RealizedEffectCandidates,
        ClassifiedOutcomes, ORFAnnotations or ORFOutcomes. A single result is
        accepted directly. Every realized outcome is retained, including
        failed ORF bindings with no effect.
    sample_id, run_id : str, optional
        Explicit source scope. Predictions without a sample are supported.
    orf_annotations : iterable of ORFAnnotation
        Standalone ORF annotations, including failed bindings with no effect.

    Returns
    -------
    orfdata.Proteome
        Available mutant sequences and a versioned graph of native results,
        baseline comparisons, assumptions and omitted outcomes. Baselines need
        not be entries. Equal sequences retain distinct source associations.

    Notes
    -----
    Requires the optional ``orf`` dependency. Reads native product and reference
    transcript properties without rerunning annotation or translation. An exact common source reference scopes the
    result; missing or mixed references remain unscoped and require explicit
    compatibility when merged with a reference proteome. Unknown taxa remain
    unknown. Unsupported nonportable extension metadata raises an error.
    """
    try:
        from orfdata import Proteome
    except ImportError as error:
        raise ImportError("to_proteome requires the optional dependency: pip install 'varcode[orf]'") from error
    defaults = {name: getattr(effects, name, None) for name in ("annotator", "annotator_version")}
    export = _Export(sample_id, run_id, defaults)
    collection = None
    if isinstance(effects, EffectCollection):
        metadata = _json(dict(defaults, sources=sorted(effects.sources), annotated_at=effects.annotated_at))
        collection = export.record("collection", _identifier("EffectCollection", metadata), payload=metadata)
    results = ((effects,) if isinstance(effects, (MutationEffect, EffectCandidate, RealizedEffectCandidate,
                                                ClassifiedOutcome, ORFAnnotation, ORFOutcome)) else effects)
    for value in results:
        event_key = export.add(value)
        if collection is not None:
            export.link(collection.key, event_key, "contains")
    for annotation in orf_annotations:
        if not isinstance(annotation, ORFAnnotation):
            raise TypeError("orf_annotations must contain ORFAnnotation objects")
        export.add(annotation)
    reference = next(iter(export.references.values())) if len(export.references) == 1 and not export.unknown_reference else None
    return Proteome.from_proteins(export.proteins, reference=reference,
                                  compatibility="strict" if reference else "mixed",
                                  records=export.records, links=export.links)
