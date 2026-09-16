# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Small adapter for Exacto's observed transcript structures and DNA/RNA links.

This is deliberately not a second variant caller or ORF finder. Existing DNA
variants anchor the import, and Exacto's ordered sequence/structure is retained.
"""

from collections import defaultdict
from contextlib import nullcontext
import csv
import gzip
import os

from .rna_evidence import RNAEvidence, make_fusion_outcome


def _rows(path, required):
    if hasattr(path, "read"):
        stream = nullcontext(path)
    else:
        name = os.fspath(path)
        opener = gzip.open if name.endswith(".gz") else open
        stream = opener(name, "rt", newline="")
    with stream as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = set(required) - set(fields)
        if missing or len(fields) != len(set(fields)):
            raise ValueError("Invalid Exacto columns; missing: %s" % sorted(missing))
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Ragged Exacto row at line %d" % reader.line_num)
            yield row


def _model_key(row):
    model_id = row["transcript_model_id"]
    if not model_id:
        raise ValueError("Empty transcript_model_id")
    refs = tuple(sorted(set(filter(None, row["reference_transcript_ids"].split(",")))))
    return model_id, refs


def _transcript(genome, transcript_id):
    if not transcript_id:
        return None
    # Ensembl stores unversioned IDs; keep the producer's full ID in evidence.
    try:
        return genome.transcript_by_id(transcript_id)
    except ValueError:
        if "." not in transcript_id:
            raise
        return genome.transcript_by_id(transcript_id.rsplit(".", 1)[0])


def _sequence(rows):
    rows = sorted(rows, key=lambda row: int(row["index"]))
    if [int(row["index"]) for row in rows] != list(range(len(rows))):
        raise ValueError("Structure indices must be unique and contiguous from zero")
    bases = []
    previous_end = None
    for row in rows:
        if row["type"] not in {"base", "event"}:
            raise ValueError("Unknown Exacto structure type: %r" % row["type"])
        if row["type"] == "event":
            if row["sequence"]:
                raise ValueError("Sequence-bearing Exacto events are not supported")
            continue
        start, end = int(row["read_start"]), int(row["read_end"])
        if start < 0 or end < start or len(row["sequence"]) != end - start + 1:
            raise ValueError("Base sequence length disagrees with inclusive read coordinates")
        if previous_end is not None and start != previous_end + 1:
            raise ValueError("Gapped or overlapping observed sequence; cannot assemble")
        previous_end = end
        bases.append(row)
    if not bases:
        raise ValueError("No observed bases for transcript model")
    return rows, bases, "".join(row["sequence"] for row in bases)


def _validate_path(bases, genome):
    """Reject multiway/back-spliced paths instead of flattening them to two ends."""
    runs = []
    previous_end = None
    for row in bases:
        chrom, strand = row["chromosome_1"], row["strand_1"]
        if (chrom != row["chromosome_2"] or strand != row["strand_2"]
                or strand not in {"+", "-"}):
            raise ValueError("A base row must describe one contig and strand")
        tid = row["transcript_id_1"]
        if tid != row["transcript_id_2"]:
            raise ValueError("A base row must have one transcript annotation")
        transcript = _transcript(genome, tid)
        if transcript is not None and chrom.removeprefix("chr") != transcript.contig.removeprefix("chr"):
            raise ValueError("Transcript annotation and base-row contig disagree")
        owner = (chrom, strand, transcript.gene_id if transcript is not None else None)
        lo, hi = sorted([int(row["position_1"]), int(row["position_2"])])
        if lo < 1:
            raise ValueError("Genomic positions must be positive")
        direction = 1 if strand == "+" else -1
        start, end = (lo, hi) if strand == "+" else (hi, lo)
        if not runs or owner != runs[-1]:
            runs.append(owner)
            if len(runs) > 2:
                raise ValueError("More than two loci/orientations in an observed fusion path")
        elif previous_end is not None and direction * (start - previous_end) <= 0 and row["kind"] != "insertion":
            raise ValueError("Overlapping or back-spliced genomic path is not a linear fusion")
        previous_end = end


def load_exacto_fusions(structures_path, integrated_path, *, variants_by_id,
                       cds_starts=None):
    """Import selected SV-linked RNA models as an RNAEvidence resolver.

    Parameters
    ----------
    structures_path, integrated_path : path or text stream
        Exacto transcript-structures and integrated-variants TSVs (optionally
        gzip-compressed). Structure ``index`` order is transcript order;
        sequences are already oriented and are NOT reverse-complemented again.
    variants_by_id : mapping
        Exacto DNA call IDs to existing StructuralVariant objects, normally
        loaded from VCF. Only these IDs are selected from the integration table.
        The explicit join avoids inventing a DNA breakpoint from an RNA splice.
    cds_starts : mapping or None
        Optional, explicitly chosen zero-based ORF starts in assembled sequence,
        keyed by ``(transcript_model_id, tuple(sorted(reference_transcript_ids)))``.
        No ORF is chosen automatically. See ``make_fusion_outcome``.

    Returns
    -------
    RNAEvidence
        ``.candidates`` retains every selected variant/model/reference-ID group.
        Pass the result to ``effects(rna_resolver=...)``. Original structure and
        integration rows, including RNA and DNA call IDs, remain in evidence.

    Notes
    -----
    Supports linear two-locus models with an annotated sense 5' anchor. A
    missing or antisense 3' partner remains TranslocationToIntergenic, never a
    guessed coding fusion. Incomplete sequence, unknown transcript IDs, circular
    paths and multi-gene (>2) models raise rather than silently losing structure.
    This does not import all Exacto variant types or primary-structure tables.
    Read/model completeness and read support are not inferred from row counts.
    """
    variants = {str(key): value for key, value in variants_by_id.items()}
    if len(variants) != len(variants_by_id):
        raise ValueError("DNA call IDs collide after conversion to strings")
    if any(not getattr(v, "is_structural", False) for v in variants.values()):
        raise ValueError("variants_by_id must contain only StructuralVariant objects")
    links = defaultdict(list)
    for row in _rows(integrated_path, ["transcript_model_id", "reference_transcript_ids",
                                       "rna_variant_call_id", "dna_variant_call_id"]):
        if row["dna_variant_call_id"] in variants:
            key = (*_model_key(row), row["dna_variant_call_id"])
            if row not in links[key]:
                links[key].append(row)
    selected = {key[:2] for key in links}
    models = defaultdict(list)
    fields = ["transcript_model_id", "reference_transcript_ids", "index",
              "read_start", "read_end", "sequence", "type", "kind", "context",
              "chromosome_1", "position_1", "strand_1", "chromosome_2",
              "position_2", "strand_2", "transcript_id_1", "transcript_id_2"]
    for row in _rows(structures_path, fields):
        key = _model_key(row)
        if key in selected:
            models[key].append(row)
    if selected - set(models):
        raise ValueError("Missing structures for linked transcript models: %s" %
                         sorted(selected - set(models)))
    candidates = []
    for (model_id, refs, dna_id), integration_rows in sorted(links.items()):
        variant = variants[dna_id]
        rows, bases, sequence = _sequence(models[(model_id, refs)])
        _validate_path(bases, variant.genome)
        if any("circular" in row["kind"].lower() or "circular" in row["context"].lower()
               for row in rows):
            raise ValueError("Circular RNA requires a separate model, not a linear fusion")
        first, last = bases[0], bases[-1]
        transcript = _transcript(variant.genome, first["transcript_id_1"])
        if transcript is None or first["strand_1"] != transcript.strand:
            raise ValueError("An annotated sense 5-prime transcript anchor is required")
        partner = _transcript(variant.genome, last["transcript_id_2"])
        transcript_ids = {row[k] for row in bases
                          for k in ["transcript_id_1", "transcript_id_2"] if row[k]}
        genes = {_transcript(variant.genome, tid).gene_id for tid in transcript_ids}
        if len(genes) > 2:
            raise ValueError("More than two genes in an observed fusion path")
        partner_status = "sense" if partner is not None else "unannotated"
        if partner is not None and last["strand_2"] != partner.strand:
            partner_status, partner = "antisense", None
        if partner is not None and partner.gene_id == transcript.gene_id:
            raise ValueError("Selected model does not identify a two-gene fusion")
        evidence = dict(
            dna_variant_call_id=dna_id,
            rna_variant_call_ids=sorted({r["rna_variant_call_id"] for r in integration_rows}),
            reference_transcript_ids=list(refs), partner_status=partner_status,
            exacto_structure=rows, exacto_integration=integration_rows,
            sequence_status="observed_model_completeness_unknown")
        candidates.append(make_fusion_outcome(
            variant, transcript, sequence=sequence, transcript_model_id=model_id,
            partner_transcript=partner, source="exacto",
            cds_start=(cds_starts or {}).get((model_id, refs)),
            extra_evidence=evidence))
    return RNAEvidence(candidates)
