"""Portable transcript snapshots for the versioned extra-ORF effect graph.

Reconstruct ORFdata records through their public constructors and bind them to
the archived RNA. No catalogue, annotation database or chromosome DNA is read.
The snapshot describes this RNA, not a queryable reference genome.
"""

from dataclasses import dataclass, fields

from serializable import DataclassSerializable, Serializable


@dataclass(frozen=True)
class _ReferenceGenome(DataclassSerializable):
    metadata: dict

    def __getattr__(self, name):
        metadata = object.__getattribute__(self, "metadata")
        if name in metadata:
            return metadata[name]
        raise AttributeError(name)

    def __hash__(self):
        return hash(tuple(sorted(self.metadata.items())))

    fasta = None

    def transcripts_at_locus(self, *args, **kwargs):
        raise ValueError("An ORF archive contains transcript snapshots, not a genome annotation database")

    @classmethod
    def capture(cls, genome):
        names = ("reference_name", "annotation_name", "annotation_version", "release",
                 "taxon_id", "assembly_accession", "strain", "strain_taxon_id", "sequence_sha256")
        metadata = {name: getattr(genome, name, None) for name in names}
        species = getattr(genome, "species", None)
        metadata["species"] = getattr(species, "latin_name", species)
        return cls(metadata)


@dataclass(frozen=True)
class _Feature(DataclassSerializable):
    start: int
    end: int
    id: object = None
    name: object = None
    contig: object = None
    strand: object = None

    @property
    def exon_id(self):
        return self.id

    @property
    def length(self):
        return self.end - self.start + 1


class _ReferenceTranscript(Serializable):
    def __init__(self, data, genome, exons, gene):
        self.data = data
        self.genome = genome
        self.exons = exons
        self.gene = gene

    def __getattr__(self, name):
        data = object.__getattribute__(self, "data")
        if name in data:
            return data[name]
        raise AttributeError(name)

    def __hash__(self):
        return hash((self.id, self.genome, self.sequence))

    def spliced_offset(self, position):
        offset = 0
        for exon in sorted(self.exons, key=lambda e: e.start, reverse=self.on_backward_strand):
            if exon.start <= position <= exon.end:
                return offset + (exon.end - position if self.on_backward_strand else position - exon.start)
            offset += exon.length
        raise ValueError("Position is outside the archived transcript's exons")

    @classmethod
    def capture(cls, transcript):
        names = (
            "id", "name", "transcript_id", "transcript_name", "version", "versioned_id",
            "contig", "strand", "start", "end", "biotype", "topology", "sequence",
            "is_protein_coding", "complete", "gene_id", "gene_name",
            "protein_id", "protein_sequence",
            "contains_start_codon", "contains_stop_codon", "on_backward_strand",
            "on_forward_strand", "translation_table", "start_codon")
        data = {name: getattr(transcript, name, None) for name in names}
        # Incomplete fusion partners can lack either annotated codon. Their
        # PyEnsembl position accessors raise, so only read defined features.
        for boundary in ("start", "stop"):
            if data["contains_%s_codon" % boundary]:
                for suffix in ("positions", "spliced_offsets"):
                    name = "%s_codon_%s" % (boundary, suffix)
                    data[name] = getattr(transcript, name)
        if data["complete"]:
            for name in ("coding_sequence", "coding_sequence_position_ranges",
                         "five_prime_utr_sequence", "three_prime_utr_sequence"):
                data[name] = getattr(transcript, name, None)
        data["transcript_id"] = transcript.id
        data["transcript_name"] = transcript.name
        data["topology"] = getattr(transcript, "topology", "linear")
        gene = transcript.gene
        stored_gene = None if gene is None else _Feature(
            getattr(gene, "start", transcript.start), getattr(gene, "end", transcript.end),
            gene.id, gene.name, transcript.contig, transcript.strand)
        if gene is not None:
            data["gene_id"], data["gene_name"] = gene.id, gene.name
        exons = tuple(_Feature(
            exon.start, exon.end, getattr(exon, "id", getattr(exon, "exon_id", None)),
            getattr(exon, "name", None), transcript.contig, transcript.strand)
            for exon in transcript.exons)
        return cls(data, _ReferenceGenome.capture(transcript.genome), exons, stored_gene)


def orf_graph_types():
    from orfdata import catalog, evidence_summary, initiation_context, model, reference
    from orfdata import liftover
    from .orf import ORFAnnotation, ORFOutcome, ORFProduct
    from .mutant_transcript import MutantTranscript, ReferenceSegment, TranscriptEdit, _AssembledAllele

    # Public data records only; never serialize a catalogue, provider or cache.
    names = {
        "ORFRecord", "UORFRecord", "PeptideEvidence", "ORFEvidence", "ORFModel",
        "ReferenceContext", "ScopedTarget", "InitiationContext", "EvidenceSummary", "Liftover",
    }
    records = {getattr(module, name) for module in (
        catalog, evidence_summary, initiation_context, model, reference, liftover)
        for name in names if hasattr(module, name)}
    return tuple(records) + (ORFAnnotation, ORFOutcome, ORFProduct,
                            MutantTranscript, ReferenceSegment, TranscriptEdit, _AssembledAllele)


def record_fields(value):
    """Constructor fields, preserving nested records instead of export summaries."""
    return {field.name: getattr(value, field.name) for field in fields(value) if field.init}
