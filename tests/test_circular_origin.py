"""Transcripts joined across the origin of a circular contig.

Hepatitis B polymerase crosses the origin of the circular HBV genome and is
annotated as GenBank ``join(2307..3215,1..1623)``: two exons in transcript
order that meet at base 1. These tests place a small coding sequence the
same way on a 60-base circular contig and compare varcode with an
independent translation of the mutated circle.
"""

from types import SimpleNamespace

import pytest
from pyensembl import Genome as AnnotationGenome

from varcode import Genome, Variant
from varcode.circular_origin import origin_join_index

BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {
    a + b + c: AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)}

ORF = "ATG" + "GCT" * 5 + "TAA"  # MAAAAA*, 12 bases before the origin
PROTEIN = "MAAAAA"
LENGTH = 60

# Each strand: the exons in transcript order, the CDS, and the codons.
LAYOUTS = {
    # Sense reads 49..60, then 1..9.
    "+": dict(exons=((49, 60), (1, 9)), cds=((49, 60), (1, 6)),
              start_codon=(49, 51), stop_codon=(7, 9)),
    # Sense reads 12..1, then 60..52.
    "-": dict(exons=((1, 12), (52, 60)), cds=((1, 12), (55, 60)),
              start_codon=(10, 12), stop_codon=(52, 54)),
}


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def contig_sequence(strand):
    bases = ["T"] * LENGTH
    if strand == "+":
        bases[48:60] = ORF[:12]
        bases[0:9] = ORF[12:]
    else:
        bases[0:12] = reverse_complement(ORF[:12])
        bases[51:60] = reverse_complement(ORF[12:])
    return "".join(bases)


def sense_of_circle(sequence, strand):
    """The coding strand read across the origin, independent of varcode."""
    if strand == "+":
        return sequence[48:60] + sequence[0:9]
    return reverse_complement(sequence[0:12]) + reverse_complement(sequence[51:60])


def translate(coding):
    protein = ""
    for offset in range(0, len(coding) - 2, 3):
        protein += CODONS[coding[offset:offset + 3]]
        if protein.endswith("*"):
            break
    return protein


def gtf_lines(strand):
    layout = LAYOUTS[strand]
    shared = ('gene_id "G"; gene_name "P"; gene_biotype "protein_coding"; '
              'transcript_id "T"; transcript_name "P-201"; '
              'transcript_biotype "protein_coding";')

    def row(feature, start, end, extra=""):
        return "\t".join((
            "C", "test", feature, str(start), str(end), ".", strand, ".",
            shared + extra))

    rows = [
        row("gene", 1, LENGTH).replace(' transcript_id "T"; transcript_name "P-201"; '
                                       'transcript_biotype "protein_coding";', ""),
        row("transcript", 1, LENGTH)]
    for number, (start, end) in enumerate(layout["exons"], 1):
        rows.append(row("exon", start, end, ' exon_number "%d"; exon_id "E%d";' % (
            number, number)))
    for number, (start, end) in enumerate(layout["cds"], 1):
        rows.append(row("CDS", start, end, ' exon_number "%d"; protein_id "P1";' % number))
    rows.append(row("start_codon", *layout["start_codon"], ' exon_number "1";'))
    rows.append(row("stop_codon", *layout["stop_codon"], ' exon_number "2";'))
    return rows


@pytest.fixture(scope="module", params=["+", "-"])
def circle(request, tmp_path_factory):
    strand = request.param
    directory = tmp_path_factory.mktemp("circle")
    (directory / "annotation.gtf").write_text("\n".join(gtf_lines(strand)) + "\n")
    (directory / "transcripts.fa").write_text(">T\n%s\n" % ORF)
    (directory / "proteins.fa").write_text(">P1\n%s\n" % PROTEIN)
    annotation = AnnotationGenome(
        reference_name="circle", annotation_name="circle",
        gtf_path_or_url=str(directory / "annotation.gtf"),
        transcript_fasta_paths_or_urls=[str(directory / "transcripts.fa")],
        protein_fasta_paths_or_urls=[str(directory / "proteins.fa")],
        cache_directory_path=str(directory))
    annotation.index()
    sequence = contig_sequence(strand)
    yield SimpleNamespace(
        strand=strand, sequence=sequence,
        genome=Genome(annotation, fasta={"C": sequence}, verify=False))


def snvs(strand):
    """Every single-base change to each base of the transcript."""
    for start, end in LAYOUTS[strand]["exons"]:
        for position in range(start, end + 1):
            ref = contig_sequence(strand)[position - 1]
            for alt in "ACGT":
                if alt != ref:
                    yield position, ref, alt


def expected_effect(circle, position, ref, alt):
    sequence = circle.sequence
    mutated = sequence[:position - 1] + alt + sequence[position:]
    protein = translate(sense_of_circle(mutated, circle.strand))
    if protein == PROTEIN + "*":
        return {"Silent"}, None
    if not protein.startswith("M"):
        return {"StartLoss", "AlternateStartCodon"}, None
    if not protein.endswith("*"):
        return {"StopLoss"}, None
    index = next(i for i, (a, b) in enumerate(zip(PROTEIN + "*", protein)) if a != b)
    return ({"Substitution", "PrematureStop"},
            "p.%s%d%s" % (PROTEIN[index], index + 1, protein[index]))


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_every_snv_matches_translation_of_the_circle(circle, annotator):
    wrong = []
    for position, ref, alt in snvs(circle.strand):
        classes, description = expected_effect(circle, position, ref, alt)
        effect = Variant("C", position, ref, alt, circle.genome).effects(
            annotator=annotator).top_priority_effect()
        name = type(effect).__name__
        if name not in classes or (
                description is not None and effect.short_description != description):
            wrong.append((position, ref, alt, sorted(classes), description,
                          name, effect.short_description))
    assert not wrong


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_bases_at_the_origin_are_not_splice_sites(circle, annotator):
    exons = LAYOUTS[circle.strand]["exons"]
    # The last base before the origin and the first after it.
    if circle.strand == "+":
        positions = (exons[0][1] - 1, exons[0][1], exons[1][0], exons[1][0] + 1)
    else:
        positions = (exons[0][0] + 1, exons[0][0], exons[1][1], exons[1][1] - 1)
    for position in positions:
        ref = circle.sequence[position - 1]
        alt = "C" if ref != "C" else "A"
        effect = Variant("C", position, ref, alt, circle.genome).effects(
            annotator=annotator).top_priority_effect()
        assert type(effect).__name__ in ("Silent", "Substitution"), (position, effect)
        # The default annotator records disrupted donors and acceptors here.
        assert getattr(effect, "splice_signal", None) is None, (position, effect)


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_indels_after_the_origin(circle, annotator):
    # Codon 6 (GCT) lies after the origin: 4..6 forward, 57..55 reverse.
    if circle.strand == "+":
        in_frame, frameshift = (3, "TGCT", "T"), (4, "GC", "G")
    else:
        in_frame, frameshift = (54, "AAGC", "A"), (55, "AG", "A")
    for position, ref, _ in (in_frame, frameshift):
        assert circle.sequence[position - 1:position - 1 + len(ref)] == ref
    effect = Variant("C", *in_frame, circle.genome).effects(
        annotator=annotator).top_priority_effect()
    assert type(effect).__name__ == "Deletion"
    assert effect.short_description in ("p.A5del", "p.A6del")
    effect = Variant("C", *frameshift, circle.genome).effects(
        annotator=annotator).top_priority_effect()
    assert type(effect).__name__ in ("FrameShift", "FrameShiftTruncation")


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_stretch_between_the_ends_is_outside_the_transcript(circle, annotator):
    exons = LAYOUTS[circle.strand]["exons"]
    if circle.strand == "+":
        stretch = range(exons[1][1] + 1, exons[0][0])
    else:
        stretch = range(exons[0][1] + 1, exons[1][0])
    for position in (stretch[0], stretch[1], stretch[-2], stretch[-1]):
        ref = circle.sequence[position - 1]
        variant = Variant("C", position, ref, "C" if ref != "C" else "A", circle.genome)
        assert variant.transcripts == []
        assert [type(effect).__name__ for effect in variant.effects(
            annotator=annotator)] == ["Intergenic"]


def test_origin_join_index_finds_the_one_step_back(circle):
    transcript = circle.genome.transcript_by_id("T")
    assert origin_join_index(transcript) == 0


@pytest.mark.parametrize("strand,exons", [
    ("+", ((40, 50), (5, 9))),         # steps back, but not onto base 1
    ("+", ((40, 45), (1, 5), (20, 25), (1, 3))),  # steps back twice
    ("-", ((40, 50), (55, 60))),      # reverse strand, not from base 1
])
def test_other_steps_back_are_refused(strand, exons):
    transcript = SimpleNamespace(
        id="T", strand=strand, on_backward_strand=strand == "-",
        exons=[SimpleNamespace(start=start, end=end) for start, end in exons])
    with pytest.raises(ValueError, match="step back"):
        origin_join_index(transcript)


def test_linear_transcripts_have_no_origin_join():
    for strand, exons in (("+", ((1, 5), (9, 12))), ("-", ((9, 12), (1, 5)))):
        transcript = SimpleNamespace(
            id="T", strand=strand, on_backward_strand=strand == "-",
            exons=[SimpleNamespace(start=start, end=end) for start, end in exons])
        assert origin_join_index(transcript) is None
