"""Transcripts joined across the origin of a circular contig.

Hepatitis B polymerase crosses the origin of the circular HBV genome and is
annotated as GenBank ``join(2307..3215,1..1623)``: exons in transcript order
that meet at base 1. These tests place a three-exon transcript the same way
on a 400-base circle, on each strand: exon 1 and exon 2 are spliced, and
exon 2 meets exon 3 at the origin. The stretch between the transcript's 3'
and 5' ends is wider than the layout flanks. Varcode is compared with an
independent reading of the mutated circle, and on the real HBV genome.
"""

import random
from pathlib import Path
from types import SimpleNamespace

import pytest
from pyensembl import Genome as AnnotationGenome

from varcode import (
    Genome,
    GermlineContext,
    StructuralVariant,
    Variant,
    VariantCollection,
)
from varcode.circular_origin import TranscriptFrame, in_transcript_order
from varcode.germline import Completeness
from varcode.genomic_layout import GenomicLayout
from varcode.splice_graph import SpliceSiteKey, SpliceSiteStatus, splice_axes
from varcode.transcript_layout import build_exon_runs

BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {
    a + b + c: AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)}
SENSE_CODONS = sorted(codon for codon, aa in CODONS.items() if aa != "*")
ANNOTATORS = ["transcript_model", "fast"]

LENGTH = 400
# Forward strand, transcript order: exon 1 and exon 2 are spliced across
# the intron 331..339; exon 2 ends at the contig's last base and exon 3
# continues from base 1.
FORWARD_EXONS = ((300, 330), (340, 400), (1, 60))
UTR5, CODING = 6, 135  # 44 residues and the stop codon


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate(coding):
    protein = ""
    for offset in range(0, len(coding) - 2, 3):
        protein += CODONS[coding[offset:offset + 3]]
        if protein.endswith("*"):
            break
    return protein


def mirror(position):
    return LENGTH + 1 - position


def build_circle(strand):
    rng = random.Random(686)
    cdna_length = sum(end - start + 1 for start, end in FORWARD_EXONS)
    coding = "ATG" + "".join(rng.choice(SENSE_CODONS) for _ in range(43)) + "TAA"
    cdna = ("".join(rng.choice("ACGT") for _ in range(UTR5)) + coding
            + "".join(rng.choice("ACGT") for _ in range(cdna_length - UTR5 - CODING)))
    bases = [rng.choice("ACGT") for _ in range(LENGTH)]
    offset = 0
    for start, end in FORWARD_EXONS:
        for position in range(start, end + 1):
            bases[position - 1] = cdna[offset]
            offset += 1
    bases[330:332] = "GT"   # donor at 331..332
    bases[337:339] = "AG"   # acceptor at 338..339
    forward = "".join(bases)
    if strand == "+":
        return forward, FORWARD_EXONS, cdna
    exons = tuple((mirror(end), mirror(start)) for start, end in FORWARD_EXONS)
    return reverse_complement(forward), exons, cdna


def cdna_blocks(exons, first, last):
    """Forward-strand blocks of transcript bases ``first..last`` (zero-based)."""
    blocks, offset = [], 0
    for start, end in exons:
        length = end - start + 1
        low, high = max(first, offset), min(last, offset + length - 1)
        if low <= high:
            blocks.append((start + low - offset, start + high - offset))
        offset += length
    return blocks


def gtf_rows(strand, exons, gene="G", transcript="T", protein="P1",
             utr5=UTR5, coding=CODING, contig="C"):
    forward_exons = exons if strand == "+" else tuple(
        (mirror(end), mirror(start)) for start, end in exons)
    shared = ('gene_id "%s"; gene_name "%s"; gene_biotype "protein_coding"; '
              'transcript_id "%s"; transcript_name "%s-201"; '
              'transcript_biotype "protein_coding";' % (gene, gene, transcript, gene))

    def row(feature, start, end, frame=".", extra=""):
        return "\t".join((contig, "test", feature, str(start), str(end), ".", strand,
                          str(frame), shared + extra))

    def blocks(first, last):
        found = cdna_blocks(forward_exons, first, last)
        if strand == "+":
            return found
        return [(mirror(end), mirror(start)) for start, end in found]

    rows = [row("transcript", min(s for s, _ in exons), max(e for _, e in exons))]
    for number, (start, end) in enumerate(exons, 1):
        rows.append(row("exon", start, end, extra=' exon_number "%d"; exon_id "%s.E%d";' % (
            number, transcript, number)))
    done = 0
    for start, end in blocks(utr5, utr5 + coding - 4):
        rows.append(row("CDS", start, end, frame=(3 - done % 3) % 3,
                        extra=' protein_id "%s";' % protein))
        done += end - start + 1
    for start, end in blocks(utr5, utr5 + 2):
        rows.append(row("start_codon", start, end, frame=0))
    for start, end in blocks(utr5 + coding - 3, utr5 + coding - 1):
        rows.append(row("stop_codon", start, end, frame=0))
    return rows


def annotation(directory, rows, transcripts, proteins):
    (directory / "annotation.gtf").write_text("\n".join(rows) + "\n")
    (directory / "transcripts.fa").write_text(
        "".join(">%s\n%s\n" % pair for pair in transcripts.items()))
    (directory / "proteins.fa").write_text(
        "".join(">%s\n%s\n" % pair for pair in proteins.items()))
    genome = AnnotationGenome(
        reference_name="circle", annotation_name=directory.name,
        gtf_path_or_url=str(directory / "annotation.gtf"),
        transcript_fasta_paths_or_urls=[str(directory / "transcripts.fa")],
        protein_fasta_paths_or_urls=[str(directory / "proteins.fa")],
        cache_directory_path=str(directory))
    genome.index()
    return genome


PARTNER = "ATGGCC" + "GGT" * 20 + "TGA"  # a linear partner gene on contig D, 21..89


@pytest.fixture(scope="module", params=["+", "-"])
def circle(request, tmp_path_factory):
    strand = request.param
    sequence, exons, cdna = build_circle(strand)
    protein = translate(cdna[UTR5:])[:-1]
    partner_contig = "T" * 20 + PARTNER + "T" * 20
    partner_rows = [
        "\t".join(("D", "test", feature, str(start), str(end), ".", "+", ".",
                   'gene_id "H"; gene_name "H"; gene_biotype "protein_coding"; '
                   'transcript_id "TH"; transcript_name "H-201"; '
                   'transcript_biotype "protein_coding";' + extra))
        for feature, start, end, extra in (
            ("transcript", 21, 89, ""),
            ("exon", 21, 89, ' exon_number "1"; exon_id "TH.E1";'),
            ("CDS", 21, 86, ' protein_id "PH";'),
            ("start_codon", 21, 23, ""),
            ("stop_codon", 87, 89, ""))]
    genome = annotation(
        tmp_path_factory.mktemp("circle"), gtf_rows(strand, exons) + partner_rows,
        {"T": cdna, "TH": PARTNER}, {"P1": protein, "PH": translate(PARTNER)[:-1]})
    yield SimpleNamespace(
        strand=strand, sequence=sequence, exons=exons, cdna=cdna, protein=protein,
        genome=Genome(genome, fasta={"C": sequence, "D": partner_contig}, verify=False))


def to_strand(circle, position):
    """A forward-strand design coordinate on this circle's strand."""
    return position if circle.strand == "+" else mirror(position)


def read_circle(circle, sequence):
    """The transcript's bases on a mutated circle of unchanged length."""
    pieces = [sequence[start - 1:end] for start, end in circle.exons]
    if circle.strand == "-":
        pieces = [reverse_complement(piece) for piece in pieces]
    return "".join(pieces)


def cdna_offset(circle, position):
    offset = 0
    for start, end in circle.exons:
        if start <= position <= end:
            return offset + (position - start if circle.strand == "+" else end - position)
        offset += end - start + 1
    return None


def expected_snv(circle, position, alt):
    offset = cdna_offset(circle, position)
    if offset < UTR5:
        return {"FivePrimeUTR"}, None
    if offset >= UTR5 + CODING:
        return {"ThreePrimeUTR"}, None
    mutated = circle.sequence[:position - 1] + alt + circle.sequence[position:]
    protein = translate(read_circle(circle, mutated)[UTR5:])
    reference = circle.protein + "*"
    if protein == reference:
        return {"Silent"}, None
    if not protein.startswith("M"):
        return {"StartLoss", "AlternateStartCodon"}, None
    if not protein.endswith("*"):
        return {"StopLoss"}, None
    index = next(i for i, (a, b) in enumerate(zip(reference, protein)) if a != b)
    return ({"Substitution", "PrematureStop"},
            "p.%s%d%s" % (reference[index], index + 1, protein[index]))


def on_transcript(variant, annotator, **kwargs):
    """The effect on transcript T; for an SV with alternatives (such as
    cryptic exons), its most likely one."""
    effects = [effect for effect in variant.effects(annotator=annotator, **kwargs)
               if getattr(getattr(effect, "transcript", None), "id", None) == "T"]
    assert len(effects) == 1, effects
    effect = effects[0]
    if type(effect).__name__ == "StructuralVariantEffect":
        return effect.most_likely_effect
    return effect


def splice_free_positions(circle):
    """Transcript bases at least three from the spliced junction."""
    junction = {to_strand(circle, p) for p in (328, 329, 330, 340, 341, 342)}
    return [p for start, end in circle.exons for p in range(start, end + 1)
            if p not in junction]


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_every_snv_matches_translation_of_the_circle(circle, annotator):
    wrong = []
    for position in splice_free_positions(circle):
        ref = circle.sequence[position - 1]
        for alt in "ACGT".replace(ref, ""):
            classes, description = expected_snv(circle, position, alt)
            effect = on_transcript(Variant("C", position, ref, alt, circle.genome), annotator)
            name = type(effect).__name__
            if name not in classes or (
                    description is not None and effect.short_description != description):
                wrong.append((position, ref, alt, sorted(classes), description,
                              name, effect.short_description))
    assert not wrong


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_the_origin_is_not_a_splice_site(circle, annotator):
    for design in (398, 399, 400, 1, 2, 3):
        position = to_strand(circle, design)
        ref = circle.sequence[position - 1]
        alt = "ACGT".replace(ref, "")[0]
        effect = on_transcript(Variant("C", position, ref, alt, circle.genome), annotator)
        assert type(effect).__name__ in (
            "Silent", "Substitution", "PrematureStop"), (position, effect)
        assert getattr(effect, "splice_signal", None) is None, (position, effect)


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_exons_meeting_at_the_origin_are_skipped_together(circle, annotator):
    # The acceptor of exon 2. Exon 2 and exon 3 are one exon on the circle,
    # and it ends the transcript, so no exon is skipped.
    position = to_strand(circle, 339)
    ref = circle.sequence[position - 1]
    effect = on_transcript(
        Variant("C", position, ref, "C" if ref != "C" else "T", circle.genome), annotator)
    candidates = getattr(effect, "candidates", None) or ()
    assert candidates
    assert not any(type(candidate.effect).__name__ == "ExonSkipping"
                   for candidate in candidates)
    assert not any(
        option.mechanism == "exon_skip"
        for candidate in candidates
        for hypothesis in (getattr(candidate, "hypotheses", None) or ())
        for _, option in getattr(getattr(hypothesis, "splice_plan", None), "choices", ()))


def test_splice_axes_skip_both_exons_meeting_at_an_origin():
    exons = [SimpleNamespace(start=s, end=e)
             for s, e in ((100, 150), (300, 400), (1, 50), (80, 90))]
    transcript = SimpleNamespace(
        contig="C", strand="+", on_backward_strand=False, exons=exons, start=1, end=400)
    run_keys = tuple((number, 1) for number in range(1, 5))
    status = SpliceSiteStatus(
        SpliceSiteKey("acceptor", 2, 1), (2, 1), "lost", "acceptor motif changed", ())
    axis, = splice_axes(transcript, (status,), run_keys)
    skip, = [option for option in axis.options if option.mechanism == "exon_skip"]
    assert skip.skipped_runs == ((2, 1), (3, 1))


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_indels_after_the_origin(circle, annotator):
    # Exon 3 holds codons 30..44; delete design bases 16..18, then one base.
    anchor = to_strand(circle, 15) if circle.strand == "+" else to_strand(circle, 19)
    ref = circle.sequence[anchor - 1:anchor + 3]
    effect = on_transcript(Variant("C", anchor, ref, ref[0], circle.genome), annotator)
    assert type(effect).__name__ in ("Deletion", "ComplexSubstitution"), effect
    effect = on_transcript(Variant("C", anchor, ref[:2], ref[0], circle.genome), annotator)
    assert type(effect).__name__ in ("FrameShift", "FrameShiftTruncation"), effect


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_stretch_between_the_ends_is_outside_the_transcript(circle, annotator):
    for design in (61, 62, 120, 200, 298, 299):
        position = to_strand(circle, design)
        ref = circle.sequence[position - 1]
        variant = Variant("C", position, ref, "ACGT".replace(ref, "")[0], circle.genome)
        assert variant.transcripts == []
        assert variant.gene_names == [] and variant.gene_ids == []
        assert [type(e).__name__ for e in variant.effects(annotator=annotator)] == [
            "Intergenic"]


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_allele_reaching_past_the_flank(circle, annotator):
    # Deletes the stop codon and runs 70 bases into the stretch between ends.
    low, high = sorted((to_strand(circle, 40), to_strand(circle, 130)))
    ref = circle.sequence[low - 2:high]
    effect = on_transcript(Variant("C", low - 1, ref, ref[0], circle.genome), annotator)
    assert type(effect).__name__ not in ("Failure", "Unresolved"), effect


def test_germline_between_the_ends_is_outside_the_layout(circle):
    germline_position = to_strand(circle, 200)
    germline_ref = circle.sequence[germline_position - 1]
    germline = Variant("C", germline_position, germline_ref,
                       "ACGT".replace(germline_ref, "")[0], circle.genome)
    context = GermlineContext.from_variants(
        VariantCollection([germline]), completeness=Completeness.COMPLETE)
    position = to_strand(circle, 15)
    ref = circle.sequence[position - 1]
    alt = "ACGT".replace(ref, "")[0]
    classes, description = expected_snv(circle, position, alt)
    effect = on_transcript(Variant("C", position, ref, alt, circle.genome),
                           "transcript_model", germline=context)
    assert type(effect).__name__ in classes, effect


def sv(circle, kind, low, high, **kwargs):
    start, end = sorted((to_strand(circle, low), to_strand(circle, high)))
    return StructuralVariant("C", start, kind, end=end, genome=circle.genome, **kwargs)


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_deletions_on_the_circle(circle, annotator):
    # In frame inside exon 3: design 17..22 holds codons 35 and 36.
    effect = on_transcript(sv(circle, "DEL", 17, 22), annotator)
    assert type(effect).__name__ == "Deletion", effect
    # From the stretch between the ends into exon 1, taking the start codon.
    effect = on_transcript(sv(circle, "DEL", 250, 310), annotator)
    assert type(effect).__name__ == "StartLoss", effect
    # Through the stretch between the ends, taking both ends.
    effect = on_transcript(sv(circle, "DEL", 50, 310), annotator)
    assert type(effect).__name__ == "StartLoss", effect
    # Entirely between the ends.
    effects = sv(circle, "DEL", 100, 200).effects(annotator=annotator)
    assert [type(e).__name__ for e in effects] == ["Intergenic"]


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_deletion_across_the_origin(circle, annotator):
    # 398..400 and 1..3: six bases, in frame. Written past the contig's end
    # and as a breakend junction that skips the origin.
    records = [
        StructuralVariant("C", 398, "DEL", end=403, genome=circle.genome),
        StructuralVariant("C", 397, "BND", alt="N[C:4[", genome=circle.genome)]
    descriptions = set()
    for record in records:
        effect = on_transcript(record, annotator)
        assert type(effect).__name__ in ("Deletion", "ComplexSubstitution"), (
            record, effect)
        descriptions.add(effect.short_description)
    assert len(descriptions) == 1, descriptions


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_duplication_after_the_origin(circle, annotator):
    effect = on_transcript(sv(circle, "DUP", 17, 22), annotator)
    assert type(effect).__name__ == "Insertion", effect


def test_fusion_from_the_circle(circle):
    # Keep the transcript through 30 bases of exon 3 and join the partner on
    # contig D from its third codon (D:27).
    position = to_strand(circle, 30)
    # On the reverse strand the circle keeps its right side, joined to the
    # reverse complement of the partner's right side.
    alt = "N[D:27[" if circle.strand == "+" else "[D:27[N"
    record = StructuralVariant("C", position, "BND", alt=alt, genome=circle.genome)
    effects = record.effects(annotator="fast")
    fusions = [e for e in effects
               if type(e).__name__ == "GeneFusion" and e.transcript.id == "T"]
    assert fusions, list(effects)
    kept = 31 + 61 + 30  # exon 1, exon 2 and 30 bases of exon 3
    assert fusions[0].mutant_transcript.cdna_sequence.startswith(circle.cdna[:kept])


def test_breakend_between_the_ends_finds_no_partner_on_the_circle(circle):
    record = StructuralVariant(
        "D", 30, "BND", alt="N[C:%d[" % to_strand(circle, 200), genome=circle.genome)
    assert not any(
        getattr(getattr(effect, "partner_transcript", None), "id", None) == "T"
        for effect in record.effects(annotator="fast"))


def test_frame_of_an_origin_join():
    transcript = SimpleNamespace(
        contig="C", strand="+", on_backward_strand=False, start=1, end=400,
        exons=[SimpleNamespace(start=s, end=e) for s, e in FORWARD_EXONS])
    frame = TranscriptFrame(transcript)
    assert frame.contig_length == 400
    assert frame.exon_intervals == ((300, 330), (340, 400), (401, 460))
    assert frame.images(50, 310) == ((50, 310), (450, 710))
    assert frame.images(100, 200) == ((100, 200),)
    assert not frame.overlaps(61, 299) and frame.overlaps(60, 61)
    assert frame.origin_partners(2) == frame.origin_partners(3) == (2, 3)
    assert frame.origin_partners(1) == (1,)


def test_terminally_redundant_transcript_keeps_both_copies():
    # Like HBV pregenomic RNA, the 3' end runs past the 5' end on the circle.
    transcript = SimpleNamespace(
        contig="C", strand="+", on_backward_strand=False, start=1, end=400,
        exons=[SimpleNamespace(start=350, end=400), SimpleNamespace(start=1, end=370)])
    frame = TranscriptFrame(transcript)
    assert frame.exon_intervals == ((350, 400), (401, 770))
    assert frame.images(360, 360) == ((360, 360), (760, 760))
    layout = GenomicLayout.from_transcript(transcript)
    runs = build_exon_runs(transcript, layout)
    assert [run.length for run in runs] == [51, 370]


@pytest.mark.parametrize("blocks,backward,expected", [
    # Across an origin, as annotated.
    (((300, 330), (340, 400), (1, 60)), False, ((300, 330), (340, 400), (1, 60))),
    (((71, 101), (1, 61), (341, 400)), True, ((71, 101), (1, 61), (341, 400))),
    # Reverse-strand exons numbered left to right are put in transcript order.
    (((100, 116), (200, 206)), True, ((200, 206), (100, 116))),
    (((5, 9), (20, 30)), False, ((5, 9), (20, 30))),
])
def test_transcript_order(blocks, backward, expected):
    assert in_transcript_order(blocks, backward) == expected


def test_mixed_isoforms_leave_the_stretch_between_ends(tmp_path):
    sequence, exons, cdna = build_circle("+")
    protein = translate(cdna[UTR5:])[:-1]
    rows = gtf_rows("+", exons)
    # A second, linear isoform of the same gene on exons 1 and 2 only.
    linear = [
        row.replace('transcript_id "T"', 'transcript_id "T2"').replace('"T.E', '"T2.E')
        for row in rows
        if "\tCDS\t" not in row and "_codon\t" not in row and "\t1\t60\t" not in row]
    linear = [row.replace("\ttranscript\t1\t400\t", "\ttranscript\t300\t400\t")
              for row in linear]
    genome = Genome(
        annotation(tmp_path, rows + linear, {"T": cdna, "T2": cdna[:92]}, {"P1": protein}),
        fasta={"C": sequence}, verify=False)
    variant = Variant("C", 200, sequence[199], "ACGT".replace(sequence[199], "")[0], genome)
    assert variant.gene_ids == []
    assert [type(e).__name__ for e in variant.effects()] == ["Intergenic"]


HBV = Path(__file__).parent / "data" / "HBV_X04615.1.fa"
HBV_GENES = {
    # GenBank X04615.1: polymerase joined across the origin, and S.
    "P": ((2307, 3215), (1, 1623)),
    "S": ((155, 835),),
}


@pytest.fixture(scope="module")
def hbv(tmp_path_factory):
    sequence = "".join(HBV.read_text().splitlines()[1:])
    rows, transcripts, proteins = [], {}, {}
    for gene, exons in HBV_GENES.items():
        coding = "".join(sequence[start - 1:end] for start, end in exons)
        transcripts["T" + gene] = coding
        proteins["P" + gene] = translate(coding)[:-1]
        rows += gtf_rows("+", exons, gene=gene, transcript="T" + gene,
                         protein="P" + gene, utr5=0, coding=len(coding),
                         contig="X04615.1")
    genome = annotation(tmp_path_factory.mktemp("hbv"), rows, transcripts, proteins)
    return SimpleNamespace(sequence=sequence, proteins=proteins, genome=Genome(
        genome, fasta={"X04615.1": sequence}, verify=False))


@pytest.mark.parametrize("annotator", ANNOTATORS)
@pytest.mark.parametrize("gene,positions", [
    ("P", range(2307, 3215, 23)),   # before the origin
    ("P", range(1, 1623, 41)),      # after it
    ("S", range(155, 835, 17)),
])
def test_hbv_polymerase_and_surface_on_the_real_genome(hbv, annotator, gene, positions):
    exons = HBV_GENES[gene]
    reference = hbv.proteins["P" + gene] + "*"
    assert gene != "P" or len(reference) == 844
    wrong = []
    for position in positions:
        ref = hbv.sequence[position - 1]
        alt = "ACGT".replace(ref, "")[position % 3]
        mutated = hbv.sequence[:position - 1] + alt + hbv.sequence[position:]
        protein = translate("".join(mutated[start - 1:end] for start, end in exons))
        effect, = [e for e in Variant("X04615.1", position, ref, alt, hbv.genome).effects(
            annotator=annotator) if e.transcript.id == "T" + gene]
        name = type(effect).__name__
        if protein == reference:
            ok = name == "Silent"
        elif not protein.startswith("M"):
            ok = name in ("StartLoss", "AlternateStartCodon")
        elif not protein.endswith("*"):
            ok = name == "StopLoss"
        else:
            index = next(i for i, (a, b) in enumerate(zip(reference, protein)) if a != b)
            ok = effect.short_description == "p.%s%d%s" % (
                reference[index], index + 1, protein[index])
        if not ok:
            wrong.append((position, ref, alt, name, effect.short_description))
    assert not wrong
