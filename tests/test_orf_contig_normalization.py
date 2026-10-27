"""Source contig spelling must survive normalized variant/RNA inputs (#636)."""

from dataclasses import replace

import pytest
from pyensembl import cached_release

from varcode import Variant, annotate_orf_variants, annotate_orf_rna_product
from varcode.effect_hypotheses import RealizedTranscriptProduct
from varcode.nucleotides import reverse_complement
from .test_uorf import Transcript


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("source,normalized", [("Mito", "MITO"), ("mtdna", "MTDNA"), ("chrM", "MT")])
def test_source_names_match_normalized_inputs_without_losing_declared_code(source, normalized, reverse):
    transcript = Transcript(reverse=reverse)
    transcript.contig = normalized
    transcript.genome = cached_release(115, species="yeast")
    # Synthetic ORF under the yeast code; its mutation has different meaning
    # under table 3 (CTG -> Thr) than under standard/table 2 (CTG -> Leu).
    record = replace(transcript.record(), contig=source, genome_build="R64-1-1",
                     taxon_id=4932, species="saccharomyces_cerevisiae", translation_table=3)
    reference = record.bind(transcript)
    ref, alt = "CCC", "CTG"
    if reverse:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    variant = Variant(source, min(transcript.positions[9:12]), ref, alt, transcript.genome,
                      convert_ucsc_contig_names=True)
    assert variant.contig == normalized
    result = annotate_orf_variants([variant], reference)
    assert result.reference.model.contig == source
    assert result.reference.model.translation_table == 3
    assert result.outcomes[0].mutant.protein_sequence == "MKTG"
    assert result.outcomes[0].mutant.reference_evidence["mutant_translation_evidence"] == "not_assessed"
    product = RealizedTranscriptProduct(
        transcript, transcript.sequence[:9] + "CTG" + transcript.sequence[12:], None,
        evidence={"base_origins": [(normalized, p, transcript.strand) for p in transcript.positions]})
    observed = annotate_orf_rna_product(reference, product)
    assert observed.product.protein_sequence == "MKTG"
    wrong = Variant("2", variant.start, variant.ref, variant.alt, transcript.genome)
    with pytest.raises(ValueError, match="contig mismatch"):
        annotate_orf_variants([wrong], reference)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("release,source,normalized", [
    (75, "chrUn_gl000216", "GL000216.1"),
    (81, "chrUn_KI270438v1", "KI270438.1"),
])
def test_ucsc_accession_aliases_match_orf_variants_and_rna_origins(release, source, normalized, reverse):
    transcript = Transcript(reverse=reverse)
    transcript.contig = source
    transcript.genome = cached_release(release)
    record = replace(transcript.record(), contig=source,
                     genome_build=transcript.genome.reference_name)
    reference = record.bind(transcript)
    ref, alt = "CCC", "CTG"
    if reverse:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    variant = Variant(source, min(transcript.positions[9:12]), ref, alt, transcript.genome,
                      convert_ucsc_contig_names=True)
    assert variant.contig == normalized
    result = annotate_orf_variants([variant], reference)
    assert result.reference.model.contig == source
    assert result.outcomes[0].mutant.protein_sequence == "MKLG"
    product = RealizedTranscriptProduct(
        transcript, transcript.sequence[:9] + "CTG" + transcript.sequence[12:], None,
        evidence={"base_origins": [(normalized, p, transcript.strand) for p in transcript.positions]})
    observed = annotate_orf_rna_product(reference, product)
    assert observed.product.protein_sequence == "MKLG"
    wrong_version = Variant(normalized.rsplit(".", 1)[0] + ".2", variant.start,
                            variant.ref, variant.alt, transcript.genome)
    with pytest.raises(ValueError, match="contig mismatch"):
        annotate_orf_variants([wrong_version], reference)
