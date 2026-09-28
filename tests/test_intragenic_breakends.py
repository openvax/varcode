"""A breakend junction with both ends in one gene (#551).

For a transcript read across it, the junction fixes the product: the
DEL, DUP or INV its kept sides describe. A caller-labeled breakend pair
is typed as that event, so a lone record must annotate the gene's
transcripts exactly as its labeled pair does, while its effects still
report the record itself.
"""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant
from varcode.effects import GeneFusion, Intronic

from .test_structural_variant_annotator import _pair_direct_breakends

# Release 95, as for the other esvee tests. CI installs it.
genome = cached_release(95)


def _summary(effect):
    primary = getattr(effect, "most_likely_effect", effect)
    model = getattr(effect, "mutant_transcript", None)
    return (type(effect).__name__, type(primary).__name__, effect.short_description,
            getattr(model, "cdna_sequence", None),
            getattr(model, "mutant_protein_sequence", None))


def test_atp5mg_kmt2a_junction_is_an_intron_of_the_readthrough_transcript():
    """openvax-v2's ATP5MG--KMT2A junction is exactly intron 1 of
    Ensembl's ATP5MG-KMT2A readthrough transcript AP001267.5-202."""
    record = StructuralVariant("chr11", 118_401_717, "BND", alt="N[chr11:118468775[",
                               genome=genome, convert_ucsc_contig_names=True)
    deletion = StructuralVariant("11", 118_401_717, "DEL", end=118_468_774,
                                 affected_start=118_401_718, genome=genome)
    readthrough = genome.transcript_by_id("ENST00000648261")
    assert (readthrough.exons[0].end, readthrough.exons[1].start) == (118_401_717, 118_468_775)
    effects = {e.transcript.transcript_id: e for e in record.effects()}
    assert isinstance(effects["ENST00000648261"], Intronic)
    # AP001267.5-201 ends inside the junction's span, so it loses its 3' end.
    other = effects["ENST00000649464"]
    assert _summary(other) == _summary(deletion.effect_on_transcript(other.transcript))
    assert other.mutant_transcript.evidence["sv_type"] == "BND"
    assert other.mutant_transcript.evidence["junction_event"] == "DEL"
    # ATP5MG's gene doesn't reach KMT2A, so its isoforms still fuse.
    assert {type(e) for e in effects.values() if e.transcript.gene_name == "ATP5MG"
            and e.transcript.is_protein_coding} == {GeneFusion}
    assert all(e.variant is record for e in effects.values())


def _records(contig, low, high, kind):
    """The two breakend records a caller writes for a DEL, DUP or INV
    between ``low`` and ``high``, labeled with it."""
    low_alt, high_alt = {
        "DEL": ("N[%s:%d[" % (contig, high), "]%s:%d]N" % (contig, low)),
        "DUP": ("]%s:%d]N" % (contig, high), "N[%s:%d[" % (contig, low)),
        "INV-left": ("N]%s:%d]" % (contig, high), "N]%s:%d]" % (contig, low)),
        "INV-right": ("[%s:%d[N" % (contig, high), "[%s:%d[N" % (contig, low)),
    }[kind]
    label = kind.split("-")[0]
    a = StructuralVariant(contig, low, "BND", alt=low_alt, genome=genome,
                          info={"SVTYPE": label, "MATEID": "b"})
    b = StructuralVariant(contig, high, "BND", alt=high_alt, genome=genome,
                          info={"SVTYPE": label, "MATEID": "a"})
    return a, b


@pytest.mark.parametrize("transcript_id", [
    "ENST00000003084",  # CFTR-201, forward strand
    "ENST00000357654",  # BRCA1-203, reverse strand
])
@pytest.mark.parametrize("kind", ["DEL", "DUP", "INV-left", "INV-right"])
@pytest.mark.parametrize("span", ["exonic", "across_exons"])
def test_lone_record_annotates_like_its_labeled_pair(transcript_id, kind, span):
    transcript = genome.transcript_by_id(transcript_id)
    exons = sorted(transcript.exons, key=lambda exon: exon.start)
    if span == "exonic":
        low = exons[4].start + 10
        high = low + 22
    else:
        low = (exons[3].end + exons[4].start) // 2
        high = (exons[6].end + exons[7].start) // 2
    a, b = _records(transcript.contig, low, high, kind)
    typed = _pair_direct_breakends(a, b)
    assert typed.sv_type == kind.split("-")[0]
    expected = _summary(typed.effect_on_transcript(transcript))
    for record in (a, b):
        effect = record.effect_on_transcript(transcript)
        assert effect.variant is record
        assert _summary(effect) == expected
