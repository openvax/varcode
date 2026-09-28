"""The openvax-v2 SV targets Varcode's other esvee tests lack.

openvax-v2, the OpenVax libraries' shared Sid test data (iskandr/osteosarc#56),
lists 12 SV targets. Varcode already tests SV0009, SV0324 and SV0384 through
``osteosarc_esvee_somatic.vcf``, and SV0030, SV0178 and SV0186 through esvee
alleles written into ``test_fusion_insertions.py`` and
``test_osteosarc_fusions.py``. ``tests/data/osteosarc_esvee_shared_svs.vcf``
holds the somatic esvee records for the other six from IPISRC044_tumor_T2_ucla,
each with its mate, with only the FORMAT and sample columns removed:

- SV0055, a chr9 deletion (23175/23378), and SV0499, a chr2 deletion
  (4820/4844);
- SV0175, a chr6 duplication (15481/15487), and SV0461, a chr20 duplication
  (43635/43762);
- SV0402, a chr20 inversion (43636/43709);
- Isovar's DLG5 deletion (26313/26315): the junction chr10:77850914 to
  77930461 with a 24-base insertion, which openvax-v2's recipe now records
  (openvax-v1 had DRAGEN's nominal breakends, 7-9 bases inward). Isovar finds
  that junction sequence, as PURPLE reports it, in every DRAGEN assembled
  contig.

LINX, run on the same calls, names a gene at one breakend each of SV0055
(ATP8B5P), SV0175 (GABBR1) and SV0499 (IMMT) and none at the others. Varcode
must find each named gene at its breakend and no protein-coding transcript
where LINX names none.
"""

from pyensembl import cached_release

from varcode import load_vcf
from varcode.effects import Failure, StructuralVariantEffect
from varcode.transforms import pair_breakends

from .data import data_path

# Release 95, as for the other esvee tests. CI installs it.
ensembl_grch38 = cached_release(95)

# record ID -> (contig, position, mate contig, mate position, esvee SVTYPE)
RECORDS = {
    "23175": ("9", 22_510_295, "9", 35_428_809, "DEL"),
    "23378": ("9", 35_428_809, "9", 22_510_295, "DEL"),
    "4820": ("2", 84_414_540, "2", 86_146_851, "DEL"),
    "4844": ("2", 86_146_851, "2", 84_414_540, "DEL"),
    "26313": ("10", 77_850_914, "10", 77_930_461, "DEL"),
    "26315": ("10", 77_930_461, "10", 77_850_914, "DEL"),
    "15481": ("6", 29_141_010, "6", 29_630_820, "DUP"),
    "15487": ("6", 29_630_820, "6", 29_141_010, "DUP"),
    "43635": ("20", 52_859_903, "20", 58_208_346, "DUP"),
    "43762": ("20", 58_208_346, "20", 52_859_903, "DUP"),
    "43636": ("20", 52_863_939, "20", 55_876_748, "INV"),
    "43709": ("20", 55_876_748, "20", 52_863_939, "INV"),
}
# Genes LINX names at a breakend of an openvax-v2 catalogue target; it names
# none at the other catalogue breakends. (DLG5's records have no LINX row.)
LINX_BREAKEND_GENES = {"23378": {"ATP8B5P"}, "15487": {"GABBR1"}, "4844": {"IMMT"}}
DLG5_RECORDS = {"26313", "26315"}
DLG5_201_ID = "ENST00000372391"


def _records():
    return load_vcf(
        data_path("osteosarc_esvee_shared_svs.vcf"),
        genome=ensembl_grch38,
        parse_structural_variants=True)


def _genes(effects, protein_coding=False):
    return {effect.transcript.gene_name for effect in effects
            if getattr(effect, "transcript", None) is not None
            and (not protein_coding or effect.transcript.biotype == "protein_coding")}


def test_records_load_as_mated_breakends():
    vc = _records()
    assert {v.sv_type for v in vc} == {"BND"}
    assert {
        vc.metadata[v]["id"]: (
            v.contig, v.start, v.mate_contig, v.mate_start, v.info.get("svtype"))
        for v in vc} == RECORDS


def test_mates_pair_into_the_callers_events():
    assert {(event.sv_type, event.contig, event.start, event.end)
            for event in pair_breakends(_records())} == {
        ("DEL", "9", 22_510_295, 35_428_808),
        ("DEL", "2", 84_414_540, 86_146_850),
        ("DEL", "10", 77_850_914, 77_930_460),
        ("DUP", "6", 29_141_009, 29_630_820),
        ("DUP", "20", 52_859_902, 58_208_346),
        ("INV", "20", 52_863_939, 55_876_748),
    }


def test_breakend_genes_agree_with_linx():
    vc = _records()
    for sv in vc:
        record_id = vc.metadata[sv]["id"]
        if record_id in DLG5_RECORDS:
            continue
        effects = sv.effects()
        named = LINX_BREAKEND_GENES.get(record_id, set())
        assert named <= _genes(effects), record_id
        assert _genes(effects, protein_coding=True) <= named, record_id


def test_dlg5_deletion_removes_the_dlg5_start_and_leaves_its_product_unresolved():
    """The deletion takes DLG5-201's start codon and 5' end but keeps its
    3' exons, so no mutant product follows from the DNA alone; Isovar
    reaches the same conclusion from the RNA."""
    (deletion,) = [event for event in pair_breakends(_records()) if event.contig == "10"]
    dlg5_201 = ensembl_grch38.transcript_by_id(DLG5_201_ID)
    assert dlg5_201.on_backward_strand
    assert all(deletion.start < p <= deletion.end
               for p in dlg5_201.start_codon_positions)
    assert deletion.start < dlg5_201.end <= deletion.end
    assert dlg5_201.start <= deletion.start
    (effect,) = [effect for effect in deletion.effects()
                 if getattr(effect, "transcript", None) is not None
                 and effect.transcript.id == DLG5_201_ID]
    assert isinstance(effect, StructuralVariantEffect)
    assert effect.short_description == "unresolved-structural-transcript"
    records = _records()
    by_id = {records.metadata[sv]["id"]: sv for sv in records}
    assert _genes(by_id["26313"].effects()) == {"DLG5"}
    assert _genes(by_id["26315"].effects()) == set()


def test_every_record_and_event_annotates_without_failure():
    vc = _records()
    for sv in list(vc) + list(pair_breakends(vc)):
        effects = sv.effects()
        assert len(effects) > 0
        assert not any(isinstance(effect, Failure) for effect in effects)
