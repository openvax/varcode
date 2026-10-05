"""RNA phasing from complete templates in the published openvax-v2 bundle.

These are observed allele combinations, not independently validated tumor
haplotypes. See tests/README.md for the sources, manual fragment audit and
filtering policy. No alignments or bases are synthesized for these tests.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path

import pytest
from pyensembl import cached_release

from varcode import MolecularPhaseResolver, RNAReadPhasingSource, Variant


SHORT_SOURCE = "2024.06.11.bostongene.align.tcga.protocol.dr32Aligned.sorted"
LONG_SOURCE = "IPISRC044_T1_sclrs_ONT.tagged"
SHORT_NTF3 = SHORT_SOURCE + ".NTF3-chr12-5494381-compound"
LONG_NTF3 = LONG_SOURCE + ".NTF3-chr12-5494381-compound"
LONG_EXOC4 = LONG_SOURCE + ".EXOC4-chr7-133274996"
MEMBERS = (SHORT_NTF3, LONG_NTF3, LONG_EXOC4)
MANIFEST_SHA256 = "b42dbce529cb5724153f035fca3783a41e673b31610a30851b2ef96d4ba43734"

NTF3 = ("12", 5494381, "A", "G")
NTF3_COMPOUND_PARTNER = ("12", 5494382, "G", "T")
# Read-observed companion alleles; no germline/somatic status is assumed.
NTF3_COMPANION = ("12", 5494466, "G", "A")
EXOC4 = ("7", 133274996, "G", "T")
EXOC4_COMPANION = ("7", 133895694, "G", "A")


@pytest.fixture(scope="module")
def shared_phasing_reads(tmp_path_factory):
    osteosarc = pytest.importorskip("osteosarc", reason="install .[test-data]")
    pytest.importorskip("pysam", reason="install .[rna]")
    root = Path(osteosarc.fetch_bundle("openvax-v2"))
    manifest_bytes = (root / "manifest.json").read_bytes()
    assert hashlib.sha256(manifest_bytes).hexdigest() == MANIFEST_SHA256
    recipe = json.loads((root / "recipe.json").read_text())
    assert recipe["selection"]["supplementary"] is True
    compound = recipe["targets"]["NTF3-chr12-5494381-compound"]
    assert (compound["assembly"], compound["contig"], compound["position"],
            compound["ref"], compound["alt"]) == ("GRCh38", "chr12", 5494381, "AG", "GT")
    paths = osteosarc.export_bundle(
        "openvax-v2", tmp_path_factory.mktemp("sid-rna-phasing"), members=MEMBERS)
    return paths, root, json.loads(manifest_bytes)


@pytest.mark.parametrize("member,records,templates,supplementary", [
    (SHORT_NTF3, 12, 6, 0),
    (LONG_NTF3, 15, 15, 0),
    (LONG_EXOC4, 40, 38, 2),
])
def test_shared_phasing_inputs_keep_complete_templates(
        shared_phasing_reads, member, records, templates, supplementary):
    from osteosarc.records import read_records

    paths, root, manifest = shared_phasing_reads
    exported = list(read_records(paths[member]))
    identities = Counter(record.digest for record in exported)
    assert identities == Counter(manifest["members"][member]["records"])
    selected = {record.template for record in exported}
    source = manifest["members"][member]["source"]
    # Compare to every record of each selected RG/QNAME in the bundle's source,
    # rather than refetching a locus and losing mates or distant split pieces.
    complete = Counter(record.digest for record in read_records(
        root / manifest["sources"][source]["bam"]) if record.template in selected)
    assert identities == complete
    assert len(exported) == records
    assert len(selected) == templates
    assert sum(record.read.is_supplementary for record in exported) == supplementary
    if member == SHORT_NTF3:
        for template in selected:
            assert {record.read.flag & 0xc0 for record in exported
                    if record.template == template} == {0x40, 0x80}


@pytest.mark.parametrize("member,alleles,counts,support,phase", [
    (SHORT_NTF3, (NTF3, NTF3_COMPOUND_PARTNER), (3, 0, 0, 3), (3, 3), True),
    (LONG_NTF3, (NTF3, NTF3_COMPOUND_PARTNER), (9, 0, 0, 5), (9, 9), True),
    (LONG_NTF3, (NTF3, NTF3_COMPANION), (0, 9, 3, 2), (9, 3), False),
    (LONG_EXOC4, (EXOC4, EXOC4_COMPANION), (1, 0, 1, 0), (19, 3), None),
])
def test_real_sid_rna_fragment_evidence(
        shared_phasing_reads, member, alleles, counts, support, phase):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81)) for allele in alleles)
    source = RNAReadPhasingSource(
        str(paths[member]), variants=variants,
        # Q20/MAPQ20, five-base edge exclusion and SAM flag filters remain
        # at their defaults. Allow a 5% allele error rate for the ONT data.
        phasing_error_rate=0.01 if member == SHORT_NTF3 else 0.05)
    try:
        assert source._phase_counts(*variants) == counts
        assert tuple(source.supports_variant(v) for v in variants) == support
        assert all(source.has_evidence(v) for v in variants)
        assert source.in_cis(*variants) is phase
        assert source.in_cis(*reversed(variants)) is phase
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(*variants) is phase
        for variant, partner in (variants, tuple(reversed(variants))):
            assert tuple(resolver.phased_partners(variant)) == ((partner,) if phase is True else ())
    finally:
        source.close()


def test_real_sid_cis_call_requires_enough_fragments(shared_phasing_reads):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81))
                     for allele in (NTF3, NTF3_COMPOUND_PARTNER))
    source = RNAReadPhasingSource(
        str(paths[LONG_NTF3]), variants=variants, min_alt_reads=10, phasing_error_rate=0.05)
    try:
        assert tuple(source.supports_variant(v) for v in variants) == (9, 9)
        assert not any(source.has_evidence(v) for v in variants)
        assert source.in_cis(*variants) is None
        assert MolecularPhaseResolver(source).in_cis(*variants) is None
    finally:
        source.close()


@pytest.mark.parametrize("skip_supplementary", [True, False])
def test_real_sid_supplementary_records_obey_the_evidence_filter(
        shared_phasing_reads, skip_supplementary):
    import pysam

    paths, _, _ = shared_phasing_reads
    variant = Variant(*EXOC4, genome=cached_release(81))
    with pysam.AlignmentFile(paths[LONG_EXOC4], "rb") as bam:
        supplementary = [read for read in bam if read.is_supplementary]
    assert len(supplementary) == 2
    source = RNAReadPhasingSource(
        str(paths[LONG_EXOC4]), skip_supplementary=skip_supplementary, phasing_error_rate=0.05)
    try:
        assert [source._read_allele(read, variant) for read in supplementary] == (
            [None, None] if skip_supplementary else ["ref", "ref"])
        assert source.supports_variant(variant) == 19
    finally:
        source.close()
