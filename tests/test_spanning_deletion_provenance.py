# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for VCF ``*`` spanning-deletion placeholder provenance (#264).

Covers:

* Placeholder-only rows retain row/sample/ALT-index provenance without creating
  an annotatable :class:`~varcode.Variant` or protein consequence.
* Mixed-ALT rows retain both the annotatable variant (with original ALT index
  and ``vcf_alt_alleles``) and the ``*`` :class:`~varcode.SpanningDeletionRecord`.
* Upstream deletion resolution links only when source evidence uniquely
  resolves it (sites-only, sample genotype match, phased haplotype match).
* Filtered, absent, genotype-mismatched, trans-phased, and ambiguous upstream
  deletions remain explicitly ``"unresolved"`` or ``"ambiguous"`` with
  ``upstream_deletion is None``.
* ``max_variants`` does not count ``*`` placeholders.
* Serialization (``to_json``/``from_json``, ``pickle``), CSV export header
  accounting, ``spanning_deletions_dataframe()``, and set/transform preservation.
"""

from __future__ import annotations

import pickle

import pytest
from pyensembl import cached_release

from varcode import SpanningDeletionRecord, VariantCollection, load_vcf
from varcode.csv_helpers import read_metadata_header
from varcode.transforms import left_align_indels


ENSEMBL81 = cached_release(81)


def _write_vcf(tmp_path, name, body):
    path = tmp_path / name
    path.write_text(body)
    return str(path)


def test_placeholder_only_row_resolves_unique_upstream_deletion_and_no_effect(tmp_path):
    # VCF deletion at 117559590 REF=ATCTT -> ALT=A deletes 117559591..117559594.
    # Downstream placeholder-only row at 117559592 has ALT=* and GT=0/1.
    path = _write_vcf(
        tmp_path,
        "placeholder_only.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "7\t117559590\tdel1\tATCTT\tA\t90\tPASS\tDP=50\tGT\t0/1\n"
        "7\t117559592\tstar1\tC\t*\t80\tPASS\tDP=48\tGT\t0/1\n",
    )
    vc = load_vcf(path, genome=ENSEMBL81)

    # Only the actual deletion is an annotatable Variant; * creates no effect.
    assert len(vc) == 1
    assert vc[0].is_deletion
    effects = vc.effects()
    assert all(e.variant == vc[0] for e in effects)

    assert len(vc.spanning_deletion_records) == 1
    rec = vc.spanning_deletion_records[0]
    assert isinstance(rec, SpanningDeletionRecord)
    assert rec.contig == "7"
    assert rec.original_contig == "7"
    assert rec.start == 117559592
    assert rec.ref == "C"
    assert rec.alt_allele_index == 0
    assert rec.gt_allele_index == 1
    assert rec.vcf_alt_alleles == ("*",)
    assert rec.id == "star1"
    assert rec.qual == 80.0
    assert rec.filter == ()
    assert rec.info == {"DP": 48}
    assert rec.sample_info == {"tumor": {"GT": "0/1"}}
    assert rec.sample_gt_carriers() == ("tumor",)
    assert rec.upstream_deletion == vc[0]
    assert rec.upstream_deletion_status == "resolved"


def test_mixed_alt_row_preserves_genotype_indices_and_resolves_phased_deletion(tmp_path):
    # Upstream deletion on hap 0 (1|0, PS=100); mixed-ALT row at 117559592 has
    # ALT=G,* with GT=2|1 (so '*' is allele index 2 on hap 0, 'G' is allele 1 on hap 1).
    path = _write_vcf(
        tmp_path,
        "mixed_alt.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "7\t117559590\tdel1\tATCTT\tA\t90\tPASS\t.\tGT:PS\t1|0:100\n"
        "7\t117559592\tmixed1\tC\tG,*\t75\tPASS\t.\tGT:PS\t2|1:100\n",
    )
    vc = load_vcf(path, genome=ENSEMBL81)

    assert len(vc) == 2
    deletion = next(v for v in vc if v.is_deletion)
    snv = next(v for v in vc if v.is_snv)
    assert vc.metadata[snv]["alt_allele_index"] == 0
    assert vc.metadata[snv]["vcf_alt_alleles"] == ("G", "*")

    assert len(vc.spanning_deletion_records) == 1
    rec = vc.spanning_deletion_records[0]
    assert rec.alt_allele_index == 1
    assert rec.gt_allele_index == 2
    assert rec.vcf_alt_alleles == ("G", "*")
    assert rec.sample_gt_carriers() == ("tumor",)
    assert rec.upstream_deletion == deletion
    assert rec.upstream_deletion_status == "resolved"


def test_filtered_absent_trans_and_ambiguous_upstream_deletions(tmp_path):
    # 1. Filtered upstream deletion (LowQual dropped under only_passing=True) -> unresolved
    filtered_path = _write_vcf(
        tmp_path,
        "filtered_del.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "7\t117559590\tdel_lowq\tATCTT\tA\t10\tLowQual\t.\tGT\t0/1\n"
        "7\t117559592\tstar1\tC\t*\t80\tPASS\t.\tGT\t0/1\n",
    )
    vc_filt = load_vcf(filtered_path, genome=ENSEMBL81, only_passing=True)
    assert len(vc_filt) == 0
    assert len(vc_filt.spanning_deletion_records) == 1
    assert vc_filt.spanning_deletion_records[0].upstream_deletion is None
    assert vc_filt.spanning_deletion_records[0].upstream_deletion_status == "unresolved"

    # 2. Phased in TRANS with the upstream deletion -> unresolved (does not mislink!)
    trans_path = _write_vcf(
        tmp_path,
        "trans_del.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "7\t117559590\tdel1\tATCTT\tA\t90\tPASS\t.\tGT:PS\t1|0:100\n"
        "7\t117559592\tstar1\tC\t*\t80\tPASS\t.\tGT:PS\t0|1:100\n",
    )
    vc_trans = load_vcf(trans_path, genome=ENSEMBL81)
    assert vc_trans.spanning_deletion_records[0].upstream_deletion is None
    assert vc_trans.spanning_deletion_records[0].upstream_deletion_status == "unresolved"

    # 3. Sample genotype mismatch (sample carries * but is 0/0 at upstream deletion) -> unresolved
    mismatch_path = _write_vcf(
        tmp_path,
        "mismatch_del.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "7\t117559590\tdel1\tATCTT\tA\t90\tPASS\t.\tGT\t0/1\t0/0\n"
        "7\t117559592\tstar1\tC\t*\t80\tPASS\t.\tGT\t0/0\t0/1\n",
    )
    vc_mm = load_vcf(mismatch_path, genome=ENSEMBL81)
    assert vc_mm.spanning_deletion_records[0].upstream_deletion is None
    assert vc_mm.spanning_deletion_records[0].upstream_deletion_status == "unresolved"

    # 4. Two distinct overlapping upstream deletions in the same sample -> ambiguous
    ambig_path = _write_vcf(
        tmp_path,
        "ambig_del.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "7\t117559589\tdel_a\tGATCTT\tG\t90\tPASS\t.\tGT\t0/1\n"
        "7\t117559590\tdel_b\tATCTT\tA\t90\tPASS\t.\tGT\t0/1\n"
        "7\t117559592\tstar1\tC\t*\t80\tPASS\t.\tGT\t0/1\n",
    )
    vc_amb = load_vcf(ambig_path, genome=ENSEMBL81)
    assert vc_amb.spanning_deletion_records[0].upstream_deletion is None
    assert vc_amb.spanning_deletion_records[0].upstream_deletion_status == "ambiguous"


def test_symbolic_del_resolves_spanning_star_when_sv_parsing_enabled(tmp_path):
    path = _write_vcf(
        tmp_path,
        "sv_del_star.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\tsvdel\tA\t<DEL>\t90\tPASS\tEND=200\n"
        "1\t150\tstar1\tC\t*\t80\tPASS\t.\n",
    )
    vc_sv = load_vcf(path, genome=ENSEMBL81, parse_structural_variants=True)
    assert len(vc_sv) == 1
    assert vc_sv.spanning_deletion_records[0].upstream_deletion == vc_sv[0]
    assert vc_sv.spanning_deletion_records[0].upstream_deletion_status == "resolved"

    # Without parse_structural_variants=True, the SV is skipped so * is unresolved.
    with pytest.warns(UserWarning, match="symbolic/breakend"):
        vc_no_sv = load_vcf(path, genome=ENSEMBL81, parse_structural_variants=False)
    assert len(vc_no_sv) == 0
    assert vc_no_sv.spanning_deletion_records[0].upstream_deletion is None
    assert vc_no_sv.spanning_deletion_records[0].upstream_deletion_status == "unresolved"


def test_max_variants_does_not_count_spanning_star(tmp_path):
    path = _write_vcf(
        tmp_path,
        "limit_star.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t117559590\tdel1\tATCTT\tA\t90\tPASS\t.\n"
        "7\t117559591\tstar1\tT\t*\t80\tPASS\t.\n"
        "7\t117559592\tstar2\tC\t*\t80\tPASS\t.\n"
        "7\t117559595\tsnv1\tT\tG\t80\tPASS\t.\n",
    )
    vc = load_vcf(path, genome=ENSEMBL81, max_variants=2)
    assert len(vc) == 2
    assert len(vc.spanning_deletion_records) == 2
    assert all(r.upstream_deletion_status == "resolved" for r in vc.spanning_deletion_records)


def test_serialization_csv_accounting_dataframe_and_transform_preservation(tmp_path):
    path = _write_vcf(
        tmp_path,
        "roundtrip_star.vcf",
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "chr7\t117559590\tdel1\tATCTT\tA\t90\tPASS\t.\tGT\t0/1\n"
        "chr7\t117559592\tstar1\tC\t*\t80\tPASS\t.\tGT\t0/1\n",
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    assert len(vc.spanning_deletion_records) == 1
    assert vc.spanning_deletion_records[0].contig == "7"
    assert vc.spanning_deletion_records[0].original_contig == "chr7"

    # JSON round-trip
    restored_json = VariantCollection.from_json(vc.to_json())
    assert restored_json == vc
    assert restored_json.exactly_equal(vc)
    assert restored_json.spanning_deletion_records == vc.spanning_deletion_records
    assert restored_json.spanning_deletion_records[0].upstream_deletion == vc[0]

    # Pickle round-trip
    restored_pkl = pickle.loads(pickle.dumps(vc))
    assert restored_pkl.exactly_equal(vc)
    assert restored_pkl.spanning_deletion_records == vc.spanning_deletion_records

    # CSV export accounting: * is not an annotatable row, header records count,
    # and spanning_deletions_dataframe() exposes structured accounting.
    csv_path = tmp_path / "variants.csv"
    vc.to_csv(str(csv_path))
    header = read_metadata_header(str(csv_path))
    assert header["spanning_deletion_records"] == "1"
    reloaded_csv = VariantCollection.from_csv(str(csv_path))
    assert len(reloaded_csv) == 1

    sdf = vc.spanning_deletions_dataframe()
    assert len(sdf) == 1
    assert sdf.iloc[0]["contig"] == "7"
    assert sdf.iloc[0]["original_contig"] == "chr7"
    assert sdf.iloc[0]["gt_allele_index"] == 1
    assert sdf.iloc[0]["upstream_deletion_status"] == "resolved"
    assert sdf.iloc[0]["samples_carrying_star"] == "tumor"

    # Set ops and transforms preserve spanning_deletion_records
    assert vc.union(vc).spanning_deletion_records == vc.spanning_deletion_records
    assert vc.intersection(vc).spanning_deletion_records == vc.spanning_deletion_records
    assert left_align_indels(vc).spanning_deletion_records == vc.spanning_deletion_records
