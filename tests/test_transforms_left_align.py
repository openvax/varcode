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

"""Tests for :func:`varcode.transforms.left_align_indels` (#369).

Covers the test matrix from PR #374:

* Pass-throughs — SNVs, MNVs, complex variants, already-canonical
  indels.
* Real shifts — homopolymer deletions/insertions, STR repeats,
  multi-base payloads.
* Bounded shifts — chromosome start, exon boundary (partial-shift
  flag), uncovered intronic positions.
* Metadata + provenance — source_variants population, original_start
  in info dict, source-path metadata re-keyed under shifted variants.
* Composition — idempotence, mixed VC with SNVs + indels + paired
  BNDs, empty VC.
"""

import pytest
from pyensembl import cached_release

from varcode import Genome, Variant, VariantCollection
from varcode.transforms import left_align_indels


# --------------------------------------------------------------------
# Test helpers — DictBackedFasta adapter (same shape as
# tests/test_genome_sequence.py).
# --------------------------------------------------------------------


class _DictBackedFasta:
    """``{contig: sequence_string}`` adapter shaped like pyfaidx.Fasta.

    ``offset`` is the 1-based genome position of the first character
    of each contig's sequence."""
    def __init__(self, contig_to_seq, *, offset=1):
        self._d = dict(contig_to_seq)
        self._offset = offset

    def __getitem__(self, contig):
        return _Slicer(self._d.get(contig, ""), self._offset)


class _Slicer:
    def __init__(self, seq, offset):
        self.seq = seq
        self._offset = offset

    def __getitem__(self, s):
        lo = max(0, (s.start or 0) - self._offset + 1)
        hi = max(0, (s.stop or 0) - self._offset + 1)
        return _Span(self.seq[lo:hi])


class _Span:
    def __init__(self, seq):
        self.seq = seq


def _genome_with_fasta(contig_to_seq, offset=1):
    """Build a :class:`varcode.Genome` wrapping pyensembl release 81
    with an attached :class:`_DictBackedFasta`. verify=False because
    the test fixtures don't agree with pyensembl transcript bases."""
    return Genome(
        cached_release(81),
        fasta=_DictBackedFasta(contig_to_seq, offset=offset),
        verify=False)


# --------------------------------------------------------------------
# Pass-throughs
# --------------------------------------------------------------------


def test_snv_passes_through():
    g = _genome_with_fasta({"1": "A" * 100}, offset=1)
    v = Variant("1", 50, "A", "T", genome=g)
    vc = VariantCollection([v])
    out = left_align_indels(vc)
    assert out[0] is v
    assert out[0].source_variants == ()


def test_mnv_passes_through():
    g = _genome_with_fasta({"1": "A" * 100}, offset=1)
    v = Variant("1", 50, "AT", "GC", genome=g)
    vc = VariantCollection([v])
    out = left_align_indels(vc)
    assert out[0] is v
    assert out[0].source_variants == ()


def test_complex_variant_passes_through():
    """Length-mismatched AND content-different (``ATG->GCC``) — not a
    clean indel. Pass through unchanged."""
    g = _genome_with_fasta({"1": "A" * 100}, offset=1)
    v = Variant("1", 50, "ATG", "GCC", genome=g)
    vc = VariantCollection([v])
    out = left_align_indels(vc)
    assert out[0] is v


def test_canonical_deletion_passes_through():
    """Already-canonical deletion in a non-repeat region — no shift."""
    # Reference at pos 49 = 'C', pos 50 = 'T', pos 51 = 'G'. Deleting T
    # at pos 50 has no equivalent position because pos 49 != T.
    # ``Variant("1", 50, "CT", "C", g)`` normalizes to ref='T',
    # alt='', start=51 (varcode's deletion-start points at the
    # deleted T). Reference layout: C at pos 50 + T at pos 51 +
    # G at pos 52..onward. The algorithm checks reference[50]='C'
    # vs payload[-1]='T'; mismatch -> no shift.
    g = _genome_with_fasta({"1": "C" * 50 + "T" + "G" * 49}, offset=1)
    v = Variant("1", 50, "CT", "C", genome=g)
    assert (v.ref, v.alt, v.start) == ("T", "", 51)
    vc = VariantCollection([v])
    out = left_align_indels(vc)
    assert out[0] is v
    assert out[0].source_variants == ()


def test_canonical_insertion_passes_through():
    """Already-canonical insertion in a non-repeat region — no shift."""
    g = _genome_with_fasta({"1": "C" * 49 + "T" + "G" * 50}, offset=1)
    v = Variant("1", 50, "T", "TA", genome=g)  # insert A after pos 50
    vc = VariantCollection([v])
    out = left_align_indels(vc)
    assert out[0] is v


# --------------------------------------------------------------------
# Real shifts
# --------------------------------------------------------------------


def test_deletion_in_homopolymer_shifts_to_leftmost():
    """Reference AAAAAAAAAA at pos 1-10. Delete 'A' at pos 10 (varcode
    normalizes 'AA'->'A' at pos 9 to ref='A',alt='',start=10).
    Should shift to the leftmost equivalent: start=1."""
    g = _genome_with_fasta({"1": "A" * 10}, offset=1)
    v = Variant("1", 9, "AA", "A", genome=g)
    assert v.start == 10  # post-normalization

    out = left_align_indels(VariantCollection([v]))
    o = out[0]
    assert o.ref == "A" and o.alt == ""
    assert o.start == 1
    assert o.source_variants == (v,)


def test_insertion_in_homopolymer_shifts_to_leftmost():
    """Reference AAAAA at pos 5-9, with non-A flanking. Insert 'A'
    after pos 9 should shift to insertion at pos 5 (the leftmost
    anchor where the inserted A is equivalent)."""
    g = _genome_with_fasta({"1": "T" + "A" * 5 + "T"}, offset=4)
    # offset=4 means: pos 4='T', pos 5..9='A', pos 10='T'.
    v = Variant("1", 9, "A", "AA", genome=g)  # insert A after pos 9
    out = left_align_indels(VariantCollection([v]))
    o = out[0]
    assert o.ref == "" and o.alt == "A"
    assert o.start == 4  # algorithm stops when reference[5]='A' but
    # check at iter 6 is at reference[4]='T'; T != A; stop. Wait, let me
    # trace: start=9, check reference[9]='A' (since pos 9 = seq[5] = 'A'),
    # match, shift to start=8. check reference[8]='A', shift. ... shift
    # to start=5. check reference[5]='A', shift to start=4. check
    # reference[4]='T', T != A, stop. So final start=4. ✓


def test_str_repeat_deletion_shifts_to_first_unit():
    """Reference: ATATATATAT (5 AT repeats) at pos 1-10. Delete 'AT'
    at pos 9-10 should shift to delete 'AT' at pos 1-2."""
    g = _genome_with_fasta({"1": "AT" * 5}, offset=1)
    v = Variant("1", 8, "TAT", "T", genome=g)
    # CT- repeat deletion. Varcode normalizes 'TAT'->'T' (anchor T)
    # so ref='AT' alt='' start=9.
    assert (v.ref, v.alt, v.start) == ("AT", "", 9)

    out = left_align_indels(VariantCollection([v]))
    o = out[0]
    assert o.ref == "AT" and o.alt == ""
    assert o.start == 1  # leftmost AT in the repeat
    assert o.source_variants == (v,)


def test_bounded_at_chromosome_start():
    """Shift bounded by start > 1 — even in an infinite homopolymer,
    the algorithm stops at pos 1."""
    # Reference A at pos 1 only, then arbitrary.
    g = _genome_with_fasta({"1": "A" + "T" * 100}, offset=1)
    v = Variant("1", 1, "AA", "A", genome=g)
    # Hmm — this normalizes to ref='A',alt='',start=2 (delete the second
    # base which is T per our reference). That's not a valid call for
    # left-alignment because there's no run.
    # Use a 2-A homopolymer at start instead.
    g2 = _genome_with_fasta({"1": "AA" + "T" * 100}, offset=1)
    v2 = Variant("1", 1, "AA", "A", genome=g2)
    # ref='A',alt='',start=2. Try to shift to 1.
    out = left_align_indels(VariantCollection([v2]))
    o = out[0]
    # Shifted from start=2 to start=1; can't go further (start > 1 false).
    assert o.start == 1
    assert o.source_variants == (v2,)


def test_multi_base_insertion_in_homopolymer():
    """Multi-base payload shift: insert 'AAA' in a long A-run."""
    g = _genome_with_fasta({"1": "T" + "A" * 20 + "T"}, offset=1)
    # pos 1='T', pos 2..21='A', pos 22='T'.
    v = Variant("1", 21, "A", "AAAA", genome=g)
    # Normalized: ref='',alt='AAA',start=21.
    out = left_align_indels(VariantCollection([v]))
    o = out[0]
    assert o.ref == "" and o.alt == "AAA"
    # Shift continues as long as reference matches payload[-1]='A'.
    # Stops when reference at check_pos is not 'A'. Reference[1]='T'.
    # So shift continues through positions 21..2 (all A), stops at
    # check_pos=1 (T != A). Final start = 1 (varcode insertion anchor).
    assert o.start == 1
    assert o.source_variants == (v,)


# --------------------------------------------------------------------
# Idempotence, mixed VC, empty VC
# --------------------------------------------------------------------


def test_idempotent_on_canonical_input():
    """Running left_align_indels twice produces the same VC the
    second time — the first call leaves variants at canonical
    positions, so the second call finds nothing to shift."""
    g = _genome_with_fasta({"1": "A" * 10 + "T" * 5}, offset=1)
    v = Variant("1", 9, "AA", "A", genome=g)
    once = left_align_indels(VariantCollection([v]))
    twice = left_align_indels(once)
    assert once[0].start == twice[0].start
    assert once[0].ref == twice[0].ref
    assert once[0].alt == twice[0].alt
    # The second call had nothing to shift, so it returns the input
    # collection unchanged.
    assert twice is once


def test_mixed_vc_only_indels_shift():
    g = _genome_with_fasta({"1": "A" * 20 + "T" * 80}, offset=1)
    snv = Variant("1", 50, "T", "G", genome=g)
    indel = Variant("1", 19, "AA", "A", genome=g)
    mnv = Variant("1", 60, "TT", "GG", genome=g)
    vc = VariantCollection([snv, indel, mnv])
    out = left_align_indels(vc)
    # Same length; SNV and MNV pass through; indel shifts.
    assert len(out) == 3
    snv_out = [v for v in out if v.is_snv][0]
    indel_out = [v for v in out if v.is_indel][0]
    assert snv_out is snv
    assert indel_out.source_variants == (indel,)
    assert indel_out.start == 1


def test_empty_vc_passes_through():
    g = _genome_with_fasta({"1": "A" * 10}, offset=1)
    # Skip the genome here — empty VC doesn't need one.
    vc = VariantCollection([])
    out = left_align_indels(vc)
    assert len(out) == 0


# --------------------------------------------------------------------
# Metadata + provenance
# --------------------------------------------------------------------


def test_original_start_recorded_in_metadata_info():
    """The original (pre-shift) start position is recorded in the
    shifted variant's metadata info dict for round-trip and
    debugging."""
    g = _genome_with_fasta({"1": "A" * 10}, offset=1)
    v = Variant("1", 9, "AA", "A", genome=g)
    metadata = {"synth": {v: {"id": "v1", "info": {"AC": 1}}}}
    vc = VariantCollection(
        variants=[v],
        sources={"synth"},
        source_to_metadata_dict=metadata)
    out = left_align_indels(vc)
    shifted = out[0]
    assert shifted is not v
    meta = out.source_to_metadata_dict["synth"][shifted]
    assert meta["info"]["original_start"] == 10  # pre-shift start
    assert meta["info"]["AC"] == 1  # other info fields preserved
    # The original variant is no longer keyed in the metadata dict.
    assert v not in out.source_to_metadata_dict["synth"]


def test_unshifted_variant_metadata_passes_through_by_reference():
    """Variants that don't shift keep their original metadata entry
    by reference — no allocation, no info-dict copy."""
    g = _genome_with_fasta({"1": "ACGT" * 25}, offset=1)
    v = Variant("1", 50, "A", "T", genome=g)  # SNV; no shift possible
    original_info = {"AC": 1}
    metadata = {"synth": {v: {"id": "v1", "info": original_info}}}
    vc = VariantCollection(
        variants=[v],
        sources={"synth"},
        source_to_metadata_dict=metadata)
    out = left_align_indels(vc)
    assert out[0] is v
    # Same metadata entry (no shifting happened, so no copy needed).
    assert out.source_to_metadata_dict["synth"][v] is metadata["synth"][v]


# --------------------------------------------------------------------
# Composition with pair_breakends
# --------------------------------------------------------------------


def test_composes_with_pair_breakends():
    """Mixed VC containing SNVs, indels, and a paired BND. After
    left_align_indels + pair_breakends: indel shifts, BNDs collapse
    into one combined SV, SNV unchanged."""
    from varcode.transforms import pair_breakends
    from varcode import load_vcf
    import tempfile

    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"m\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1000\tsnv\tA\tT\t100\tPASS\t.\n"
        "1\t1100\tindel\tAA\tA\t100\tPASS\t.\n"
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\n"
    )
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False) as f:
        f.write(body)
        path = f.name

    # Use a real pyensembl genome (no FASTA) — left_align will fall
    # through to transcript cDNA. The chr1 indel won't shift (Tier 2
    # only covers exonic positions, and pos 1100 is intergenic), but
    # the test still verifies composition works.
    vc = load_vcf(path, genome="GRCh38", parse_structural_variants=True)
    assert len(vc) == 4

    vc = left_align_indels(vc)
    vc = pair_breakends(vc)

    # 1 SNV + 1 indel (unshifted, intergenic with no FASTA) + 1 combined BND.
    assert len(vc) == 3
    bnds = [v for v in vc if getattr(v, "sv_type", None) == "BND"]
    assert len(bnds) == 1
    assert len(bnds[0].source_variants) == 2


# --------------------------------------------------------------------
# Tier 2 coverage — exonic indels shift, intronic ones don't
# --------------------------------------------------------------------


def test_intronic_indel_without_fasta_passes_through():
    """An indel in an intron without a FASTA has no reference
    coverage (Tier 2 returns empty); the algorithm exits immediately
    and the variant passes through unchanged."""
    # Use real pyensembl genome. CFTR (chr7, + strand) — pick an
    # intronic position.
    from pyensembl import cached_release
    g = cached_release(81)
    t = g.transcript_by_id("ENST00000003084")
    intron_pos = t.exons[0].end + 100  # 100bp into the first intron
    v = Variant("7", intron_pos, "AA", "A", genome=g)
    vc = VariantCollection([v])
    out = left_align_indels(vc)
    # No shift — variant passes through.
    assert out[0] is v
    assert out[0].source_variants == ()


def test_partial_shift_flag_set_when_coverage_runs_out():
    """When the algorithm shifts at least once but stops because the
    reference returns empty (Tier 2 exon boundary, or off-end-of-FASTA),
    info['left_align_partial'] = True flags the result."""
    # Reference covers only positions 5-10 (all A); positions <5 are
    # uncovered. A deletion in the run can shift from pos 10 down to
    # pos 5, then stops at the coverage boundary.
    g = _genome_with_fasta({"1": "A" * 6}, offset=5)
    v = Variant("1", 9, "AA", "A", genome=g)
    # Normalized: ref='A', alt='', start=10.
    metadata = {"synth": {v: {"id": "v1", "info": {}}}}
    vc = VariantCollection(
        variants=[v],
        sources={"synth"},
        source_to_metadata_dict=metadata)
    out = left_align_indels(vc)
    shifted = out[0]
    assert shifted.start == 5  # shifted to coverage boundary
    info = out.source_to_metadata_dict["synth"][shifted]["info"]
    assert info.get("left_align_partial") is True
    assert info["original_start"] == 10
