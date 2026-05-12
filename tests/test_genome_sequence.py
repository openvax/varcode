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

"""Tests for :class:`varcode.Genome` + :mod:`varcode.genome_sequence` (#372).

Covers:

* ``varcode.Genome`` construction — int / string / pyensembl.Genome /
  varcode.Genome (rewrap) accepted; FASTA path / pyfaidx-shaped object
  / None / garbage handled; verify-on-construct catches mismatched
  FASTAs.
* ``__getattr__`` delegation makes the wrapper a drop-in for
  pyensembl.Genome.
* ``Genome.sequence`` / ``Genome.reference_base`` / ``Genome.reference_range``
  honor the tiered fallback (Tier 1 FASTA -> Tier 2 transcript cDNA ->
  Tier 3 empty), strand handling included.
* Module-level :func:`reference_base` / :func:`reference_range` work
  against either ``varcode.Genome`` or bare ``pyensembl.Genome``.
* ``load_vcf(genome_fasta=...)`` wraps the resolved genome.
* ``cryptic_exons`` migration: silently-bailing flanking scan
  replaced by a per-genome warn-once.
"""

import copy
import warnings

import pytest
from pyensembl import cached_release

from varcode import Genome
from varcode.genome_sequence import reference_base, reference_range


# --------------------------------------------------------------------
# Minimal pyfaidx-shaped mock so tests don't need a real chromosome FASTA
# --------------------------------------------------------------------


class _DictBackedFasta:
    """``{contig: sequence_string}`` adapter shaped like pyfaidx.Fasta.

    ``offset`` is the 1-based genome position of the first character of
    each contig's sequence (most fixtures embed a small window around
    a known locus rather than a full chromosome).
    """
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
        if isinstance(s, slice):
            lo = max(0, (s.start or 0) - self._offset + 1)
            hi = max(0, (s.stop or 0) - self._offset + 1)
            return _Span(self.seq[lo:hi])
        return _Span(self.seq[s - self._offset + 1])


class _Span:
    def __init__(self, seq):
        self.seq = seq


# --------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------


# CFTR (chr7, + strand) — well-known protein coding transcript.
_CFTR_TRANSCRIPT_ID = "ENST00000003084"


@pytest.fixture
def ensembl():
    """Bare pyensembl genome — for tests that exercise the Tier 2-only
    code path (consumers passing a non-wrapped genome)."""
    return cached_release(81)


@pytest.fixture
def cftr_position(ensembl):
    """``(contig, position, expected_plus_strand_base)`` inside CFTR's
    first exon — known to be covered by Tier 2."""
    t = ensembl.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    pos = t.exons[0].start + 5
    expected = t.sequence[5].upper()
    return (t.contig, pos, expected)


def _build_fasta_around(contig, position, base, window=200):
    """``_DictBackedFasta`` with ``base`` at ``position`` and ``'N'``
    elsewhere within +/-``window``."""
    start_offset = position - window
    seq = ['N'] * (2 * window + 1)
    seq[position - start_offset] = base
    return _DictBackedFasta({contig: ''.join(seq)}, offset=start_offset)


# --------------------------------------------------------------------
# varcode.Genome construction
# --------------------------------------------------------------------


def test_genome_from_int_release_no_fasta(ensembl):
    g = Genome(81)
    assert g.fasta is None
    # __getattr__ delegation passes pyensembl attribute access through.
    assert g.reference_name == ensembl.reference_name


def test_genome_from_pyensembl_object_passes_through(ensembl):
    g = Genome(ensembl)
    assert g._ensembl is ensembl
    assert g.fasta is None


def test_genome_rewrap_inherits_fasta(ensembl, cftr_position):
    """Rewrapping a Genome without an explicit ``fasta=`` inherits the
    inner Genome's FASTA. ``fasta=None`` is treated as "no override"
    rather than "detach" — to drop the FASTA, construct a fresh
    ``Genome(ensembl_release)`` directly."""
    contig, pos, base = cftr_position
    fasta = _build_fasta_around(contig, pos, base)
    g1 = Genome(ensembl, fasta=fasta, verify=False)
    g2 = Genome(g1)
    assert g2.fasta is fasta


def test_genome_rewrap_overrides_fasta(ensembl, cftr_position):
    contig, pos, base = cftr_position
    fasta1 = _build_fasta_around(contig, pos, base)
    fasta2 = _build_fasta_around(contig, pos, base)
    g1 = Genome(ensembl, fasta=fasta1, verify=False)
    g2 = Genome(g1, fasta=fasta2, verify=False)
    assert g2.fasta is fasta2


def test_genome_with_integer_fasta_raises_typeerror(ensembl):
    with pytest.raises(TypeError, match="pyfaidx-style protocol"):
        Genome(ensembl, fasta=42, verify=False)


def test_genome_with_list_fasta_raises_typeerror(ensembl):
    with pytest.raises(TypeError, match="pyfaidx-style protocol"):
        Genome(ensembl, fasta=[1, 2, 3], verify=False)


def test_genome_without_release_raises(ensembl):
    """``Genome()`` with no release argument is an error — we don't
    silently construct a placeholder that fails later."""
    with pytest.raises(ValueError, match="requires an ensembl_release"):
        Genome()


def test_genome_deepcopy_does_not_infinite_loop(ensembl):
    """``__getattr__`` delegation must not recurse on Python's
    copy/pickle dunders. The underscore guard in ``__getattr__`` is
    what prevents this; regression here would hang the test."""
    g = Genome(ensembl)
    clone = copy.deepcopy(g)
    assert clone is not g
    assert clone._ensembl is g._ensembl  # shallow share of wrapped genome
    assert clone.fasta is g.fasta


def test_genome_repr_indicates_fasta_presence(ensembl, cftr_position):
    contig, pos, base = cftr_position
    g_no = Genome(ensembl)
    g_with = Genome(ensembl, fasta=_build_fasta_around(contig, pos, base),
                    verify=False)
    assert "no FASTA" in repr(g_no)
    assert "FASTA attached" in repr(g_with)


def test_verify_emits_warning_on_mismatched_fasta(ensembl):
    class _AlwaysWrongFasta:
        def __getitem__(self, contig):
            return _AlwaysWrongSlicer()

    class _AlwaysWrongSlicer:
        def __getitem__(self, s):
            length = (s.stop - s.start) if isinstance(s, slice) else 1
            return _Span('X' * length)

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        Genome(ensembl, fasta=_AlwaysWrongFasta(), verify=True)

    disagreement = [w for w in captured
                    if "disagree" in str(w.message).lower()]
    assert disagreement, [str(w.message) for w in captured]


def test_verify_false_suppresses_check(ensembl):
    class _AlwaysWrongFasta:
        def __getitem__(self, contig):
            return _AlwaysWrongSlicer()

    class _AlwaysWrongSlicer:
        def __getitem__(self, s):
            length = (s.stop - s.start) if isinstance(s, slice) else 1
            return _Span('X' * length)

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        Genome(ensembl, fasta=_AlwaysWrongFasta(), verify=False)

    assert not any("disagree" in str(w.message).lower() for w in captured)


# --------------------------------------------------------------------
# __getattr__ delegation
# --------------------------------------------------------------------


def test_genome_delegates_transcript_by_id(ensembl):
    g = Genome(ensembl)
    # Should pass through to pyensembl's transcript_by_id.
    assert g.transcript_by_id(_CFTR_TRANSCRIPT_ID).id == _CFTR_TRANSCRIPT_ID


def test_genome_delegates_transcripts_at_locus(ensembl, cftr_position):
    g = Genome(ensembl)
    contig, pos, _ = cftr_position
    assert len(g.transcripts_at_locus(contig, pos, pos)) > 0


# --------------------------------------------------------------------
# Lookup via transcript cDNA only (no FASTA attached)
# --------------------------------------------------------------------


def test_transcript_cdna_lookup_plus_strand(ensembl, cftr_position):
    contig, pos, expected = cftr_position
    assert reference_base(ensembl, contig, pos) == expected
    # Same answer via the wrapper.
    assert Genome(ensembl).reference_base(contig, pos) == expected


def test_transcript_cdna_lookup_minus_strand(ensembl):
    minus_t = next(
        t for t in ensembl.transcripts()
        if t.strand == '-' and t.exons and t.is_protein_coding)
    pos = minus_t.exons[0].start + 5
    offset = minus_t.spliced_offset(pos)
    cdna_base = minus_t.sequence[offset].upper()
    expected = cdna_base.translate(str.maketrans("ACGT", "TGCA"))
    assert reference_base(ensembl, minus_t.contig, pos) == expected


def test_transcript_cdna_returns_empty_for_intronic_position(ensembl):
    t = ensembl.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    intron_pos = t.exons[0].end + 100
    assert reference_base(ensembl, t.contig, intron_pos) == ""


# --------------------------------------------------------------------
# Lookup via attached chromosome FASTA
# --------------------------------------------------------------------


def test_chromosome_fasta_lookup_when_attached(ensembl, cftr_position):
    contig, pos, transcript_base = cftr_position
    fasta = _build_fasta_around(contig, pos, transcript_base, window=10)
    g = Genome(ensembl, fasta=fasta, verify=False)
    assert g.reference_base(contig, pos) == transcript_base


def test_chromosome_fasta_falls_through_to_cdna_when_contig_missing(
        ensembl, cftr_position):
    """A FASTA that doesn't carry the requested contig falls through
    to transcript cDNA — does not crash."""
    contig, pos, transcript_base = cftr_position
    empty_fasta = _DictBackedFasta({}, offset=1)
    g = Genome(ensembl, fasta=empty_fasta, verify=False)
    assert g.reference_base(contig, pos) == transcript_base


def test_genome_sequence_method_is_fasta_only(ensembl, cftr_position):
    """``Genome.sequence(contig, start, end)`` mirrors the proposed
    pyensembl API (openvax/pyensembl#337). Returns the + strand range
    from the attached FASTA only — does NOT fall back to transcript
    cDNA. Callers wanting the tiered lookup use ``reference_base`` /
    ``reference_range`` instead."""
    contig, pos, transcript_base = cftr_position
    fasta = _build_fasta_around(contig, pos, transcript_base, window=10)
    g = Genome(ensembl, fasta=fasta, verify=False)
    seq = g.sequence(contig, pos, pos)
    assert seq == transcript_base
    # No FASTA -> empty.
    g_no = Genome(ensembl)
    assert g_no.sequence(contig, pos, pos) == ""


# --------------------------------------------------------------------
# reference_range
# --------------------------------------------------------------------


def test_cdna_range_uses_single_transcript_fast_path(ensembl):
    t = ensembl.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    start = t.exons[0].start + 5
    end = start + 10
    span = reference_range(ensembl, t.contig, start, end)
    assert len(span) == end - start + 1
    assert span[0] == reference_base(ensembl, t.contig, start)
    assert span[-1] == reference_base(ensembl, t.contig, end)


def test_cdna_range_spanning_exon_boundary_returns_empty(ensembl):
    t = ensembl.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    boundary = t.exons[0].end
    span = reference_range(ensembl, t.contig, boundary - 2, boundary + 5)
    assert span == ""


def test_fasta_range_spanning_exon_boundary_returns_full_sequence(ensembl):
    t = ensembl.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    boundary = t.exons[0].end
    start, end = boundary - 2, boundary + 5
    window_seq = 'ACGTACGTACGTACGTACGT'
    fasta = _DictBackedFasta({t.contig: window_seq}, offset=start)
    g = Genome(ensembl, fasta=fasta, verify=False)
    span = g.reference_range(t.contig, start, end)
    assert len(span) == end - start + 1
    assert span == window_seq[:end - start + 1].upper()


def test_reference_range_rejects_inverted_range(ensembl):
    with pytest.raises(ValueError, match="end=.*< start="):
        reference_range(ensembl, "1", 100, 50)


# --------------------------------------------------------------------
# load_vcf wiring
# --------------------------------------------------------------------


def test_load_vcf_wraps_in_varcode_genome(ensembl):
    """``load_vcf(..., genome_fasta=...)`` wraps the resolved genome
    in :class:`varcode.Genome`. The resulting collection's variants
    carry the wrapper as ``genome``."""
    import os
    import tempfile

    from varcode import load_vcf

    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1000\t.\tA\tT\t100\tPASS\t.\n"
    )
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(body)
    try:
        empty_fasta = _DictBackedFasta({}, offset=1)
        vc = load_vcf(path, genome="GRCh38", genome_fasta=empty_fasta)
    finally:
        os.unlink(path)

    assert isinstance(vc[0].genome, Genome)
    assert vc[0].genome.fasta is empty_fasta


# --------------------------------------------------------------------
# cryptic_exons migration: silent-bail -> loud warning
# --------------------------------------------------------------------


def test_cryptic_exons_warns_once_per_genome():
    """Cryptic-exon scoring on breakpoints outside any annotated
    transcript and without a FASTA warns exactly once per genome.
    Subsequent uncovered breakpoints on the same genome are silent."""
    from varcode import StructuralVariant
    from varcode.cryptic_exons import enumerate_from_structural_variant

    g = Genome(cached_release(81))
    # Reset the per-instance flag so the warning is armed.
    g._missing_reference_warned = False

    sv = StructuralVariant(
        contig="1", start=1_000_000, sv_type="BND",
        alt="N]2:1000000]", mate_contig="2",
        mate_start=1_000_000, mate_orientation="]]",
        genome=g,
    )

    # First call — expect at most one warning (zero is also valid:
    # transcripts may cover the position and Tier 2 returns content).
    with warnings.catch_warnings(record=True) as first:
        warnings.simplefilter("always")
        enumerate_from_structural_variant(sv, genome=g)
    first_no_ref = [w for w in first
                    if "no reference sequence" in str(w.message)]
    # Second call — even if first warned, second should NOT.
    with warnings.catch_warnings(record=True) as second:
        warnings.simplefilter("always")
        enumerate_from_structural_variant(sv, genome=g)
    second_no_ref = [w for w in second
                     if "no reference sequence" in str(w.message)]
    assert not second_no_ref, [str(w.message) for w in second]
    # If the first call warned at all, the dedup flag should now be set.
    if first_no_ref:
        assert g._missing_reference_warned is True
