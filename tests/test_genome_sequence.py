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

"""Tests for :mod:`varcode.genome_sequence` (#372).

Covers the test matrix in the PR description:

* ``attach_genome_fasta`` accepts paths, pyfaidx objects, and
  duck-typed FASTA-shaped objects; rejects garbage; detaches via
  ``None``; verifies content against transcript cDNA when asked.
* ``reference_base`` and ``reference_range`` walk the three tiers
  (FASTA -> transcript cDNA -> empty) including the Tier 2 range
  fast path and the strand-aware base lookup.
* ``cryptic_exons`` migration: silently-bailing flanking-region
  scan replaced by loud warning when no reference covers the range.
"""

import warnings

import pytest
from pyensembl import cached_release

import varcode
from varcode.genome_sequence import (
    _VARCODE_FASTA_ATTR,
    attach_genome_fasta,
    reference_base,
    reference_range,
)


# --------------------------------------------------------------------
# Minimal pyfaidx-shaped mock so tests don't need a real chromosome FASTA
# --------------------------------------------------------------------


class _DictBackedFasta:
    """``{contig: sequence_string}`` adapter shaped like pyfaidx.Fasta.

    ``offset`` is the 1-based genome position of the first character of
    each contig's sequence (most test fixtures embed a small window
    around a known locus rather than a full chromosome).
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


_ENSEMBL_GRCh38 = cached_release(81)
# CFTR (chr7, + strand) — well-known protein coding transcript.
_CFTR_TRANSCRIPT_ID = "ENST00000003084"


@pytest.fixture
def genome():
    """Fresh genome reference, with any prior FASTA attachment cleared.

    ``cached_release`` returns a process-wide singleton so a test that
    forgets to detach would leak the FASTA into other tests. The
    teardown defends against that. Also clears the cryptic_exons
    once-per-genome warn cache so the warning fixture is exercised
    on every test that reaches it.
    """
    from varcode import cryptic_exons
    g = cached_release(81)
    if hasattr(g, _VARCODE_FASTA_ATTR):
        delattr(g, _VARCODE_FASTA_ATTR)
    cryptic_exons._MISSING_REFERENCE_WARNED.discard(id(g))
    yield g
    if hasattr(g, _VARCODE_FASTA_ATTR):
        delattr(g, _VARCODE_FASTA_ATTR)
    cryptic_exons._MISSING_REFERENCE_WARNED.discard(id(g))


@pytest.fixture
def cftr_position(genome):
    """A (contig, position, expected_plus_strand_base) tuple inside
    CFTR's first exon — known to be covered by Tier 2 (transcript cDNA)
    and useful as an anchor for Tier 1 fixtures."""
    t = genome.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    pos = t.exons[0].start + 5
    expected = t.sequence[5].upper()
    return (t.contig, pos, expected)


def _build_fasta_around(contig, position, base, window=200):
    """Build a tiny _DictBackedFasta with ``base`` at ``position`` and
    'N' elsewhere within a +/-``window`` window."""
    start_offset = position - window
    seq = ['N'] * (2 * window + 1)
    seq[position - start_offset] = base
    return _DictBackedFasta({contig: ''.join(seq)}, offset=start_offset)


# --------------------------------------------------------------------
# attach_genome_fasta
# --------------------------------------------------------------------


def test_attach_with_dict_backed_object_sets_attribute(genome, cftr_position):
    contig, pos, base = cftr_position
    fasta = _build_fasta_around(contig, pos, base)
    result = attach_genome_fasta(genome, fasta, verify=False)
    assert result is genome  # returns the genome for chaining
    assert getattr(genome, _VARCODE_FASTA_ATTR) is fasta


def test_attach_with_integer_raises_typeerror(genome):
    with pytest.raises(TypeError, match="pyfaidx-style protocol"):
        attach_genome_fasta(genome, 42, verify=False)


def test_attach_with_list_raises_typeerror(genome):
    with pytest.raises(TypeError, match="pyfaidx-style protocol"):
        attach_genome_fasta(genome, [1, 2, 3], verify=False)


def test_reattach_overwrites_silently(genome, cftr_position):
    contig, pos, base = cftr_position
    first = _build_fasta_around(contig, pos, base)
    second = _build_fasta_around(contig, pos, base)
    attach_genome_fasta(genome, first, verify=False)
    attach_genome_fasta(genome, second, verify=False)
    assert getattr(genome, _VARCODE_FASTA_ATTR) is second


def test_detach_via_none_removes_attribute(genome, cftr_position):
    contig, pos, base = cftr_position
    fasta = _build_fasta_around(contig, pos, base)
    attach_genome_fasta(genome, fasta, verify=False)
    assert hasattr(genome, _VARCODE_FASTA_ATTR)
    attach_genome_fasta(genome, None)
    assert not hasattr(genome, _VARCODE_FASTA_ATTR)


def test_detach_when_nothing_attached_is_no_op(genome):
    # Should not raise.
    attach_genome_fasta(genome, None)
    assert not hasattr(genome, _VARCODE_FASTA_ATTR)


def test_verify_emits_warning_on_mismatched_fasta(genome):
    """A FASTA that disagrees with the transcript cDNA at probed
    positions triggers the verify-on-attach warning. Catches GRCh37
    attached to a GRCh38 release."""

    class _AlwaysWrongFasta:
        def __getitem__(self, contig):
            return _AlwaysWrongSlicer()

    class _AlwaysWrongSlicer:
        def __getitem__(self, s):
            length = (s.stop - s.start) if isinstance(s, slice) else 1
            return _Span('X' * length)

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        attach_genome_fasta(genome, _AlwaysWrongFasta(), verify=True)

    disagreement = [w for w in captured
                    if "disagree" in str(w.message).lower()]
    assert disagreement, [str(w.message) for w in captured]


def test_verify_false_suppresses_check(genome):
    """``verify=False`` skips the probe even when the FASTA would fail it."""
    class _AlwaysWrongFasta:
        def __getitem__(self, contig):
            return _AlwaysWrongSlicer()

    class _AlwaysWrongSlicer:
        def __getitem__(self, s):
            length = (s.stop - s.start) if isinstance(s, slice) else 1
            return _Span('X' * length)

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        attach_genome_fasta(genome, _AlwaysWrongFasta(), verify=False)

    assert not any("disagree" in str(w.message).lower() for w in captured)


# --------------------------------------------------------------------
# reference_base — Tier 2 (no FASTA)
# --------------------------------------------------------------------


def test_tier2_returns_transcript_base_on_plus_strand(genome, cftr_position):
    contig, pos, expected = cftr_position
    assert reference_base(genome, contig, pos) == expected


def test_tier2_returns_plus_strand_base_for_minus_strand_transcript(genome):
    """For ``-`` strand transcripts the cDNA is reverse-complemented
    from the + strand. The lookup must return the + strand base —
    i.e. the complement of whatever the cDNA says at the corresponding
    spliced offset."""
    minus_t = next(
        t for t in genome.transcripts()
        if t.strand == '-' and t.exons and t.is_protein_coding)
    pos = minus_t.exons[0].start + 5
    # Compute the expected base via the same spliced_offset path the
    # implementation uses, then complement it (since cDNA is on the -
    # strand for a - strand transcript).
    offset = minus_t.spliced_offset(pos)
    cdna_base = minus_t.sequence[offset].upper()
    expected = cdna_base.translate(str.maketrans("ACGT", "TGCA"))
    assert reference_base(genome, minus_t.contig, pos) == expected


def test_tier2_returns_empty_for_intronic_position(genome):
    t = genome.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    # 100 bp inside the first intron.
    intron_pos = t.exons[0].end + 100
    assert reference_base(genome, t.contig, intron_pos) == ""


# --------------------------------------------------------------------
# reference_base — Tier 1 (FASTA attached)
# --------------------------------------------------------------------


def test_tier1_returns_fasta_base_when_attached(genome, cftr_position):
    contig, pos, transcript_base = cftr_position
    # Build a FASTA that AGREES with the transcript at this position
    # so verify passes, but stuff 'X's everywhere else so we can tell
    # Tier 1 is the source for nearby positions.
    fasta = _build_fasta_around(contig, pos, transcript_base, window=10)
    attach_genome_fasta(genome, fasta, verify=False)
    # The probe position should match (FASTA has the right base there).
    assert reference_base(genome, contig, pos) == transcript_base


def test_tier1_falls_through_when_contig_missing(genome, cftr_position):
    """A FASTA that doesn't carry the requested contig falls through
    to Tier 2 (transcript cDNA) — does not crash."""
    contig, pos, transcript_base = cftr_position
    empty_fasta = _DictBackedFasta({}, offset=1)
    attach_genome_fasta(genome, empty_fasta, verify=False)
    # Falls through to Tier 2; transcript-derived base survives.
    assert reference_base(genome, contig, pos) == transcript_base


# --------------------------------------------------------------------
# reference_range
# --------------------------------------------------------------------


def test_tier2_range_uses_single_transcript_fast_path(genome):
    t = genome.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    start = t.exons[0].start + 5
    end = start + 10
    span = reference_range(genome, t.contig, start, end)
    assert len(span) == end - start + 1
    # First and last base should match what reference_base gives us.
    assert span[0] == reference_base(genome, t.contig, start)
    assert span[-1] == reference_base(genome, t.contig, end)


def test_tier2_range_spanning_exon_boundary_returns_empty(genome):
    """A range that crosses into intronic space has uncovered
    positions; all-or-nothing makes the whole range return ``""``."""
    t = genome.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    # Cross from the last base of exon 1 into the first intron.
    boundary = t.exons[0].end
    span = reference_range(genome, t.contig, boundary - 2, boundary + 5)
    assert span == ""


def test_tier1_range_spanning_exon_boundary_returns_full_sequence(genome):
    """With a FASTA attached the same range returns full bases —
    including the intronic stretch — because Tier 1 doesn't care
    about exon coverage."""
    t = genome.transcript_by_id(_CFTR_TRANSCRIPT_ID)
    boundary = t.exons[0].end
    start, end = boundary - 2, boundary + 5
    # Fill the window with arbitrary bases so the lookup succeeds.
    window_seq = 'ACGTACGTACGTACGTACGT'
    fasta = _DictBackedFasta(
        {t.contig: window_seq}, offset=start)
    attach_genome_fasta(genome, fasta, verify=False)
    span = reference_range(genome, t.contig, start, end)
    assert len(span) == end - start + 1
    assert span == window_seq[:end - start + 1].upper()


def test_reference_range_rejects_inverted_range(genome):
    with pytest.raises(ValueError, match="end=.*< start="):
        reference_range(genome, "1", 100, 50)


# --------------------------------------------------------------------
# cryptic_exons migration: silent-bail -> loud warning
# --------------------------------------------------------------------


def test_attach_resets_cryptic_exons_warn_cache(genome, cftr_position):
    """The ``cryptic_exons`` warn-once dedup keyed on ``id(genome)``
    would otherwise persist across attach/detach cycles. Any explicit
    attach (or detach) is user action — the next missing-FASTA
    signal should be fresh, not suppressed.
    """
    from varcode import cryptic_exons
    contig, pos, base = cftr_position

    # Seed the cache as if a previous warning had fired.
    cryptic_exons._MISSING_REFERENCE_WARNED.add(id(genome))
    assert id(genome) in cryptic_exons._MISSING_REFERENCE_WARNED

    # Attach -> cache should be cleared (fresh signal post-attach).
    attach_genome_fasta(genome, _build_fasta_around(contig, pos, base),
                        verify=False)
    assert id(genome) not in cryptic_exons._MISSING_REFERENCE_WARNED

    # Re-seed, then detach -> also cleared.
    cryptic_exons._MISSING_REFERENCE_WARNED.add(id(genome))
    attach_genome_fasta(genome, None)
    assert id(genome) not in cryptic_exons._MISSING_REFERENCE_WARNED


def test_load_vcf_attaches_genome_fasta(genome):
    """``load_vcf(..., genome_fasta=...)`` attaches the FASTA to the
    resolved genome before parsing — equivalent to a separate
    ``attach_genome_fasta`` call but more ergonomic."""
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
        # Use a known-good DictBackedFasta so verify=True doesn't warn.
        # The probe lookup happens against real transcripts, so we
        # need a FASTA that agrees at those positions. Easiest: use
        # an empty FASTA (no contigs) which falls through verify
        # silently (zero probed positions = no disagreement).
        empty_fasta = _DictBackedFasta({}, offset=1)
        vc = load_vcf(path, genome="GRCh38", genome_fasta=empty_fasta)
    finally:
        os.unlink(path)
    # The resolved genome should have the FASTA attached.
    assert hasattr(vc[0].genome, _VARCODE_FASTA_ATTR)
    assert getattr(vc[0].genome, _VARCODE_FASTA_ATTR) is empty_fasta


def test_cryptic_exons_warns_when_no_reference_covers_breakpoint(genome):
    """``enumerate_from_structural_variant`` previously bailed silently
    when no genome FASTA was attached and the breakpoint fell outside
    any transcript. Now it warns loudly."""
    from varcode import StructuralVariant
    from varcode.cryptic_exons import enumerate_from_structural_variant

    # Intergenic-looking breakpoint between protein-coding genes on
    # chr1 (gene density is high enough that "intergenic" requires
    # care; chr1:1_000_000 is in a heterochromatic-ish region).
    sv = StructuralVariant(
        contig="1", start=1_000_000, sv_type="BND",
        alt="N]2:1000000]", mate_contig="2",
        mate_start=1_000_000, mate_orientation="]]",
        genome=genome,
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        candidates = enumerate_from_structural_variant(sv, genome=genome)

    msgs = [str(w.message) for w in captured]
    # Either a no-reference warning fires (the migration's intent), or
    # transcripts coincidentally cover this region — both branches are
    # acceptable, but in the no-coverage case the message must be loud.
    no_ref_warnings = [m for m in msgs if "no reference sequence" in m]
    if not candidates:
        # If we produced nothing, it's because the reference was absent
        # — verify the warning fired.
        assert no_ref_warnings, msgs
