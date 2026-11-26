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

"""Default DNA-only annotator, flat DNA hypotheses, and structural search anchors.

This module implements :class:`DNAEffectAnnotator` (registered as ``"dna"``,
the default annotator; also aliased as ``"unified"`` for 12.3.x compatibility):

* **One flat candidate shape** (:class:`DNAHypothesis`) across splicing,
  phasing, germline, and structural hypotheses — no branching between
  ``RealizedEffectCandidate`` and ``EffectCandidate`` and no 4-level unwrapping.
* **Explicit ascending genomic intron coordinates** (:class:`GenomicSpliceJunction`)
  on both ``+`` and ``-`` strands so downstream RNA tools (such as Isovar) can
  match BAM CIGAR ``N`` intervals without strand-flip bugs.
* **Deterministic DNA search anchors** (:class:`DNASearchAnchor` and
  :func:`build_dna_search_anchor`) for large insertions and structural
  rearrangements (including gene-to-intergenic translocations), giving Isovar
  the exact 5' anchor exon, exit reading frame, 3' re-entry exon or oriented
  mate search window, wild-type competing junctions, and mutant pre-mRNA
  template recipe without guessing deep/open-ended splicing from DNA alone.
* **Safe interior-exon coding fast-path** so simple SNVs, MNVs, small indels,
  non-coding, and incomplete transcripts avoid full genomic-layout construction
  while every splice-window, germline, phase, non-canonical, and structural case
  runs the full transcript model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Optional, Tuple
from weakref import WeakKeyDictionary

from serializable import DataclassSerializable

from .effect_candidates import EffectCandidate
from .effect_hypotheses import RealizedEffectCandidate
from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import STANDARD, codon_table_for_transcript
from .effects.effect_classes import (
    IncompleteTranscript,
    NoncodingTranscript,
    Silent,
)
from .effects.effect_ordering import classification_priority
from .effects.fast_path import try_fast_path_snv
from .effects.selenocysteine import reference_selenocysteine
from .effects.transcript_helpers import has_complete_coding_sequence
from .errors import ReferenceMismatchError
from .mutant_transcript import MutantTranscript
from .nucleotides import reverse_complement
from .structural_variant import affected_span
from .sv_allele_parser import breakend_sides
from .transcript_model import TranscriptModelEffectAnnotator
from .version import __version__ as _varcode_version


LARGE_INSERTION_THRESHOLD = 50
DEFAULT_MATE_SEARCH_WINDOW = 25_000


@dataclass(frozen=True)
class GenomicSpliceJunction(DataclassSerializable):
    """One exon-to-exon splice junction in explicit ascending genomic coordinates.

    Parameters
    ----------
    contig : str
        Normalized reference contig name.
    intron_start : int
        1-based inclusive genomic start of the excised intron
        (``intron_start <= intron_end`` on both ``+`` and ``-`` strands,
        matching BAM CIGAR ``N`` intervals).
    intron_end : int
        1-based inclusive genomic end of the excised intron.
    strand : str
        Transcript strand (``"+"`` or ``"-"``).
    donor_base : int
        1-based genomic coordinate of the last base of the 5' exon.
    acceptor_base : int
        1-based genomic coordinate of the first base of the 3' exon.
    """

    contig: str
    intron_start: int
    intron_end: int
    strand: str
    donor_base: int
    acceptor_base: int

    def __post_init__(self) -> None:
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")
        if self.intron_start > self.intron_end:
            raise ValueError(
                f"intron_start ({self.intron_start}) must be <= intron_end ({self.intron_end})"
            )

    @property
    def is_back_junction(self) -> bool:
        """True when the 5' donor lies downstream of the 3' acceptor (e.g. tandem DUP)."""
        if self.strand == "+":
            return self.donor_base > self.acceptor_base
        return self.donor_base < self.acceptor_base

    @classmethod
    def from_exon_join(
        cls,
        contig: str,
        donor_base: int,
        acceptor_base: int,
        strand: str,
    ) -> Optional["GenomicSpliceJunction"]:
        """Build a junction from 5'-exon donor and 3'-exon acceptor bases.

        Uses ascending 1-based genomic span ``(min + 1, max - 1)`` so both
        forward-spliced introns and tandem-duplication back-junctions (such as
        EGFR-KDD exon 25 -> exon 18) are preserved with ``intron_start <= intron_end``.
        Returns ``None`` when the two exon ends abut directly with no excised
        genomic span between them.
        """
        if strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {strand!r}")
        lo, hi = (donor_base, acceptor_base) if donor_base <= acceptor_base else (acceptor_base, donor_base)
        intron_start = lo + 1
        intron_end = hi - 1
        if intron_start > intron_end:
            return None
        return cls(
            contig=str(contig),
            intron_start=intron_start,
            intron_end=intron_end,
            strand=strand,
            donor_base=donor_base,
            acceptor_base=acceptor_base,
        )


@dataclass(frozen=True)
class DNAHypothesis(DataclassSerializable):
    """Flat, self-describing DNA-only candidate hypothesis for a variant effect.

    Replaces nested ``RealizedEffectCandidate -> ClassifiedOutcome ->
    EffectHypothesis -> SplicePlan`` traversal with one uniform record across
    point, splice, germline/phase, and structural consequences.
    """

    effect: Any
    source: str = "varcode.dna"
    uncertainty: Tuple[str, ...] = ()
    splice_mechanisms: Tuple[str, ...] = ()
    phase: Tuple[Tuple[Any, str], ...] = ()
    novel_junctions: Tuple[GenomicSpliceJunction, ...] = ()
    baseline_protein: Optional[str] = None
    mutant_protein: Optional[str] = None
    baseline_cdna: Optional[str] = None
    mutant_cdna: Optional[str] = None
    mutant_transcript: Optional[MutantTranscript] = None
    ordinal_rank: Tuple[int, ...] = ()
    probability: Optional[float] = None
    evidence: Mapping[str, Any] = field(default_factory=dict)

    @property
    def short_description(self) -> str:
        return self.effect.short_description


@dataclass(frozen=True)
class DNASearchAnchor(DataclassSerializable):
    """Deterministic DNA search geometry for resolving open-ended RNA splicing.

    For large insertions and structural rearrangements (including
    gene-to-intergenic translocations), DNA alone cannot predict deep cryptic
    splicing or polyadenylation. This record exposes the exact deterministic
    5' anchor exon, exit reading frame, 3' re-entry exon or oriented mate
    search window, competing wild-type junctions, and ordered mutant pre-mRNA
    template recipe so an upstream RNA tool (such as Isovar) knows where to
    harvest reads, what template to re-align against, and how to translate any
    observed spliced isoform.

    All genomic spans are 1-based inclusive ascending ``(start, end)``.
    """

    five_prime_transcript_id: str
    five_prime_gene_name: Optional[str]
    five_prime_contig: str
    five_prime_strand: str
    last_intact_exon_number: Optional[int]
    last_intact_exon_span: Optional[Tuple[int, int]]
    skip_fallback_exon_number: Optional[int]
    skip_fallback_exon_span: Optional[Tuple[int, int]]
    last_intact_donor_base: Optional[int]
    breakpoint_in_intron: bool
    five_prime_cds_prefix: Optional[str]
    exit_frame: Optional[int]
    three_prime_reentry_exon_number: Optional[int] = None
    three_prime_reentry_exon_span: Optional[Tuple[int, int]] = None
    three_prime_reentry_acceptor_base: Optional[int] = None
    three_prime_entry_frame: Optional[int] = None
    inserted_sequence: Optional[str] = None
    mate_search_region: Optional[Tuple[str, int, int, str]] = None
    wildtype_junctions: Tuple[GenomicSpliceJunction, ...] = ()
    template_segments: Tuple[
        Tuple[str, Optional[str], Optional[int], Optional[int], str], ...
    ] = ()

    def materialize_template(
        self,
        sequence_provider: Callable[[str, int, int], str],
    ) -> str:
        """Build the oriented 5'-to-3' mutant pre-mRNA template string.

        Parameters
        ----------
        sequence_provider : callable
            ``(contig, start, end) -> str`` returning 1-based inclusive
            plus-strand genomic sequence.

        Returns
        -------
        str
            Concatenated 5'-to-3' template sequence combining host genomic
            segments, any literal inserted DNA, and any oriented mate window.
        """
        pieces: list[str] = []
        for label, contig, start, end, strand in self.template_segments:
            if label == "inserted_sequence":
                if self.inserted_sequence:
                    pieces.append(self.inserted_sequence)
                continue
            if contig is None or start is None or end is None:
                continue
            seq = sequence_provider(contig, start, end)
            if strand == "-":
                seq = reverse_complement(seq)
            pieces.append(seq)
        return "".join(pieces)


def _exon_donor_and_acceptor(exon: Any, strand: str) -> Tuple[int, int]:
    """Return ``(acceptor_base, donor_base)`` in 1-based genomic coordinates."""
    if strand == "+":
        return exon.start, exon.end
    return exon.end, exon.start


def _cds_bounds(transcript: Any) -> Optional[Tuple[int, int]]:
    """Return inclusive spliced ``(cds_start, cds_end)``, tolerating poly-A stops."""
    if not getattr(transcript, "is_protein_coding", False):
        return None
    if not has_complete_coding_sequence(transcript):
        return None
    try:
        cds_start = min(transcript.start_codon_spliced_offsets)
    except (AttributeError, ValueError):
        return None
    try:
        cds_end = max(transcript.stop_codon_spliced_offsets)
    except (AttributeError, ValueError):
        cds_seq = getattr(transcript, "coding_sequence", None)
        if not cds_seq:
            return None
        cds_end = cds_start + len(cds_seq) - 1
    return cds_start, cds_end


def _cds_prefix_up_to(transcript: Any, genomic_pos: int) -> Optional[str]:
    """Return 5' CDS nucleotides up to and including ``genomic_pos``, or ``None``."""
    bounds = _cds_bounds(transcript)
    if bounds is None:
        return None
    cds_start, cds_end = bounds
    try:
        offset = transcript.spliced_offset(genomic_pos)
    except (KeyError, ValueError):
        return None
    if offset < cds_start:
        return None
    stop_limit = min(offset + 1, cds_end + 1)
    return str(transcript.sequence[cds_start:stop_limit])


def _cds_entry_frame_at(transcript: Any, genomic_pos: int) -> Optional[int]:
    """Return CDS offset mod 3 entering ``genomic_pos``, or ``None`` outside CDS."""
    bounds = _cds_bounds(transcript)
    if bounds is None:
        return None
    cds_start, cds_end = bounds
    try:
        offset = transcript.spliced_offset(genomic_pos)
    except (KeyError, ValueError):
        return None
    if offset < cds_start or offset > cds_end:
        return None
    return (offset - cds_start) % 3


def _canonical_junctions_around(
    transcript: Any,
    exon_numbers: Tuple[int, ...],
) -> Tuple[GenomicSpliceJunction, ...]:
    strand = "-" if transcript.on_backward_strand else "+"
    exons = transcript.exons
    seen: set[Tuple[str, int, int]] = set()
    junctions: list[GenomicSpliceJunction] = []
    for num in exon_numbers:
        for left_idx in (num - 2, num - 1, num):
            if 0 <= left_idx < len(exons) - 1:
                _, donor = _exon_donor_and_acceptor(exons[left_idx], strand)
                acceptor, _ = _exon_donor_and_acceptor(exons[left_idx + 1], strand)
                j = GenomicSpliceJunction.from_exon_join(
                    transcript.contig, donor, acceptor, strand
                )
                if j is not None:
                    key = (j.contig, j.intron_start, j.intron_end)
                    if key not in seen:
                        seen.add(key)
                        junctions.append(j)
    return tuple(junctions)


def _oriented_inserted_sequence(variant: Any, strand: str) -> Optional[str]:
    """Return 5'-to-3' transcript-oriented inserted DNA if present."""
    raw: Optional[str] = None
    if getattr(variant, "alt_assembly", None):
        raw = str(variant.alt_assembly)
    elif getattr(variant, "is_structural", False):
        symbolic = getattr(variant, "symbolic_alt", "") or ""
        if "[" in symbolic or "]" in symbolic:
            # BND alt like NACG[9:123[ or ]9:123]ACGN — strip anchor base.
            for sep in ("[", "]"):
                if sep in symbolic:
                    parts = symbolic.split(sep)
                    lead, trail = parts[0], parts[-1]
                    if lead and len(lead) > 1:
                        raw = lead[1:]
                    elif trail and len(trail) > 1:
                        raw = trail[:-1]
                    break
    else:
        alt = getattr(variant, "trimmed_alt", "")
        ref = getattr(variant, "trimmed_ref", "")
        if len(alt) > len(ref):
            raw = str(alt)
    if not raw:
        return None
    return reverse_complement(raw) if strand == "-" else raw


def build_dna_search_anchor(
    variant: Any,
    transcript: Any,
    *,
    mate_window: int = DEFAULT_MATE_SEARCH_WINDOW,
) -> DNASearchAnchor:
    """Build deterministic DNA search geometry for ``(variant, transcript)``.

    Works for both small/large explicit insertions (``Variant``) and structural
    variants (``StructuralVariant``: ``INS``, ``DEL``, ``DUP``, ``INV``, ``BND``)
    without guessing splicing inside inserted or intergenic sequence.
    """
    strand = "-" if transcript.on_backward_strand else "+"
    exons = transcript.exons
    is_sv = bool(getattr(variant, "is_structural", False))
    sv_type = getattr(variant, "sv_type", None)

    if is_sv:
        bp_low, bp_high = affected_span(variant)
    else:
        bp_low = variant.trimmed_base1_start
        bp_high = variant.trimmed_base1_end

    # Breakpoint entering from the 5' side of the transcript:
    # on '+' strand that's bp_low; on '-' strand that's bp_high.
    five_prime_bp = bp_high if strand == "-" else bp_low
    three_prime_bp = bp_low if strand == "-" else bp_high

    last_intact_idx: Optional[int] = None
    exonic_hit_idx: Optional[int] = None
    for idx, ex in enumerate(exons):
        if ex.start <= five_prime_bp <= ex.end:
            exonic_hit_idx = idx
            last_intact_idx = idx
            break
        if strand == "+":
            if ex.end < five_prime_bp:
                last_intact_idx = idx
        else:
            if ex.start > five_prime_bp:
                last_intact_idx = idx

    breakpoint_in_intron = exonic_hit_idx is None
    last_intact_num = (last_intact_idx + 1) if last_intact_idx is not None else None
    last_intact_span = (
        (exons[last_intact_idx].start, exons[last_intact_idx].end)
        if last_intact_idx is not None
        else None
    )

    fallback_idx = (
        last_intact_idx - 1
        if last_intact_idx is not None and last_intact_idx > 0
        else None
    )
    fallback_num = (fallback_idx + 1) if fallback_idx is not None else None
    fallback_span = (
        (exons[fallback_idx].start, exons[fallback_idx].end)
        if fallback_idx is not None
        else None
    )

    if last_intact_idx is None:
        last_donor: Optional[int] = None
        cds_prefix: Optional[str] = None
    elif breakpoint_in_intron:
        _, last_donor = _exon_donor_and_acceptor(exons[last_intact_idx], strand)
        cds_prefix = _cds_prefix_up_to(transcript, last_donor)
    else:
        last_donor = five_prime_bp
        cds_prefix = _cds_prefix_up_to(transcript, five_prime_bp)

    exit_frame = len(cds_prefix) % 3 if cds_prefix is not None else None

    # 3' re-entry exon on the host transcript (for intragenic insertions/local SVs,
    # not inter-locus BND translocations).
    is_bnd = sv_type == "BND" or bool(getattr(variant, "mate_contig", None))
    reentry_idx: Optional[int] = None
    if not is_bnd:
        for idx, ex in enumerate(exons):
            if strand == "+":
                if ex.start > three_prime_bp:
                    reentry_idx = idx
                    break
            else:
                if ex.end < three_prime_bp:
                    reentry_idx = idx
                    break

    reentry_num = (reentry_idx + 1) if reentry_idx is not None else None
    reentry_span = (
        (exons[reentry_idx].start, exons[reentry_idx].end)
        if reentry_idx is not None
        else None
    )
    if reentry_idx is not None:
        reentry_acceptor, _ = _exon_donor_and_acceptor(exons[reentry_idx], strand)
        reentry_frame = _cds_entry_frame_at(transcript, reentry_acceptor)
    else:
        reentry_acceptor = None
        reentry_frame = None

    inserted_seq = _oriented_inserted_sequence(variant, strand)

    # Oriented mate search region for BND / translocation events.
    mate_region: Optional[Tuple[str, int, int, str]] = None
    if is_bnd and getattr(variant, "mate_contig", None) and getattr(variant, "mate_start", None):
        symbolic = getattr(variant, "symbolic_alt", "") or ""
        mate_strand = "+"
        if "[" in symbolic or "]" in symbolic:
            _, mate_side = breakend_sides(symbolic)
            mate_strand = "+" if mate_side == "right" else "-"
        mate_pos = int(variant.mate_start)
        if mate_strand == "+":
            mate_region = (
                str(variant.mate_contig),
                mate_pos,
                mate_pos + max(1, mate_window) - 1,
                "+",
            )
        else:
            mate_region = (
                str(variant.mate_contig),
                max(1, mate_pos - max(1, mate_window) + 1),
                mate_pos,
                "-",
            )

    # Build ordered mutant pre-mRNA template recipe.
    template_segments: list[
        Tuple[str, Optional[str], Optional[int], Optional[int], str]
    ] = []
    if last_intact_idx is not None:
        anchor_start_exon = exons[fallback_idx] if fallback_idx is not None else exons[last_intact_idx]
        if strand == "+":
            up_start = anchor_start_exon.start
            up_end = five_prime_bp
        else:
            up_start = five_prime_bp
            up_end = anchor_start_exon.end
        if up_start <= up_end:
            template_segments.append(
                ("five_prime_host", str(transcript.contig), up_start, up_end, strand)
            )
    if inserted_seq:
        template_segments.append(("inserted_sequence", None, None, None, strand))
    if mate_region is not None:
        template_segments.append(
            ("mate_search_window", mate_region[0], mate_region[1], mate_region[2], mate_region[3])
        )
    elif reentry_idx is not None:
        end_exon = exons[min(len(exons) - 1, reentry_idx + 1)]
        if strand == "+":
            down_start = three_prime_bp + (0 if not is_sv else 1)
            down_end = end_exon.end
        else:
            down_start = end_exon.start
            down_end = max(1, three_prime_bp - (0 if not is_sv else 1))
        if down_start <= down_end:
            template_segments.append(
                ("three_prime_host", str(transcript.contig), down_start, down_end, strand)
            )

    anchor_exon_nums = tuple(
        n for n in (last_intact_num, reentry_num) if n is not None
    )
    wt_junctions = _canonical_junctions_around(transcript, anchor_exon_nums)

    return DNASearchAnchor(
        five_prime_transcript_id=str(transcript.id),
        five_prime_gene_name=getattr(transcript, "gene_name", None),
        five_prime_contig=str(transcript.contig),
        five_prime_strand=strand,
        last_intact_exon_number=last_intact_num,
        last_intact_exon_span=last_intact_span,
        skip_fallback_exon_number=fallback_num,
        skip_fallback_exon_span=fallback_span,
        last_intact_donor_base=last_donor,
        breakpoint_in_intron=breakpoint_in_intron,
        five_prime_cds_prefix=cds_prefix,
        exit_frame=exit_frame,
        three_prime_reentry_exon_number=reentry_num,
        three_prime_reentry_exon_span=reentry_span,
        three_prime_reentry_acceptor_base=reentry_acceptor,
        three_prime_entry_frame=reentry_frame,
        inserted_sequence=inserted_seq,
        mate_search_region=mate_region,
        wildtype_junctions=wt_junctions,
        template_segments=tuple(template_segments),
    )


def _needs_search_anchor(variant: Any) -> bool:
    if getattr(variant, "is_structural", False):
        return getattr(variant, "sv_type", None) in ("INS", "DEL", "DUP", "INV", "BND")
    if getattr(variant, "alt_assembly", None):
        return True
    alt = getattr(variant, "trimmed_alt", "")
    ref = getattr(variant, "trimmed_ref", "")
    return (len(alt) - len(ref)) >= LARGE_INSERTION_THRESHOLD


def _reference_exon_joins(transcript: Any) -> set[Tuple[str, int, int]]:
    """Canonical ``(contig, donor_base, acceptor_base)`` pairs for ``transcript``."""
    strand = "-" if getattr(transcript, "on_backward_strand", False) else "+"
    exons = getattr(transcript, "exons", ())
    contig = str(getattr(transcript, "contig", ""))
    joins: set[Tuple[str, int, int]] = set()
    for left, right in zip(exons, exons[1:]):
        _, donor = _exon_donor_and_acceptor(left, strand)
        acceptor, _ = _exon_donor_and_acceptor(right, strand)
        joins.add((contig, donor, acceptor))
    return joins


def _convert_junction_signature(
    junction_signature: Tuple[Tuple[Any, ...], ...],
    strand: str,
    reference_joins: Optional[set[Tuple[str, int, int]]] = None,
) -> Tuple[GenomicSpliceJunction, ...]:
    result: list[GenomicSpliceJunction] = []
    for entry in junction_signature:
        if len(entry) != 4:
            continue
        c1, donor, c2, acceptor = entry
        if str(c1) != str(c2):
            continue
        contig = str(c1)
        donor_i, acceptor_i = int(donor), int(acceptor)
        if reference_joins and (contig, donor_i, acceptor_i) in reference_joins:
            continue
        j = GenomicSpliceJunction.from_exon_join(contig, donor_i, acceptor_i, strand)
        if j is not None:
            result.append(j)
    return tuple(result)


def _infer_uncertainty_axes(effect: Any) -> Tuple[str, ...]:
    axes: set[str] = set()
    realized = getattr(effect, "realized_candidates", None)
    if realized is not None and len(realized) > 1:
        has_splice = False
        has_phase = False
        seen_phases: set[Tuple[Any, ...]] = set()
        seen_mechs: set[Tuple[str, ...]] = set()
        for cand in realized:
            for outcome in cand.outcomes:
                hyp = outcome.hypothesis
                phase_key = tuple(
                    (str(v), state) for v, state in hyp.phase
                )
                seen_phases.add(phase_key)
                if hyp.splice_plan is not None and hyp.splice_plan.choices:
                    mechs = tuple(opt.mechanism for _, opt in hyp.splice_plan.choices)
                    seen_mechs.add(mechs)
                    if any(m != "normal" for m in mechs):
                        has_splice = True
        if len(seen_phases) > 1:
            has_phase = True
        if len(seen_mechs) > 1 or has_splice:
            axes.add("splicing")
        if has_phase:
            axes.add("phasing")
    else:
        cands = getattr(effect, "candidates", ()) or ()
        if len(cands) > 1:
            if getattr(effect.variant, "is_structural", False):
                axes.add("structural")
            else:
                axes.add("splicing")
    return tuple(sorted(axes))


def _flatten_to_dna_hypotheses(effect: Any, transcript: Any) -> Tuple[DNAHypothesis, ...]:
    strand = "-" if getattr(transcript, "on_backward_strand", False) else "+"
    axes = _infer_uncertainty_axes(effect)
    realized = getattr(effect, "realized_candidates", None)

    if realized:
        ref_joins = _reference_exon_joins(transcript)
        hypotheses: list[DNAHypothesis] = []
        for cand in realized:
            assert isinstance(cand, RealizedEffectCandidate)
            mechs: list[str] = []
            phase_items: list[Tuple[Any, str]] = []
            seen_phase: set[Tuple[str, str]] = set()
            junctions: list[GenomicSpliceJunction] = []
            seen_j: set[Tuple[str, int, int, int, int]] = set()
            for outcome in cand.outcomes:
                hyp = outcome.hypothesis
                if hyp.splice_plan is not None:
                    for _, opt in hyp.splice_plan.choices:
                        if opt.mechanism not in mechs:
                            mechs.append(opt.mechanism)
                for v, state in hyp.phase:
                    pk = (str(v), state)
                    if pk not in seen_phase:
                        seen_phase.add(pk)
                        phase_items.append((v, state))
                for j in _convert_junction_signature(
                    outcome.mutant.junction_signature, strand, ref_joins
                ):
                    jk = (j.contig, j.intron_start, j.intron_end, j.donor_base, j.acceptor_base)
                    if jk not in seen_j:
                        seen_j.add(jk)
                        junctions.append(j)
            first = cand.outcomes[0]
            mt = first.mutant.to_mutant_transcript() if first.mutant is not None else None
            hyp_record = DNAHypothesis(
                effect=cand.effect,
                uncertainty=axes,
                splice_mechanisms=tuple(mechs) if mechs else ("normal",),
                phase=tuple(phase_items),
                novel_junctions=tuple(junctions),
                baseline_protein=first.baseline.protein_sequence,
                mutant_protein=first.mutant.protein_sequence,
                baseline_cdna=first.baseline.cdna_sequence,
                mutant_cdna=first.mutant.cdna_sequence,
                mutant_transcript=mt,
                ordinal_rank=cand.ordinal_key,
                probability=cand.probability,
                evidence={
                    "n_merged_outcomes": len(cand.outcomes),
                    "phase_mass": cand.phase_mass,
                },
            )
            # Non-field cache for in-memory proteome/ORF-discovery interop;
            # excluded from dataclass equality, repr, and JSON serialization.
            object.__setattr__(hyp_record, "_paired_outcomes", cand.outcomes)
            hypotheses.append(hyp_record)
        return tuple(hypotheses)

    legacy_cands = getattr(effect, "candidates", ()) or ()
    if legacy_cands:
        out: list[DNAHypothesis] = []
        for idx, cand in enumerate(legacy_cands):
            inner = cand.effect if isinstance(cand, EffectCandidate) else cand
            legacy_mt: Optional[MutantTranscript] = getattr(
                inner, "mutant_transcript", None
            )
            out.append(
                DNAHypothesis(
                    effect=inner,
                    source=cand.source if isinstance(cand, EffectCandidate) else "varcode.dna",
                    uncertainty=axes,
                    mutant_protein=getattr(inner, "mutant_protein_sequence", None)
                    or (legacy_mt.mutant_protein_sequence if legacy_mt is not None else None),
                    mutant_cdna=legacy_mt.cdna_sequence if legacy_mt is not None else None,
                    mutant_transcript=legacy_mt,
                    ordinal_rank=(idx,),
                    evidence=dict(cand.evidence) if isinstance(cand, EffectCandidate) else {},
                )
            )
        return tuple(out)

    bare_mt: Optional[MutantTranscript] = getattr(effect, "mutant_transcript", None)
    return (
        DNAHypothesis(
            effect=effect,
            uncertainty=axes,
            baseline_protein=getattr(effect, "original_protein_sequence", None),
            mutant_protein=getattr(effect, "mutant_protein_sequence", None)
            or (bare_mt.mutant_protein_sequence if bare_mt is not None else None),
            mutant_cdna=bare_mt.cdna_sequence if bare_mt is not None else None,
            mutant_transcript=bare_mt,
        ),
    )


_TRANSCRIPT_FAST_PATH_CACHE: WeakKeyDictionary[Any, Any] = WeakKeyDictionary()
_STANDARD_FORWARD = STANDARD.forward_table
_STANDARD_STOPS = STANDARD.stop_codons
_ACGT = frozenset("ACGT")
_FAST_PATH_MAX_ALLELE_BP = 20


def _transcript_fast_path_info(
    transcript: Any,
) -> Optional[Tuple[str, Tuple[Tuple[int, int], ...], int, str, str]]:
    """Cached per-transcript safety check for the interior coding fast-path.

    Returns ``(strand, exon_bounds, cds_start, sequence, protein_sequence)``
    when ``transcript`` uses the standard genetic code, starts with ``ATG``,
    has no circular origin wrap, no selenocysteine, and no ambiguous ``X``
    residues; otherwise returns ``None``.
    """
    try:
        cached = _TRANSCRIPT_FAST_PATH_CACHE.get(transcript, _MISSING)
    except TypeError:
        cached = _MISSING
    if cached is not _MISSING:
        return cached

    result: Optional[Tuple[str, Tuple[Tuple[int, int], ...], int, str, str]] = None
    strand = "-" if transcript.on_backward_strand else "+"
    exons = transcript.exons
    wraps_origin = any(ex.start > ex.end for ex in exons)
    if not wraps_origin and len(exons) > 1:
        if strand == "+":
            wraps_origin = any(
                exons[i + 1].start < exons[i].start for i in range(len(exons) - 1)
            )
        else:
            wraps_origin = any(
                exons[i + 1].end > exons[i].end for i in range(len(exons) - 1)
            )
    if (
        not wraps_origin
        and not reference_selenocysteine(transcript)
        and codon_table_for_transcript(transcript) is STANDARD
    ):
        cds_start = min(transcript.start_codon_spliced_offsets)
        seq = str(transcript.sequence)
        if seq[cds_start : cds_start + 3].upper() == "ATG":
            ref_protein = str(transcript.protein_sequence)
            if "X" not in ref_protein:
                exon_bounds = tuple((int(ex.start), int(ex.end)) for ex in exons)
                result = (strand, exon_bounds, cds_start, seq, ref_protein)

    try:
        _TRANSCRIPT_FAST_PATH_CACHE[transcript] = result
    except TypeError:
        pass
    return result


_MISSING = object()


def _wrap_unified_simple_effect(
    effect: Any,
    variant: Any,
    baseline_protein: Optional[str],
    mutant_protein: Optional[str],
) -> Any:
    effect.variants = (variant,)
    hypothesis = DNAHypothesis(
        effect=effect,
        uncertainty=(),
        splice_mechanisms=("normal",),
        baseline_protein=baseline_protein,
        mutant_protein=mutant_protein,
        ordinal_rank=(),
        probability=1.0,
    )
    effect.candidates = (hypothesis,)
    effect.most_likely_candidate = hypothesis
    effect.highest_priority_candidate = hypothesis
    effect.probability = 1.0
    effect.mechanism_rank = ()
    return effect


def _translate_standard_until_stop(dna: str) -> Optional[Tuple[str, bool]]:
    """Translate ``dna`` in-frame under the standard genetic code until a stop.

    Returns ``(protein, found_stop)`` or ``None`` if any codon is ambiguous.
    """
    aas: list[str] = []
    limit = len(dna) - (len(dna) % 3)
    for i in range(0, limit, 3):
        codon = dna[i : i + 3].upper()
        if codon in _STANDARD_STOPS:
            return "".join(aas), True
        aa = _STANDARD_FORWARD.get(codon)
        if aa is None:
            return None
        aas.append(aa)
    return "".join(aas), False


def _try_unified_fast_path(variant: Any, transcript: Any) -> Any:
    """Safe fast-path for noncoding, incomplete, and interior coding small variants.

    Handles:
    * ``NoncodingTranscript`` and ``IncompleteTranscript`` (matching
      ``transcript_model`` directly),
    * Interior-exon SNVs via single-codon lookup,
    * Interior-exon MNVs and small indels (<20 bp) strictly outside splice
      windows and start/stop codons via local translation +
      :func:`classify_from_protein_diff`.

    Returns ``None`` whenever any guard does not hold, falling through to the
    full transcript model.
    """
    if getattr(variant, "is_structural", False) or getattr(variant, "alt_assembly", None):
        return None
    if not getattr(transcript, "is_protein_coding", False):
        return _wrap_unified_simple_effect(
            NoncodingTranscript(variant, transcript), variant, None, None
        )
    if not has_complete_coding_sequence(transcript):
        return _wrap_unified_simple_effect(
            IncompleteTranscript(variant, transcript), variant, None, None
        )

    ref = getattr(variant, "trimmed_ref", "")
    alt = getattr(variant, "trimmed_alt", "")
    if len(ref) >= _FAST_PATH_MAX_ALLELE_BP or len(alt) >= _FAST_PATH_MAX_ALLELE_BP:
        return None
    if alt and not (_ACGT.issuperset(alt.upper())):
        return None

    info = _transcript_fast_path_info(transcript)
    if info is None:
        return None
    strand, exon_bounds, cds_start, seq, ref_protein = info

    v_start = variant.trimmed_base1_start
    v_end = variant.trimmed_base1_end
    orig_start = variant.start
    orig_end = variant.end
    g_start = min(v_start, orig_start)
    g_end = max(v_end, orig_end)
    last_idx = len(exon_bounds) - 1

    # Require the entire edit (including VCF padding anchor) to sit strictly
    # inside one exon and strictly outside canonical splice windows.
    hit_exon = False
    for idx, (ex_start, ex_end) in enumerate(exon_bounds):
        if ex_start <= g_start and g_end <= ex_end:
            if strand == "+":
                if idx > 0 and g_start <= ex_start + (1 if not ref else 0):
                    return None
                if idx < last_idx and g_end >= ex_end - (3 if not ref else 2):
                    return None
            else:
                if idx > 0 and g_end >= ex_end - (1 if not ref else 0):
                    return None
                if idx < last_idx and g_start <= ex_start + (3 if not ref else 2):
                    return None
            hit_exon = True
            break
    if not hit_exon:
        return None

    # Verify reference bases against the transcript (including anchor bases on
    # pure insertions/deletions) so ReferenceMismatchError is always raised.
    orig_ref = getattr(variant, "ref", "")
    if orig_ref:
        try:
            orig_cdna_start = transcript.spliced_offset(
                orig_end if strand == "-" else orig_start
            )
        except (KeyError, ValueError):
            return None
        orig_cdna_ref = (
            reverse_complement(orig_ref).upper() if strand == "-" else orig_ref.upper()
        )
        observed = seq[orig_cdna_start : orig_cdna_start + len(orig_cdna_ref)].upper()
        if observed != orig_cdna_ref:
            raise ReferenceMismatchError(
                variant,
                transcript,
                expected_ref=observed,
                observed_ref=orig_cdna_ref,
                transcript_offset=orig_cdna_start,
                genome_start=orig_start,
                genome_end=orig_end,
            )

    # Map the trimmed edit to 0-based half-open [cdna_start, cdna_end) in transcript cDNA.
    try:
        if not ref:
            anchor_off = transcript.spliced_offset(v_start)
            cdna_start = anchor_off + 1 if strand == "+" else anchor_off
            cdna_end = cdna_start
        else:
            cdna_start = transcript.spliced_offset(v_end if strand == "-" else v_start)
            cdna_end = cdna_start + len(ref)
    except (KeyError, ValueError):
        return None

    cds_edit_start = cdna_start - cds_start
    cds_edit_end = cdna_end - cds_start
    cds_len = len(ref_protein) * 3
    if cds_edit_start < 3 or cds_edit_end >= cds_len:
        return None

    cdna_ref = (reverse_complement(ref) if strand == "-" else ref).upper()
    cdna_alt = (reverse_complement(alt) if strand == "-" else alt).upper()

    # Case 1: O(1) single-codon SNV path.
    if len(cdna_ref) == 1 and len(cdna_alt) == 1:
        effect = try_fast_path_snv(
            variant=variant,
            transcript=transcript,
            trimmed_cdna_ref=cdna_ref,
            trimmed_cdna_alt=cdna_alt,
            sequence_from_start_codon=seq[cds_start:],
            cds_offset=cds_edit_start,
        )
        if effect is None:
            return None
        if isinstance(effect, Silent):
            effect.mutant_protein_sequence = ref_protein
        return _wrap_unified_simple_effect(
            effect, variant, ref_protein, effect.mutant_protein_sequence
        )

    # Case 2: Interior-exon MNV or small indel (<20 bp), translated locally
    # and classified through the exact same classify_from_protein_diff used by
    # transcript_model.
    first_codon = cds_edit_start // 3
    window_cds_start = first_codon * 3
    length_delta = len(cdna_alt) - len(cdna_ref)

    if length_delta % 3 == 0:
        window_cds_end = ((cds_edit_end + 2) // 3) * 3
        window_dna = (
            seq[cds_start + window_cds_start : cdna_start]
            + cdna_alt
            + seq[cdna_end : cds_start + window_cds_end]
        )
        translated = _translate_standard_until_stop(window_dna)
        if translated is None:
            return None
        window_aa, hit_stop = translated
        if hit_stop:
            # Pure codon-aligned insertions creating an immediate stop use edit
            # provenance in classify_from_protein_diff to keep aa_ref=""; defer
            # this rare case to _full for exact parity.
            if not cdna_ref:
                return None
            mut_protein = ref_protein[:first_codon] + window_aa
        else:
            mut_protein = (
                ref_protein[:first_codon] + window_aa + ref_protein[window_cds_end // 3 :]
            )
        if mut_protein == ref_protein:
            # Defer rare synonymous MNVs/indels so _unchanged_product_effect
            # computes exact aa_pos/aa_ref.
            return None
        effect = classify_from_protein_diff(
            variant=variant,
            transcript=transcript,
            ref_protein=ref_protein,
            mut_protein=mut_protein,
            length_delta=0,
            frame_change_aa_offset=first_codon,
            stop_observed=True,
        )
    else:
        tail_dna = (
            seq[cds_start + window_cds_start : cdna_start] + cdna_alt + seq[cdna_end:]
        )
        translated = _translate_standard_until_stop(tail_dna)
        if translated is None:
            return None
        tail_aa, found_stop = translated
        mut_protein = ref_protein[:first_codon] + tail_aa
        if mut_protein == ref_protein:
            return None
        effect = classify_from_protein_diff(
            variant=variant,
            transcript=transcript,
            ref_protein=ref_protein,
            mut_protein=mut_protein,
            length_delta=length_delta,
            frame_change_aa_offset=first_codon,
            stop_observed=found_stop,
        )

    return _wrap_unified_simple_effect(
        effect, variant, ref_protein, effect.mutant_protein_sequence
    )


def finalize_dna_effect(effect: Any, variant: Any, transcript: Any) -> Any:
    """Attach flat :class:`DNAHypothesis` candidates and optional search anchor."""
    from .effects.effect_classes import MultiOutcomeEffect

    if effect is NotImplemented:
        return effect
    hypotheses = _flatten_to_dna_hypotheses(effect, transcript)
    effect.candidates = hypotheses
    if not isinstance(effect, MultiOutcomeEffect):
        effect.most_likely_candidate = hypotheses[0]
        effect.highest_priority_candidate = max(
            hypotheses, key=lambda h: classification_priority(h.effect)
        )
    # Drop nested realized_candidates so iter_candidate_paths and
    # get_effect_candidates yield the flat DNAHypothesis tuple uniformly.
    if "realized_candidates" in effect.__dict__:
        del effect.__dict__["realized_candidates"]
    if variant is not None and _needs_search_anchor(variant):
        effect.dna_search_anchor = build_dna_search_anchor(variant, transcript)
    return effect


class DNAEffectAnnotator:
    """Default DNA-only annotator combining a safe coding fast-path with flat hypotheses.

    Registered as ``"dna"`` (also aliased as ``"unified"`` for 12.3.x
    compatibility). Keeps ``fast`` and ``transcript_model`` intact as explicit
    opt-in annotators.
    """

    name = "dna"
    version = _varcode_version

    def __init__(self) -> None:
        self._full = TranscriptModelEffectAnnotator()

    def _finalize(self, effect: Any, variant: Any, transcript: Any) -> Any:
        return finalize_dna_effect(effect, variant, transcript)

    def annotate_on_transcript(self, variant: Any, transcript: Any) -> Any:
        fast_result = _try_unified_fast_path(variant, transcript)
        if fast_result is not None:
            return fast_result
        result = self._full.annotate_on_transcript(variant, transcript)
        return self._finalize(result, variant, transcript)

    def annotate_with_context(
        self,
        variant: Any,
        transcript: Any,
        germline_ctx: Any,
        phase_resolver: Any = None,
    ) -> Any:
        result = self._full.annotate_with_context(
            variant, transcript, germline_ctx, phase_resolver
        )
        return self._finalize(result, variant, transcript)

    def annotate_haplotype(
        self,
        variants: Any,
        transcript: Any,
        germline_ctx: Any = None,
        phase_resolver: Any = None,
    ) -> Any:
        result = self._full.annotate_haplotype(
            variants, transcript, germline_ctx, phase_resolver
        )
        primary = variants[0] if variants else None
        return self._finalize(result, primary, transcript)


# Backwards-compatible alias for varcode 12.3.x callers.
UnifiedDNAEffectAnnotator = DNAEffectAnnotator
