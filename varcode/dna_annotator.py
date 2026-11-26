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

"""Default DNA-only annotator returning ``MutationEffect``, ``EffectCandidate``, and ``MutantTranscript``.

:class:`DNAEffectAnnotator` (registered as ``"dna"``, the default annotator):

* **Guarded interior-exon fast-path** for non-coding transcripts, incomplete
  transcripts, interior coding SNVs, and interior coding MNVs / small indels
  (``< 20`` bp) strictly outside splice windows and start/stop codons.
* **One flat candidate shape** (:class:`~varcode.EffectCandidate`) across
  point, splice, germline/phase, and structural consequences, carrying
  ``mechanisms``, ``phase``, ``uncertainty``, and ``probability``.
* **Ascending genomic splice junctions** on :class:`~varcode.MutantTranscript`
  (``novel_junctions`` and ``back_junctions`` as 1-based inclusive
  ``(contig, start, end, strand)`` tuples with ``start <= end`` on both
  ``"+"`` and ``"-"`` strands, directly matching BAM CIGAR ``N`` intervals).
* **Open rearrangement templates** on :class:`~varcode.MutantTranscript`
  (``exit_frame``, ``template_sequence(sequence_provider)``, and breakpoint
  anchor metadata in ``evidence``) for large insertions (``>= 50`` bp) and
  structural rearrangements without a supplied ``alt_assembly``.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Optional, Tuple
from weakref import WeakKeyDictionary

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
DEFAULT_MATE_SEARCH_WINDOW_BP = 25_000
_FAST_PATH_MAX_ALLELE_BP = 20

JunctionTuple = Tuple[str, int, int, str]
TemplateSegment = Tuple[str, str, int, int, str]


# ---------------------------------------------------------------------------
# Junction & open-template geometry helpers
# ---------------------------------------------------------------------------


def _exon_donor_and_acceptor(exon: Any, strand: str) -> Tuple[int, int]:
    """Return ``(acceptor_base, donor_base)`` in 1-based genomic coordinates."""
    if strand == "+":
        return int(exon.start), int(exon.end)
    return int(exon.end), int(exon.start)


def _exon_join_to_junction(
    contig: str,
    donor_base: int,
    acceptor_base: int,
    strand: str,
) -> Optional[Tuple[str, int, int, str, bool]]:
    """Convert exon-edge coordinates to ``(contig, start, end, strand, is_back)``.

    Coordinates are always 1-based inclusive ascending (``start <= end``) on
    both ``"+"`` and ``"-"`` strands.
    """
    is_back = (
        (strand == "+" and donor_base > acceptor_base)
        or (strand == "-" and donor_base < acceptor_base)
    )
    if is_back:
        start = min(donor_base, acceptor_base) + 1
        end = max(donor_base, acceptor_base) - 1
    elif strand == "+":
        start = donor_base + 1
        end = acceptor_base - 1
    else:
        start = acceptor_base + 1
        end = donor_base - 1
    if start > end:
        return None
    return (str(contig), int(start), int(end), strand, is_back)


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


def _extract_novel_and_back_junctions(
    junction_signature: Tuple[Tuple[Any, ...], ...],
    strand: str,
    reference_joins: set[Tuple[str, int, int]],
) -> Tuple[Tuple[JunctionTuple, ...], Tuple[JunctionTuple, ...]]:
    novel: list[JunctionTuple] = []
    back: list[JunctionTuple] = []
    seen_novel: set[JunctionTuple] = set()
    seen_back: set[JunctionTuple] = set()
    for entry in junction_signature:
        if len(entry) != 4:
            continue
        c1, donor, c2, acceptor = entry
        if str(c1) != str(c2):
            continue
        contig = str(c1)
        donor_i, acceptor_i = int(donor), int(acceptor)
        if (contig, donor_i, acceptor_i) in reference_joins:
            continue
        converted = _exon_join_to_junction(contig, donor_i, acceptor_i, strand)
        if converted is None:
            continue
        c, s, e, st, is_back = converted
        jt: JunctionTuple = (c, s, e, st)
        if is_back:
            if jt not in seen_back:
                seen_back.add(jt)
                back.append(jt)
        else:
            if jt not in seen_novel:
                seen_novel.add(jt)
                novel.append(jt)
    return tuple(novel), tuple(back)


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
) -> Tuple[JunctionTuple, ...]:
    strand = "-" if transcript.on_backward_strand else "+"
    exons = transcript.exons
    seen: set[JunctionTuple] = set()
    junctions: list[JunctionTuple] = []
    for num in exon_numbers:
        for left_idx in (num - 2, num - 1, num):
            if 0 <= left_idx < len(exons) - 1:
                _, donor = _exon_donor_and_acceptor(exons[left_idx], strand)
                acceptor, _ = _exon_donor_and_acceptor(exons[left_idx + 1], strand)
                converted = _exon_join_to_junction(
                    transcript.contig, donor, acceptor, strand
                )
                if converted is not None:
                    c, s, e, st, _ = converted
                    jt: JunctionTuple = (c, s, e, st)
                    if jt not in seen:
                        seen.add(jt)
                        junctions.append(jt)
    return tuple(junctions)


def _needs_open_template(variant: Any) -> bool:
    if getattr(variant, "alt_assembly", None):
        return False
    if getattr(variant, "is_structural", False):
        return getattr(variant, "sv_type", None) in ("INS", "DEL", "DUP", "INV", "BND")
    alt = getattr(variant, "trimmed_alt", "")
    ref = getattr(variant, "trimmed_ref", "")
    return (len(alt) - len(ref)) >= LARGE_INSERTION_THRESHOLD


def _build_open_template_evidence(
    variant: Any,
    transcript: Any,
    *,
    mate_window_bp: int = DEFAULT_MATE_SEARCH_WINDOW_BP,
) -> dict[str, Any]:
    """Compute deterministic search-template metadata for an SV or large insertion."""
    strand = "-" if transcript.on_backward_strand else "+"
    exons = transcript.exons
    contig = str(transcript.contig)

    sv_type = getattr(variant, "sv_type", None)
    is_bnd = sv_type == "BND" or getattr(variant, "mate_contig", None) is not None

    if getattr(variant, "is_structural", False):
        span_start, span_end = affected_span(variant)
        if span_end is None:
            span_end = span_start
    else:
        span_start = variant.trimmed_base1_start
        span_end = variant.trimmed_base1_end

    if strand == "+":
        five_prime_bp = span_start
        three_prime_bp = span_end
    else:
        five_prime_bp = span_end
        three_prime_bp = span_start

    # Locate last intact 5' exon.
    last_intact_idx: Optional[int] = None
    breakpoint_in_intron = True
    for idx, ex in enumerate(exons):
        if ex.start <= five_prime_bp <= ex.end:
            last_intact_idx = idx
            breakpoint_in_intron = False
            break
        if strand == "+" and ex.end < five_prime_bp:
            last_intact_idx = idx
        elif strand == "-" and ex.start > five_prime_bp:
            last_intact_idx = idx

    last_intact_num = last_intact_idx + 1 if last_intact_idx is not None else None
    skip_fallback_num = (
        last_intact_num - 1
        if last_intact_num is not None and last_intact_num > 1
        else None
    )

    if last_intact_idx is not None:
        _, last_donor = _exon_donor_and_acceptor(exons[last_intact_idx], strand)
        cds_prefix = (
            _cds_prefix_up_to(transcript, last_donor)
            if breakpoint_in_intron
            else _cds_prefix_up_to(transcript, five_prime_bp)
        )
    else:
        cds_prefix = None
    exit_frame = len(cds_prefix) % 3 if cds_prefix is not None else None

    # Locate 3' re-entry exon on the same transcript (not applicable for BND).
    three_prime_reentry_idx: Optional[int] = None
    if not is_bnd:
        for idx, ex in enumerate(exons):
            if strand == "+":
                if ex.start > three_prime_bp:
                    three_prime_reentry_idx = idx
                    break
                if ex.start <= three_prime_bp <= ex.end:
                    three_prime_reentry_idx = idx
                    break
            else:
                if ex.end < three_prime_bp:
                    three_prime_reentry_idx = idx
                    break
                if ex.start <= three_prime_bp <= ex.end:
                    three_prime_reentry_idx = idx
                    break

    three_prime_reentry_num = (
        three_prime_reentry_idx + 1 if three_prime_reentry_idx is not None else None
    )
    three_prime_entry_frame: Optional[int] = None
    if three_prime_reentry_idx is not None:
        reentry_acceptor, _ = _exon_donor_and_acceptor(
            exons[three_prime_reentry_idx], strand
        )
        three_prime_entry_frame = _cds_entry_frame_at(transcript, reentry_acceptor)

    # Inserted sequence (literal bases on point insertion or SV alt_assembly/info).
    inserted_seq: Optional[str] = None
    if not getattr(variant, "is_structural", False):
        if len(variant.trimmed_alt) > len(variant.trimmed_ref):
            inserted_seq = str(variant.trimmed_alt).upper()
    else:
        info = getattr(variant, "info", None) or {}
        if isinstance(info, dict) and isinstance(info.get("SVINSSEQ"), str):
            inserted_seq = info["SVINSSEQ"].upper()

    # Mate search region for BND / translocations.
    mate_search_region: Optional[JunctionTuple] = None
    if is_bnd and getattr(variant, "mate_contig", None) is not None:
        mate_contig = str(variant.mate_contig)
        mate_pos = int(variant.mate_start)
        mate_strand = "+"
        sym_alt = getattr(variant, "symbolic_alt", "") or ""
        sides = breakend_sides(sym_alt) if sym_alt else None
        if sides is not None:
            _, mate_side = sides
            mate_strand = "+" if mate_side == "right" else "-"
        if mate_strand == "+":
            m_start = mate_pos
            m_end = mate_pos + mate_window_bp - 1
        else:
            m_start = max(1, mate_pos - mate_window_bp + 1)
            m_end = mate_pos
        mate_search_region = (mate_contig, m_start, m_end, mate_strand)

    # Wild-type junctions competing with novel mutant splicing near breakpoint.
    relevant_exons = tuple(
        n
        for n in (last_intact_num, three_prime_reentry_num)
        if n is not None
    )
    wt_junctions = _canonical_junctions_around(transcript, relevant_exons)

    # Ordered 5'->3' mutant pre-mRNA template segments.
    template_segments: list[TemplateSegment] = []
    if last_intact_idx is not None:
        anchor_start_idx = max(0, last_intact_idx - 1)
        if strand == "+":
            host_s = int(exons[anchor_start_idx].start)
            host_e = int(five_prime_bp)
        else:
            host_s = int(five_prime_bp)
            host_e = int(exons[anchor_start_idx].end)
        if host_s <= host_e:
            template_segments.append(("five_prime_host", contig, host_s, host_e, strand))

    if inserted_seq:
        template_segments.append(
            ("inserted_sequence", contig, int(five_prime_bp), int(five_prime_bp), strand)
        )

    if mate_search_region is not None:
        mc, ms, me, mstrand = mate_search_region
        template_segments.append(("mate_search_window", mc, ms, me, mstrand))
    elif three_prime_reentry_idx is not None:
        reentry_end_idx = min(len(exons) - 1, three_prime_reentry_idx + 1)
        if strand == "+":
            tail_s = int(three_prime_bp) + (1 if not inserted_seq else 0)
            tail_e = int(exons[reentry_end_idx].end)
        else:
            tail_s = int(exons[reentry_end_idx].start)
            tail_e = int(three_prime_bp) - (1 if not inserted_seq else 0)
        if tail_s <= tail_e:
            template_segments.append(("three_prime_host", contig, tail_s, tail_e, strand))

    evidence: dict[str, Any] = {
        "last_intact_exon": last_intact_num,
        "skip_fallback_exon": skip_fallback_num,
        "three_prime_reentry_exon": three_prime_reentry_num,
        "exit_frame": exit_frame,
        "three_prime_entry_frame": three_prime_entry_frame,
        "breakpoint_in_intron": breakpoint_in_intron,
        "wildtype_junctions": wt_junctions,
        "mate_search_region": mate_search_region,
        "template_segments": tuple(template_segments),
    }
    if inserted_seq is not None:
        evidence["inserted_sequence"] = inserted_seq
    return evidence


def _is_literal_large_insertion(variant: Any) -> bool:
    if getattr(variant, "is_structural", False):
        if getattr(variant, "sv_type", None) != "INS":
            return False
        info = getattr(variant, "info", None) or {}
        return isinstance(info, dict) and isinstance(info.get("SVINSSEQ"), str) and bool(info["SVINSSEQ"])
    alt = getattr(variant, "trimmed_alt", "")
    ref = getattr(variant, "trimmed_ref", "")
    return (len(alt) - len(ref)) >= LARGE_INSERTION_THRESHOLD


def _attach_open_template(effect: Any, variant: Any, transcript: Any) -> None:
    if variant is None or transcript is None or not _needs_open_template(variant):
        return
    existing_mt = getattr(effect, "mutant_transcript", None)
    if existing_mt is not None:
        anchor_ev = _build_open_template_evidence(variant, transcript)
        merged = {**(existing_mt.evidence or {}), **anchor_ev}
        effect.mutant_transcript = replace(existing_mt, evidence=merged)
    elif _is_literal_large_insertion(variant):
        anchor_ev = _build_open_template_evidence(variant, transcript)
        anchor_ev.setdefault("sequence_status", "open_search_template")
        effect.mutant_transcript = MutantTranscript(
            reference_transcript=transcript,
            reference_segments=(),
            cdna_sequence=None,
            mutant_protein_sequence=None,
            annotator_name="dna",
            evidence=anchor_ev,
        )


# ---------------------------------------------------------------------------
# Candidate flattening onto EffectCandidate
# ---------------------------------------------------------------------------


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
                phase_key = tuple((str(v), state) for v, state in hyp.phase)
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


def _attach_junctions_to_effect(
    effect: Any,
    novel_junctions: Tuple[JunctionTuple, ...],
    back_junctions: Tuple[JunctionTuple, ...],
) -> None:
    if not novel_junctions and not back_junctions:
        return
    mt = getattr(effect, "mutant_transcript", None)
    if mt is None:
        return
    updated_ev = dict(mt.evidence or {})
    if novel_junctions:
        updated_ev["novel_junctions"] = novel_junctions
    if back_junctions:
        updated_ev["back_junctions"] = back_junctions
    effect.mutant_transcript = replace(mt, evidence=updated_ev)


def _flatten_to_candidates(effect: Any, transcript: Any) -> Tuple[EffectCandidate, ...]:
    strand = "-" if getattr(transcript, "on_backward_strand", False) else "+"
    axes = _infer_uncertainty_axes(effect)
    realized = getattr(effect, "realized_candidates", None)

    if realized:
        ref_joins = _reference_exon_joins(transcript)
        candidates: list[EffectCandidate] = []
        for cand in realized:
            assert isinstance(cand, RealizedEffectCandidate)
            mechs: list[str] = []
            phase_items: list[Tuple[Any, str]] = []
            seen_phase: set[Tuple[str, str]] = set()
            novel_list: list[JunctionTuple] = []
            back_list: list[JunctionTuple] = []
            seen_novel: set[JunctionTuple] = set()
            seen_back: set[JunctionTuple] = set()

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
                novel_j, back_j = _extract_novel_and_back_junctions(
                    outcome.mutant.junction_signature, strand, ref_joins
                )
                for jt in novel_j:
                    if jt not in seen_novel:
                        seen_novel.add(jt)
                        novel_list.append(jt)
                for jt in back_j:
                    if jt not in seen_back:
                        seen_back.add(jt)
                        back_list.append(jt)

            _attach_junctions_to_effect(cand.effect, tuple(novel_list), tuple(back_list))

            ec = EffectCandidate(
                effect=cand.effect,
                source="varcode.dna",
                evidence={
                    "n_merged_outcomes": len(cand.outcomes),
                    "phase_mass": cand.phase_mass,
                },
                mechanisms=tuple(mechs) if mechs else ("normal",),
                phase=tuple(phase_items),
                uncertainty=axes,
            )
            # Non-field cache for in-memory proteome / ORF-discovery / RNA-refinement
            # interop; excluded from dataclass equality, repr, and JSON serialization.
            object.__setattr__(ec, "_paired_outcomes", cand.outcomes)
            candidates.append(ec)
        return tuple(candidates)

    legacy_cands = getattr(effect, "candidates", ()) or ()
    if legacy_cands:
        out: list[EffectCandidate] = []
        for cand in legacy_cands:
            if isinstance(cand, EffectCandidate):
                out.append(
                    replace(
                        cand,
                        uncertainty=cand.uncertainty or axes,
                    )
                )
            else:
                out.append(
                    EffectCandidate(
                        effect=cand,
                        source="varcode.dna",
                        uncertainty=axes,
                    )
                )
        return tuple(out)

    return (
        EffectCandidate(
            effect=effect,
            source="varcode.dna",
            uncertainty=axes,
        ),
    )


# ---------------------------------------------------------------------------
# Guarded interior-exon fast-path
# ---------------------------------------------------------------------------


_TRANSCRIPT_FAST_PATH_CACHE: WeakKeyDictionary[Any, Any] = WeakKeyDictionary()
_STANDARD_FORWARD = STANDARD.forward_table
_STANDARD_STOPS = STANDARD.stop_codons
_ACGT = frozenset("ACGT")
_MISSING = object()


def _transcript_fast_path_info(
    transcript: Any,
) -> Optional[Tuple[str, Tuple[Tuple[int, int], ...], int, str, str]]:
    """Cached per-transcript safety check for the interior coding fast-path."""
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


def _wrap_simple_dna_effect(effect: Any, variant: Any) -> Any:
    effect.variants = (variant,)
    candidate = EffectCandidate(
        effect=effect,
        source="varcode.dna",
        mechanisms=("normal",),
        uncertainty=(),
    )
    effect.candidates = (candidate,)
    effect.most_likely_candidate = candidate
    effect.highest_priority_candidate = candidate
    effect.probability = 1.0
    effect.mechanism_rank = ()
    return effect


def _translate_standard_until_stop(dna: str) -> Optional[Tuple[str, bool]]:
    """Translate ``dna`` in-frame under the standard genetic code until a stop."""
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


def _try_dna_fast_path(variant: Any, transcript: Any) -> Any:
    """Safe fast-path for noncoding, incomplete, and interior coding small variants."""
    if getattr(variant, "is_structural", False) or getattr(variant, "alt_assembly", None):
        return None
    if not getattr(transcript, "is_protein_coding", False):
        return _wrap_simple_dna_effect(NoncodingTranscript(variant, transcript), variant)
    if not has_complete_coding_sequence(transcript):
        return _wrap_simple_dna_effect(IncompleteTranscript(variant, transcript), variant)

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
        return _wrap_simple_dna_effect(effect, variant)

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
            if not cdna_ref:
                return None
            mut_protein = ref_protein[:first_codon] + window_aa
        else:
            mut_protein = (
                ref_protein[:first_codon]
                + window_aa
                + ref_protein[window_cds_end // 3 :]
            )
        if mut_protein == ref_protein:
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

    return _wrap_simple_dna_effect(effect, variant)


# ---------------------------------------------------------------------------
# Public annotator
# ---------------------------------------------------------------------------


def finalize_dna_effect(effect: Any, variant: Any, transcript: Any) -> Any:
    """Attach flat :class:`~varcode.EffectCandidate` tuples and open template metadata."""
    from .effects.effect_classes import MultiOutcomeEffect

    if effect is NotImplemented:
        return effect
    _attach_open_template(effect, variant, transcript)
    candidates = _flatten_to_candidates(effect, transcript)
    effect.candidates = candidates
    if not isinstance(effect, MultiOutcomeEffect):
        effect.most_likely_candidate = candidates[0]
        effect.highest_priority_candidate = max(
            candidates, key=lambda c: classification_priority(c.effect)
        )
    if "realized_candidates" in effect.__dict__:
        del effect.__dict__["realized_candidates"]
    return effect


class DNAEffectAnnotator:
    """Default DNA-only annotator combining a safe coding fast-path with flat candidates.

    Registered as ``"dna"`` (the default annotator). Returns standard
    :class:`~varcode.MutationEffect`, :class:`~varcode.EffectCandidate`, and
    :class:`~varcode.MutantTranscript` objects across point, splice,
    germline/phase, and structural variants. Keeps ``fast`` and
    ``transcript_model`` intact as explicit opt-in annotators.
    """

    name = "dna"
    version = _varcode_version

    def __init__(self) -> None:
        self._full = TranscriptModelEffectAnnotator()

    def annotate_on_transcript(self, variant: Any, transcript: Any) -> Any:
        fast_result = _try_dna_fast_path(variant, transcript)
        if fast_result is not None:
            return fast_result
        result = self._full.annotate_on_transcript(variant, transcript)
        return finalize_dna_effect(result, variant, transcript)

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
        return finalize_dna_effect(result, variant, transcript)

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
        return finalize_dna_effect(result, primary, transcript)
