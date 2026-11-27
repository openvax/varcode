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


"""DNA-only splice mechanism ordering, normal-splicing, and exon-skipping builders."""

from .circular_origin import origin_partners
from .effect_candidates import EffectCandidate
from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import codon_table_for_transcript, translate_sequence
from .effects.effect_classes import (
    CrypticAcceptor,
    CrypticDonor,
    ExonSkipping,
    ExonicSpliceSite,
    IntronRetention,
    IntronicSpliceSite,
    NormalSplicing,
    SpliceAcceptor,
    SpliceDonor,
)
from .effects.effect_helpers import exon_length
from .effects.selenocysteine import edited_selenocysteine
from .effects.transcript_helpers import has_complete_coding_sequence, poly_a_completion
from .mutant_transcript import MutantTranscript, TranscriptEdit


# ---------------------------------------------------------------------
# DNA-only mechanism order.
#
# Hand-tuned heuristic ordering. These are not probabilities. The order
# only gives varcode a deterministic leading candidate when DNA alone
# leaves several splice mechanisms possible.
# ---------------------------------------------------------------------


# Mechanism orders are tuples of :class:`SpliceMechanismEffect`
# subclasses — class identity IS the mechanism, no parallel enum.

# Disrupted canonical donor (intronic +1/+2): exon skipping dominates,
# intron retention common, leaky normal splicing rare.
_MECHANISM_ORDER_SPLICE_DONOR = (
    ExonSkipping,
    IntronRetention,
    CrypticDonor,
    NormalSplicing,
)

# Disrupted canonical acceptor (intronic -1/-2): same rough
# distribution as donor.
_MECHANISM_ORDER_SPLICE_ACCEPTOR = (
    ExonSkipping,
    IntronRetention,
    CrypticAcceptor,
    NormalSplicing,
)

# Disrupted exonic splice site (last 3 of exon): often splicing still
# proceeds (the disruption is in the exon, where the spliceosome can
# tolerate more variation), so normal splicing is more competitive.
_MECHANISM_ORDER_EXONIC_SPLICE_SITE = (
    NormalSplicing,
    ExonSkipping,
    CrypticDonor,
    IntronRetention,
)

# Intronic splice site (positions +3 to +6 or -3): less critical
# region; normal splicing dominates. Two variants — donor-side (+3 to +6)
# and acceptor-side (-3) — which differ only in the cryptic direction.
_MECHANISM_ORDER_INTRONIC_SPLICE_SITE_DONOR = (
    NormalSplicing,
    ExonSkipping,
    IntronRetention,
    CrypticDonor,
)

_MECHANISM_ORDER_INTRONIC_SPLICE_SITE_ACCEPTOR = (
    NormalSplicing,
    ExonSkipping,
    IntronRetention,
    CrypticAcceptor,
)


def intronic_splice_side_is_acceptor(splice_effect):
    """True if an IntronicSpliceSite effect is on the acceptor side of
    its nearest exon.

    The classifier emits IntronicSpliceSite for positions +3–6 (donor
    side, after exon in transcript order) and -3 (acceptor side,
    before exon in transcript order). The side determines whether a
    cryptic donor or cryptic acceptor is the relevant alternative.

    Caller contract: ``splice_effect`` must be an ``IntronicSpliceSite``
    (guaranteed by :func:`mechanism_order_for`), so ``variant`` and
    ``nearest_exon`` are always present.
    """
    exon = splice_effect.nearest_exon
    variant = splice_effect.variant
    if exon.strand == "+":
        return variant.trimmed_base1_start < exon.start
    # Reverse strand: acceptor-side intronic variants sit past the
    # genomic end of the exon (which is the 5' end in transcript order).
    return variant.trimmed_base1_end > exon.end


def mechanism_order_for(splice_effect):
    """Pick the ordered list of splice mechanisms to enumerate for a
    given splice-disruption effect.

    When a variant disrupts a splice signal, several mechanisms are
    biologically possible (exon skipping, intron retention, cryptic-
    site use, or normal splicing). DNA alone can't say which one fires,
    so :func:`enumerate_splice_outcomes` emits *all* of them as
    candidates — and this function decides the order they appear in,
    which sets the ``most_likely`` candidate. The ordering is a hand-
    tuned, DNA-only heuristic (see the ``_MECHANISM_ORDER_*`` tuples
    above); it is not a probability.

    Returns the matching ``_MECHANISM_ORDER_*`` tuple, or ``None`` if
    ``splice_effect`` is not a splice disruption — in which case the
    caller leaves the effect unwrapped (no ``SpliceOutcomeSet``).

    All four splice-disruption classes are wrapped, including
    ``IntronicSpliceSite``. Dispatch is most-specific-first because
    ``SpliceDonor`` and ``SpliceAcceptor`` are *subclasses* of
    ``IntronicSpliceSite`` (see effect_classes.py): they must be tested
    before the ``IntronicSpliceSite`` branch, otherwise a canonical
    donor/acceptor hit (+1/+2 or -1/-2) would fall through to the
    weaker +3..+6 / -3 intronic-window ordering. ``IntronicSpliceSite``
    additionally needs a runtime side check to choose between a cryptic
    donor and a cryptic acceptor, which is why it can't share a fixed
    table entry with the others.

    Skipping an exon joins the donor before it to the acceptor after
    it, so a transcript's first and last exons can't be skipped. When
    the affected exon is one of them, ``ExonSkipping`` is left out and
    the next mechanism leads (#543).
    """
    if isinstance(splice_effect, SpliceDonor):          # intron +1/+2
        order = _MECHANISM_ORDER_SPLICE_DONOR
    elif isinstance(splice_effect, SpliceAcceptor):     # intron -1/-2
        order = _MECHANISM_ORDER_SPLICE_ACCEPTOR
    elif isinstance(splice_effect, ExonicSpliceSite):   # last 3 of exon
        order = _MECHANISM_ORDER_EXONIC_SPLICE_SITE
    elif isinstance(splice_effect, IntronicSpliceSite):  # intron +3..+6 / -3
        if intronic_splice_side_is_acceptor(splice_effect):
            order = _MECHANISM_ORDER_INTRONIC_SPLICE_SITE_ACCEPTOR
        else:
            order = _MECHANISM_ORDER_INTRONIC_SPLICE_SITE_DONOR
    else:
        return None
    if _affects_terminal_exon(splice_effect):
        order = tuple(m for m in order if m is not ExonSkipping)
    return order


def _affects_terminal_exon(splice_effect):
    """True if the splice signal belongs to its transcript's first or
    last exon, counting exons joined at a circular origin as one."""
    transcript = getattr(splice_effect, "transcript", None)
    exon = affected_exon(splice_effect)
    if transcript is None or exon is None:
        return False
    exons = transcript.exons
    return any(skipped.id in (exons[0].id, exons[-1].id)
               for skipped in _skipped_exons(transcript, exon))


def _skipped_exons(transcript, exon):
    """The exons skipped with ``exon``: both exons meeting at a circular
    origin, which are one exon on the circle, otherwise ``exon`` alone."""
    exons = list(transcript.exons)
    number = next((index for index, candidate in enumerate(exons, 1)
                   if candidate.id == exon.id), None)
    if number is None:
        return (exon,)
    return tuple(exons[partner - 1] for partner in
                 origin_partners(transcript, number))


# ---------------------------------------------------------------------
# Splice-mechanism Effect construction
# ---------------------------------------------------------------------


def wrap_splice_candidate(mechanism_effect):
    """Lift a :class:`SpliceMechanismEffect` into an
    unscored :class:`EffectCandidate` with ``source="varcode"``.

    Trivial wrapper kept as a named helper for symmetry with how
    other annotators construct ``EffectCandidate`` objects.
    """
    return EffectCandidate(
        effect=mechanism_effect,
        source="varcode",
    )


def build_normal_splicing_candidate(splice_effect):
    """Normal splicing: the splice signal was hit but the spliceosome
    handles it; the protein consequence (if any) is whatever the
    underlying coding change normally produces.

    For :class:`ExonicSpliceSite`, ``alternate_effect`` is the
    underlying coding change (a :class:`Substitution` etc.). For
    intronic splice classes there is no underlying coding change —
    the variant is in/near an intron — so ``coding_effect=None``.
    """
    transcript = getattr(splice_effect, "transcript", None)
    coding_effect = getattr(splice_effect, "alternate_effect", None)
    return wrap_splice_candidate(
        NormalSplicing(
            variant=splice_effect.variant,
            transcript=transcript,
            splice_signal=splice_effect,
            coding_effect=coding_effect),
    )


def aa_fields_from_protein_diff(
        variant, transcript, mutant_transcript, length_delta):
    """Run :func:`classify_from_protein_diff` against the resolved
    mutant transcript and return its ``aa_ref`` / ``aa_alt`` /
    offset fields as a dict, plus the classified protein effect. The
    splice-mechanism Effect builders copy the AA vocab onto the
    mechanism instance and retain the classified effect under
    ``protein_effect`` so severity queries can delegate to it while
    class identity remains the splice mechanism.
    """
    ref_protein = str(transcript.protein_sequence)
    mut_protein = mutant_transcript.mutant_protein_sequence
    classified = classify_from_protein_diff(
        variant=variant,
        transcript=transcript,
        ref_protein=ref_protein,
        mut_protein=mut_protein,
        length_delta=length_delta,
        mutant_transcript=mutant_transcript,
    )
    return {
        "aa_ref": getattr(classified, "aa_ref", None),
        "aa_alt": getattr(classified, "aa_alt", None),
        "aa_mutation_start_offset": getattr(
            classified, "aa_mutation_start_offset", None),
        "aa_mutation_end_offset": getattr(
            classified, "aa_mutation_end_offset", None),
        "mutant_protein_sequence": mut_protein,
        "protein_effect": classified,
    }


def build_exon_skipping_candidate(splice_effect):
    """Exon skipping: the affected exon is excluded from the mature
    transcript.

    Builds a :class:`MutantTranscript` by deleting the exon's cDNA
    range, translates, runs :func:`classify_from_protein_diff` to
    extract the protein-level AA vocab, and copies those fields onto
    a fresh :class:`ExonSkipping` instance. Returns the
    :class:`ExonSkipping` in its "predicted" state (protein fields
    ``None``) when the protein math couldn't resolve.
    """
    transcript = getattr(splice_effect, "transcript", None)
    exon = affected_exon(splice_effect)
    if transcript is None or exon is None:
        # We can't even compute frame without the exon.
        return wrap_splice_candidate(
            ExonSkipping(
                variant=splice_effect.variant,
                transcript=transcript,
                splice_signal=splice_effect,
                affected_exon=exon,
                in_frame=False),
        )

    skipped = _skipped_exons(transcript, exon)
    skipped_length = sum(exon_length(each) for each in skipped)
    in_frame = (skipped_length % 3 == 0)
    mt = build_exon_skip_mutant_transcript(
        splice_effect.variant, transcript, exon, skipped=skipped)

    if mt is None or mt.mutant_protein_sequence is None:
        # Couldn't compute protein (incomplete transcript, exon not
        # found, etc.) — return predicted-state ExonSkipping.
        return wrap_splice_candidate(
            ExonSkipping(
                variant=splice_effect.variant,
                transcript=transcript,
                splice_signal=splice_effect,
                affected_exon=exon,
                in_frame=in_frame,
                mutant_transcript=mt),
        )

    aa_fields = aa_fields_from_protein_diff(
        splice_effect.variant, transcript, mt, length_delta=-skipped_length)
    return wrap_splice_candidate(
        ExonSkipping(
            variant=splice_effect.variant,
            transcript=transcript,
            splice_signal=splice_effect,
            affected_exon=exon,
            in_frame=in_frame,
            mutant_transcript=mt,
            **aa_fields),
    )




def build_exon_skip_mutant_transcript(variant, transcript, exon, skipped=None):
    """Build a :class:`MutantTranscript` representing an exon skip.

    Deletes the exon's cDNA range from the transcript sequence and
    translates the result. Adjusts the CDS start offset when the
    skipped exon precedes the coding region (5'UTR exon). When the
    skipped exon contains the start codon, returns a MutantTranscript
    with an empty protein so the classifier emits StartLoss.
    ``skipped`` names consecutive exons skipped together, by default
    ``(exon,)``.
    """
    if not has_complete_coding_sequence(transcript):
        return None
    skipped = tuple(skipped or (exon,))
    try:
        exon_start_in_tx = min(
            exon_start_offset_in_transcript(transcript, each)
            for each in skipped)
    except (StopIteration, ValueError):
        return None

    length = sum(exon_length(each) for each in skipped)
    full_sequence = str(transcript.sequence)

    post_skip_cdna = (
        full_sequence[:exon_start_in_tx]
        + full_sequence[exon_start_in_tx + length:]
    )

    cds_start = min(transcript.start_codon_spliced_offsets)

    # Adjust CDS start for the deletion.
    if exon_start_in_tx + length <= cds_start:
        # Exon is entirely in the 5' UTR → shift CDS start left.
        new_cds_start = cds_start - length
    elif exon_start_in_tx >= cds_start:
        # Exon is entirely after the CDS start → no shift.
        new_cds_start = cds_start
    else:
        # Exon overlaps the CDS start → start codon lost.
        new_cds_start = None

    edit = TranscriptEdit(
        cdna_start=exon_start_in_tx,
        cdna_end=exon_start_in_tx + length,
        alt_bases="",
        source_variant=variant,
    )

    mut_protein = ""
    if new_cds_start is not None and new_cds_start < len(post_skip_cdna):
        codon_table = codon_table_for_transcript(transcript)
        coding = post_skip_cdna[new_cds_start:]
        coding += poly_a_completion(transcript, coding)
        truncated = coding[:(len(coding) // 3) * 3]
        selenocysteine = {pos - new_cds_start
                          for pos in edited_selenocysteine(transcript, (edit,))}
        try:
            mut_protein = translate_sequence(
                truncated, codon_table=codon_table, to_stop=True,
                selenocysteine=selenocysteine)
        except ValueError:
            mut_protein = ""

    return MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=post_skip_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )


def affected_exon(splice_effect):
    """Best-effort extraction of the affected exon from a splice effect.

    ExonicSpliceSite carries it directly. SpliceDonor/SpliceAcceptor
    carry the nearest exon. IntronicSpliceSite likewise.
    """
    if hasattr(splice_effect, "exon"):
        return splice_effect.exon
    if hasattr(splice_effect, "nearest_exon"):
        return splice_effect.nearest_exon
    return None



def exon_start_offset_in_transcript(transcript, exon):
    """Offset (in transcript coordinates) of the first base of the
    given exon.
    """
    offset = 0
    for ex in transcript.exons:
        if ex.exon_id == exon.exon_id:
            return offset
        offset += exon_length(ex)
    raise ValueError("Exon %s not found in transcript %s" % (
        exon.exon_id, transcript.transcript_id))




__all__ = [
    "aa_fields_from_protein_diff",
    "affected_exon",
    "build_exon_skip_mutant_transcript",
    "build_exon_skipping_candidate",
    "build_normal_splicing_candidate",
    "exon_start_offset_in_transcript",
    "intronic_splice_side_is_acceptor",
    "mechanism_order_for",
    "wrap_splice_candidate",
]
