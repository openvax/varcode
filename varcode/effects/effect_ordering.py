# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import

from .effect_classes import (
    Failure,
    IncompleteTranscript,
    Intergenic,
    Intragenic,
    NoncodingTranscript,
    Intronic,
    ThreePrimeUTR,
    FivePrimeUTR,
    Silent,
    Substitution,
    Insertion,
    Deletion,
    ComplexSubstitution,
    AlternateStartCodon,
    IntronicSpliceSite,
    ExonicSpliceSite,
    StopLoss,
    SpliceDonor,
    SpliceAcceptor,
    PrematureStop,
    FrameShiftTruncation,
    StartLoss,
    FrameShift,
    ExonLoss,
)
from ..common import apply_groupby

transcript_effect_priority_list = [
    Failure,
    IncompleteTranscript,
    Intergenic,
    Intragenic,
    NoncodingTranscript,
    Intronic,
    ThreePrimeUTR,
    # mutations to the upstream 5' UTR may change the ORF (reading frame),
    # so give 5' UTR mutations higher prioriry
    FivePrimeUTR,
    Silent,
    # intronic variants near the splice boundaries but which aren't
    # the two nucleotides closest to the exon
    IntronicSpliceSite,
    # exonic variants near a splice boundary
    ExonicSpliceSite,
    # mutation in the two nucleotides immediately following an exon/intron
    # boundary
    SpliceDonor,
    # mutation in the two nucleotides immediately preceding an intron/exon
    # boundary
    SpliceAcceptor,
    Substitution,
    Insertion,
    Deletion,
    ComplexSubstitution,
    # mutation which changes the start codon from e.g. ATG > TTG that can
    # be interpreted as silent but also has some chance of causing an
    # alternative ORF
    AlternateStartCodon,
    # modification or deletion of stop codon
    StopLoss,
    # creation of a new stop codon within the coding sequence
    PrematureStop,
    # frame-shift which creates immediate stop codon, same as PrematureStop
    FrameShiftTruncation,
    # modification or deletion of a start codon
    StartLoss,
    # out-of-frame insertion or deletion
    FrameShift,
    ExonLoss,
]

transcript_effect_priority_dict = {
    transcript_effect_class: priority
    for (priority, transcript_effect_class)
    in enumerate(transcript_effect_priority_list)
}


def effect_priority(effect):
    """
    Returns the integer priority for a given transcript effect.
    """
    return transcript_effect_priority_dict.get(effect.__class__, -1)


def apply_to_field_if_exists(effect, field_name, fn, default):
    """
    Apply function to specified field of effect if it is not None,
    otherwise return default.
    """
    value = getattr(effect, field_name, None)
    if value is None:
        return default
    else:
        return fn(value)


def apply_to_transcript_if_exists(effect, fn, default):
    """
    Apply function to transcript associated with effect,
    if it exists, otherwise return default.
    """
    return apply_to_field_if_exists(
        effect=effect,
        field_name="transcript",
        fn=fn,
        default=default)


def apply_to_gene_if_exists(effect, fn, default):
    return apply_to_field_if_exists(
        effect=effect,
        field_name="gene",
        fn=fn,
        default=default)


def number_exons_in_associated_transcript(effect):
    """
    Number of exons on transcript associated with effect,
    if there is one (otherwise return 0).
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: len(t.exons),
        default=0)


def cds_length_of_associated_transcript(effect):
    """
    Length of coding sequence of transcript associated with effect,
    if there is one (otherwise return 0).
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: len(t.coding_sequence) if (t.complete and t.coding_sequence) else 0,
        default=0)


def length_of_associated_transcript(effect):
    """
    Length of spliced mRNA sequence of transcript associated with effect,
    if there is one (otherwise return 0).
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: len(t.sequence),
        default=0)


def name_of_associated_transcript(effect):
    """
    Name of transcript associated with effect,
    if there is one (otherwise return "").
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: t.name,
        default="")


def gene_id_of_associated_transcript(effect):
    """
    Ensembl gene ID of transcript associated with effect, returns
    None if effect does not have transcript.
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: t.gene_id,
        default=None)


def effect_has_complete_transcript(effect):
    """
    Parameters
    ----------
    effect : subclass of MutationEffect

    Returns True if effect has transcript and that transcript has complete CDS
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: t.complete,
        default=False)


def effect_associated_with_protein_coding_gene(effect):
    """
    Parameters
    ----------
    effect : subclass of MutationEffect

    Returns True if effect is associated with a gene and that gene
    has a protein_coding biotype.
    """
    return apply_to_gene_if_exists(
        effect=effect,
        fn=lambda g: g.biotype == "protein_coding",
        default=False)


def effect_associated_with_protein_coding_transcript(effect):
    """
    Parameters
    ----------
    effect : subclass of MutationEffect

    Returns True if effect is associated with a transcript and that transcript
    has a protein_coding biotype.
    """
    return apply_to_transcript_if_exists(
        effect=effect,
        fn=lambda t: t.biotype == "protein_coding",
        default=False)


def transcript_name_ends_with_01(effect):
    """
    Often the canonical transcript seems to be named something like
    "TP53-001" or "TP53-201", so as a last-ditch tie breaker we're
    using whether the name ends in "01".

    Parameters
    ----------
    effect : subclass of MutationEffect

    Returns bool
    """
    return name_of_associated_transcript(effect).endswith("01")


def parse_transcript_number(effect):
    """
    Try to parse the number at the end of a transcript name associated with
    an effect.

    e.g. TP53-001 returns the integer 1.

    Parameters
    ----------
    effect : subclass of MutationEffect

    Returns int
    """
    name = name_of_associated_transcript(effect)
    if "-" not in name:
        return 0
    parts = name.split("-")
    last_part = parts[-1]
    if last_part.isdigit():
        return int(last_part)
    else:
        return 0


def multi_gene_effect_sort_key(effect):
    """
    This function acts as a sort key for choosing the highest priority
    effect across multiple genes (so does not assume that effects might
    involve the same start/stop codons).

    Returns tuple with the following elements:
        1) Integer priority of the effect type.
        2) Does the associated gene have a "protein_coding" biotype?
            False if no gene is associated with effect.
        3) Does the associated transcript have a "protein_coding" biotype?
            False if no transcript is associated with effect.
        4) Is the associated transcript complete?
            False if no transcript is associated with effect.
        5) CDS length
             This value will be 0 if the effect has no associated transcript
             or if the transcript is noncoding or incomplete
        6) Total length of the transcript
             This value will be 0 intra/intergenic variants effects without
             an associated transcript.
        7) Number of exons
             This value will be 0 intra/intergenic variants effects without
             an associated transcript.
        8) If everything is the same up this point then let's use the very
           sloppy heuristic of preferring transcripts like "TP53-201" over
           "TP53-206", so anything ending with "01" is considered better.
        9) Lastly, if we end up with two transcripts like "TP53-202" and
           "TP53-203", prefer the one with the lowest number in the name.
    """

    return tuple([
        effect_priority(effect),
        effect_associated_with_protein_coding_gene(effect),
        effect_associated_with_protein_coding_transcript(effect),
        effect_has_complete_transcript(effect),
        cds_length_of_associated_transcript(effect),
        length_of_associated_transcript(effect),
        number_exons_in_associated_transcript(effect),
        transcript_name_ends_with_01(effect),
        -parse_transcript_number(effect)
    ])


def select_between_exonic_splice_site_and_alternate_effect(effect):
    """
    If the given effect is an ExonicSpliceSite then it might contain
    an alternate effect of higher priority. In that case, return the
    alternate effect. Otherwise, this acts as an identity function.
    """
    if effect.__class__ is not ExonicSpliceSite:
        return effect
    if effect.alternate_effect is None:
        return effect
    splice_priority = effect_priority(effect)
    alternate_priority = effect_priority(effect.alternate_effect)
    if splice_priority > alternate_priority:
        return effect
    else:
        return effect.alternate_effect


def keep_max_priority_effects(effects):
    """
    Given a list of effects, only keep the ones with the maximum priority
    effect type.

    Parameters
    ----------
    effects : list of MutationEffect subclasses

    Returns list of same length or shorter
    """
    priority_values = tuple(map(effect_priority, effects))
    max_priority = max(priority_values)
    return [e for (e, p) in zip(effects, priority_values) if p == max_priority]


def keep_effects_on_protein_coding_transcripts(effects):
    """
    Given a list of effects, only keep the ones associated with genes
    of a protein coding biotype.

    Parameters
    ----------
    effects : list of MutationEffect subclasses

    Returns list of same length or shorter
    """
    return [
        e for e in effects if effect_associated_with_protein_coding_transcript(e)
    ]


def keep_effects_on_protein_coding_genes(effects):
    """
    Given a list of effects, only keep the ones associated with genes
    of a protein coding biotype.

    Parameters
    ----------
    effects : list of MutationEffect subclasses

    Returns list of same length or shorter
    """
    return [e for e in effects if effect_associated_with_protein_coding_gene(e)]


def keep_effects_on_complete_transcripts(effects):
    """
    Given a list of effects, only keep the ones associated with complete
    transcripts (which have annotated start and stop codons).

    Parameters
    ----------
    effects : list of MutationEffect subclasses

    Returns list of same length or shorter
    """
    return [e for e in effects if effect_has_complete_transcript(e)]


def filter_pipeline(effects, filters):
    """
    Apply each filter to the effect list sequentially. If any filter
    returns zero values then ignore it. As soon as only one effect is left,
    return it.

    Parameters
    ----------
    effects : list of MutationEffect subclass instances

    filters : list of functions
        Each function takes a list of effects and returns a list of effects

    Returns list of effects
    """
    for filter_fn in filters:
        if len(effects) == 1:
            return effects
        filtered_effects = filter_fn(effects)
        if len(filtered_effects) == 0:
            return effects
        effects = filtered_effects
    return effects


def tie_breaking_sort_key_for_single_gene_effects(effect):
    """
    This function assumes that effects have already been filtered
    by effect priority and we now have a collision where multiple effects
    arising from different transcripts. Use this sort key to pick transcripts
    for the following features:

        1) Highest number of exons
        2) Transcript name ending with "01"
        3) Lowest number in transcript name (e.g. TP53-002 over TP53-003)

    Returns tuple.
    """

    return tuple([
        number_exons_in_associated_transcript(effect),
        transcript_name_ends_with_01(effect),
        -parse_transcript_number(effect),
    ])


def top_priority_effect_for_single_gene(effects):
    """
    For effects which are from the same gene, check to see if there
    is a canonical transcript with both the maximum length CDS
    and maximum length full transcript sequence.

    If not, then use number of exons and transcript name as tie-breaking
    features.

    Parameters
    ----------
    effects : list of MutationEffect subclass instances

    Returns single effect object
    """

    # first filter effects to keep those on
    # 1) maximum priority effects
    # 2) protein coding genes
    # 3) protein coding transcripts
    # 4) complete transcripts
    #
    # If any of these filters drop all the effects then we move on to the next
    # filtering step.
    effects = filter_pipeline(
        effects=effects,
        filters=[
            keep_max_priority_effects,
            keep_effects_on_protein_coding_genes,
            keep_effects_on_protein_coding_transcripts,
            keep_effects_on_complete_transcripts,
        ],
    )

    if len(effects) == 1:
        return effects[0]

    # compare CDS length and transcript lengths of remaining effects
    # if one effect has the maximum of both categories then return it
    cds_lengths = [cds_length_of_associated_transcript(e) for e in effects]
    max_cds_length = max(cds_lengths)

    # get set of indices of all effects with maximum CDS length
    max_cds_length_indices = {
        i
        for (i, l) in enumerate(cds_lengths)
        if l == max_cds_length
    }

    seq_lengths = [length_of_associated_transcript(e) for e in effects]
    max_seq_length = max(seq_lengths)

    # get set of indices for all effects whose associated transcript
    # has maximum sequence length
    max_seq_length_indices = {
        i
        for (i, l) in enumerate(seq_lengths)
        if l == max_seq_length
    }

    # which effects have transcripts with both the longest CDS and
    # longest full transcript sequence?
    intersection_of_indices = \
        max_cds_length_indices.intersection(max_seq_length_indices)

    n_candidates = len(intersection_of_indices)
    if n_candidates == 1:
        best_index = intersection_of_indices.pop()
        return effects[best_index]

    elif n_candidates == 0:
        # if set of max CDS effects and max sequence length effects is disjoint
        # then let's try to do the tie-breaking sort over their union
        union_of_indices = max_cds_length_indices.union(max_seq_length_indices)
        candidate_effects = [effects[i] for i in union_of_indices]

    else:
        # if multiple effects have transcripts with the max CDS length and
        # the max full sequence length then run the tie-breaking sort
        # over all these candidates
        candidate_effects = [effects[i] for i in intersection_of_indices]

    # break ties by number of exons, whether name of transcript ends if "01",
    # and all else being equal, prefer transcript names that end with lower
    # numbers
    return max(
        candidate_effects,
        key=tie_breaking_sort_key_for_single_gene_effects)


def top_priority_effect(effects):
    """
    Given a collection of variant transcript effects,
    return the top priority object. ExonicSpliceSite variants require special
    treatment since they actually represent two effects -- the splicing modification
    and whatever else would happen to the exonic sequence if nothing else gets
    changed. In cases where multiple transcripts give rise to multiple
    effects, use a variety of filtering and sorting heuristics to pick
    the canonical transcript.
    """
    if len(effects) == 0:
        raise ValueError("List of effects cannot be empty")

    effects = map(
        select_between_exonic_splice_site_and_alternate_effect,
        effects)

    effects_grouped_by_gene = apply_groupby(
        effects, fn=gene_id_of_associated_transcript, skip_none=False)

    if None in effects_grouped_by_gene:
        effects_without_genes = effects_grouped_by_gene.pop(None)
    else:
        effects_without_genes = []

    # if we had any effects associated with genes then choose one of those
    if len(effects_grouped_by_gene) > 0:
        effects_with_genes = [
            top_priority_effect_for_single_gene(gene_effects)
            for gene_effects in effects_grouped_by_gene.values()
        ]
        return max(effects_with_genes, key=multi_gene_effect_sort_key)
    else:
        # if all effects were without genes then choose the best among those
        assert len(effects_without_genes) > 0
        return max(effects_without_genes, key=multi_gene_effect_sort_key)

