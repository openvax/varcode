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
    # changing the stop codon might have more significance than
    # silent variants internal to the coding sequence
    # REMOVED: AlternateStopCodon,
    Substitution,
    Insertion,
    Deletion,
    ComplexSubstitution,
    # mutation which changes the start codon from e.g. ATG > TTG that can
    # be interpreted as silent but also has some chance of causing an
    # alternative ORF
    AlternateStartCodon,
    # intronic variants near the splice boundaries but which aren't
    # the two nucleotides closest to the exon
    IntronicSpliceSite,
    # exonic variants near a splice boundary
    ExonicSpliceSite,
    # modification or deletion of stop codon
    StopLoss,
    # mutation in the two nucleotides immediately following an exon/intron
    # boundary
    SpliceDonor,
    # mutation in the two nucleotides immediately preceding an intron/exon
    # boundary
    SpliceAcceptor,
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
    Returns the integer priority for a given transcript effect
    """
    # since intergenic variants may have a None value for their
    # highest_priority effect it simplifies other code to handle None
    # here
    if effect is None:
        return -1
    return transcript_effect_priority_dict[effect.__class__]

def effect_sort_key(effect):
    """Returns key tuple with the following fields that should be sorted
    lexicographically:
        - what is the CDS length?
            This value will be 0 if the effect has no associated transcript
            or if the transcript is noncoding or incomplete
        - what is the total length of the transcript?
            This value will be 0 intra/intergenic variants effects without
            an associated transcript.
    """
    # by default, we'll return (0, 0) for effects without an associated
    # transcript
    transcript_length = 0
    cds_length = 0
    if hasattr(effect, 'transcript') and effect.transcript is not None:
        transcript_length = len(effect.transcript)
        if effect.transcript.complete:
            cds_length = len(effect.transcript.coding_sequence)
    effect_class_priority = effect_priority(effect)
    return (effect_class_priority, cds_length, transcript_length)

def top_priority_effect(effects):
    """
    Given a collection of variant transcript effects,
    return the top priority object. In case of multiple transcript
    effects with the same priority, return the one affecting the longest
    transcript.
    """
    if len(effects) == 0:
        raise ValueError("List of effects cannot be empty")
    return max(effects, key=effect_sort_key)
