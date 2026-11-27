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

"""Public barrel re-exporting all Varcode effect classes.

Effect classes are organized by domain across:

* :mod:`varcode.effects.base` — base hierarchy and non-coding transcript effects
* :mod:`varcode.effects.splice_effects` — splice-signal and splice-mechanism effects
* :mod:`varcode.effects.coding_effects` — coding-sequence and amino-acid effects
* :mod:`varcode.effects.structural_effects` — structural, haplotype, and phase effects
"""

from .base import (
    Exonic,
    Failure,
    FivePrimeUTR,
    GermlineAlleleOverlap,
    HypothesisLimit,
    IncompleteTranscript,
    Intergenic,
    Intragenic,
    Intronic,
    MultiOutcomeEffect,
    MutationEffect,
    NoncodingTranscript,
    ThreePrimeUTR,
    TranscriptMutationEffect,
    Unresolved,
)
from .coding_effects import (
    AlternateStartCodon,
    CodingMutation,
    ComplexSubstitution,
    Deletion,
    FrameShift,
    FrameShiftTruncation,
    Insertion,
    KnownAminoAcidChange,
    NonsilentCodingMutation,
    PrematureStop,
    Silent,
    StartGain,
    StartLoss,
    StopLoss,
    Substitution,
)
from .splice_effects import (
    CrypticAcceptor,
    CrypticDonor,
    CrypticSpliceSiteEffect,
    ExonicSpliceSite,
    ExonLoss,
    ExonSkipping,
    IntronicSpliceSite,
    IntronRetention,
    NormalSplicing,
    SpliceAcceptor,
    SpliceDonor,
    SpliceMechanismEffect,
    SpliceSite,
)
from .structural_effects import (
    CrypticExonCandidate,
    GeneFusion,
    HaplotypeEffect,
    Inversion,
    LargeDeletion,
    LargeDuplication,
    PhaseCandidateSet,
    StructuralVariantEffect,
    TranslocationToIntergenic,
    TranslocationToUnannotatedContig,
)

__all__ = [
    # base & non-coding
    "MutationEffect",
    "MultiOutcomeEffect",
    "Intergenic",
    "Intragenic",
    "TranscriptMutationEffect",
    "GermlineAlleleOverlap",
    "Failure",
    "Unresolved",
    "HypothesisLimit",
    "NoncodingTranscript",
    "IncompleteTranscript",
    "FivePrimeUTR",
    "ThreePrimeUTR",
    "Intronic",
    "Exonic",
    # splice
    "SpliceSite",
    "IntronicSpliceSite",
    "SpliceDonor",
    "SpliceAcceptor",
    "SpliceMechanismEffect",
    "NormalSplicing",
    "ExonSkipping",
    "IntronRetention",
    "CrypticSpliceSiteEffect",
    "CrypticDonor",
    "CrypticAcceptor",
    "ExonLoss",
    "ExonicSpliceSite",
    # coding
    "CodingMutation",
    "Silent",
    "AlternateStartCodon",
    "NonsilentCodingMutation",
    "StartGain",
    "StartLoss",
    "KnownAminoAcidChange",
    "Substitution",
    "ComplexSubstitution",
    "Insertion",
    "Deletion",
    "PrematureStop",
    "StopLoss",
    "FrameShift",
    "FrameShiftTruncation",
    # structural, haplotype, phase
    "StructuralVariantEffect",
    "LargeDeletion",
    "LargeDuplication",
    "Inversion",
    "GeneFusion",
    "TranslocationToIntergenic",
    "TranslocationToUnannotatedContig",
    "CrypticExonCandidate",
    "HaplotypeEffect",
    "PhaseCandidateSet",
]
