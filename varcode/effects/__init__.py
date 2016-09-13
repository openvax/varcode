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

from .effect_collection import EffectCollection
from .effect_ordering import (
    effect_priority,
    effect_sort_key,
    top_priority_effect,
)
from .effect_prediction import (
    predict_variant_effects,
    predict_variant_effect_on_transcript,
    predict_variant_effect_on_transcript_or_failure,
)
from .effect_classes import (
    MutationEffect,
    TranscriptMutationEffect,
    NonsilentCodingMutation,
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

__all__ = [
    "EffectCollection",
    # effect ordering
    "effect_priority",
    "effect_sort_key",
    "top_priority_effect",

    # prediction functions
    "predict_variant_effects",
    "predict_variant_effect_on_transcript",
    "predict_variant_effect_on_transcript_or_failure",

    # effect classes
    "MutationEffect",
    "TranscriptMutationEffect",
    "Failure",
    "IncompleteTranscript",
    "Intergenic",
    "Intragenic",
    "IncompleteTranscript",
    "NoncodingTranscript",
    "ThreePrimeUTR",
    "FivePrimeUTR",
    "Intronic",
    "Silent",
    "NonsilentCodingMutation",
    "Substitution",
    "Insertion",
    "Deletion",
    "ComplexSubstitution",
    "AlternateStartCodon",
    "IntronicSpliceSite",
    "ExonicSpliceSite",
    "StopLoss",
    "SpliceDonor",
    "SpliceAcceptor",
    "PrematureStop",
    "FrameShiftTruncation",
    "StartLoss",
    "FrameShift",
    "ExonLoss",
]
