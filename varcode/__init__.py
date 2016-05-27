# Copyright (c) 2015. Mount Sinai School of Medicine
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

from .variant import Variant
from .collection import Collection
from .effect_collection import EffectCollection
from .variant_collection import VariantCollection
from .maf import load_maf, load_maf_dataframe
from .vcf import load_vcf, load_vcf_fast
from .effect_ordering import effect_priority, top_priority_effect
from .effects import (
    MutationEffect,
    TranscriptMutationEffect,
    Failure,
    Intergenic,
    Intragenic,
    IncompleteTranscript,
    NoncodingTranscript,
    Intronic,
    ThreePrimeUTR,
    FivePrimeUTR,
    Silent,
    NonsilentCodingMutation,
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

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


__all__ = [
    # basic classes
    "Variant",
    "Collection",
    "EffectCollection",
    "VariantCollection",
    "effect_priority",
    "top_priority_effect",
    # file loading
    "load_maf",
    "load_maf_dataframe",
    "load_vcf",
    "load_vcf_fast",
    # effects
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
