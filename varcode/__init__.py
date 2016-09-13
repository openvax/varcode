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
from .variant_collection import VariantCollection
from .maf import load_maf, load_maf_dataframe
from .vcf import load_vcf, load_vcf_fast
from .effects import (
    effect_priority,
    top_priority_effect,
    EffectCollection,
    MutationEffect,
    NonsilentCodingMutation,
)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = [
    # basic classes
    "Variant",
    "EffectCollection",
    "VariantCollection",
    # effects
    "effect_priority",
    "top_priority_effect",
    "MutationEffect",
    "NonsilentCodingMutation",
    # file loading
    "load_maf",
    "load_maf_dataframe",
    "load_vcf",
    "load_vcf_fast",
]
