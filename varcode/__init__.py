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

from .annotators import (
    EffectAnnotator,
    FastEffectAnnotator,
    UnsupportedVariantError,
    get_annotator,
    get_default_annotator,
    register_annotator,
    resolve_annotator,
    set_default_annotator,
    use_annotator,
)
from .errors import ReferenceMismatchError, SampleNotFoundError
from .genotype import Genotype, Zygosity
from .mutant_transcript import (
    MutantTranscript,
    TranscriptEdit,
    apply_variant_to_transcript,
)
from .splice_outcomes import (
    SpliceCandidate,
    SpliceOutcome,
    SpliceOutcomeSet,
)
from .variant import Variant
from .variant_collection import VariantCollection
from .maf import load_maf, load_maf_dataframe
from .effects import (
    effect_priority,
    top_priority_effect,
    EffectCollection,
    MultiOutcomeEffect,
    MutationEffect,
    NonsilentCodingMutation,
)
from .version import __version__


def __getattr__(name):
    # Lazy import for VCF loaders. Pulling in .vcf at package-import time
    # drags PyVCF3 into sys.modules, which in turn triggers rpy2/R
    # initialization if rpy2 happens to be installed in the environment
    # — a hundreds-of-ms surprise that every `import varcode` eats even
    # when no VCF is ever loaded. See #302.
    if name in ("load_vcf", "load_vcf_fast"):
        from . import vcf as _vcf
        value = getattr(_vcf, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "__version__",

    # basic classes
    "Variant",
    "EffectCollection",
    "VariantCollection",

    # genotype / zygosity
    "Genotype",
    "Zygosity",

    # splice outcome possibility set (openvax/varcode#262)
    "SpliceCandidate",
    "SpliceOutcome",
    "SpliceOutcomeSet",

    # MutantTranscript data model + pluggable annotators (openvax/varcode#271)
    "MutantTranscript",
    "TranscriptEdit",
    "apply_variant_to_transcript",
    "EffectAnnotator",
    "FastEffectAnnotator",
    "UnsupportedVariantError",
    "get_annotator",
    "get_default_annotator",
    "register_annotator",
    "resolve_annotator",
    "set_default_annotator",
    "use_annotator",

    # effects
    "effect_priority",
    "top_priority_effect",
    "MultiOutcomeEffect",
    "MutationEffect",
    "NonsilentCodingMutation",

    # exceptions
    "ReferenceMismatchError",
    "SampleNotFoundError",

    # file loading
    "load_maf",
    "load_maf_dataframe",
    "load_vcf",
    "load_vcf_fast",
]
