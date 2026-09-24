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

import logging.config
from importlib import resources
import sys

from .version_info import print_version_info
from .variant_args import make_variants_parser, variant_collection_from_args
from ..effects import Failure


# as_file, rather than str(), so that logging.conf is also readable when varcode
# is zip-imported (zipapp, pex, shiv, zipped egg), where files() hands back a
# zipfile.Path that open() can't use.
with resources.as_file(resources.files(__package__) / "logging.conf") as _logging_conf:
    logging.config.fileConfig(_logging_conf)
logger = logging.getLogger(__name__)

arg_parser = make_variants_parser(
    description="Annotate variants with overlapping gene names and predicted coding effects")

arg_parser.add_argument(
    "--output-csv", help="CSV output: variant coordinates, gene/transcript IDs, effect type and description, plus SV fields")

arg_parser.add_argument(
    "--one-per-variant",
    default=False,
    action="store_true",
    help=(
        "Only return highest priority effect overlapping a variant, "
        "otherwise all overlapping transcripts are returned."))

arg_parser.add_argument(
    "--only-coding",
    default=False,
    action="store_true",
    help="Filter silent and non-coding effects")

def main(args_list=None):
    """
    Script which loads variants and annotates them with overlapping genes
    and predicted coding effects.

    Example usage:
        varcode
            --vcf mutect.vcf \
            --vcf strelka.vcf \
            --maf tcga_brca.maf \
            --variant chr1 498584 C G \
            --json-variants more_variants.json
    """
    print_version_info()
    if args_list is None:
        args_list = sys.argv[1:]

    args = arg_parser.parse_args(args_list)
    variants = variant_collection_from_args(args)
    effects = variants.effects(raise_on_error=not args.skip_errors)
    # Annotation failures are diagnostics, not biological effects to rank or
    # discard with coding filters. Keep every failed transcript visible.
    failures = [effect for effect in effects if isinstance(effect, Failure)]
    effects = effects.filter(lambda effect: not isinstance(effect, Failure))
    if args.only_coding:
        effects = effects.drop_silent_and_noncoding()
    if args.one_per_variant:
        variant_to_effect_dict = effects.top_priority_effect_per_variant()
        effects = effects.clone_with_new_elements(list(variant_to_effect_dict.values()))

    effects = effects.clone_with_new_elements(list(effects) + failures)

    effects_dataframe = effects.to_dataframe()
    if args.skip_errors:
        effects_dataframe["annotation_error"] = [getattr(e, "error", None) for e in effects]
    logger.info('\n%s', effects)
    if args.output_csv:
        effects_dataframe.to_csv(args.output_csv, index=False)
