# PYTHON_ARGCOMPLETE_OK
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

import sys

import argcomplete

from .variant_args import (
    add_variant_args,
    check_output_directory,
    configure_cli_logging,
    exit_on_user_error,
    make_variants_parser,
    variant_collection_from_args,
    write_cli_dataframe,
)


def add_effects_args(parser):
    """Attach ``varcode effects`` options onto ``parser``."""
    parser.add_argument(
        "--output-csv", "--output", "-o", dest="output_csv",
        help="CSV output: variant coordinates, gene/transcript IDs, effect type and description, plus SV fields",
    )
    parser.add_argument(
        "--extra-orfs", action="append", metavar="EVIDENCE",
        help="Add ORFdata consequences: all, ribosome, ms, proteomics or immunopeptidomics. "
             "Repeat to include ORFs with any selected evidence. Uses the ORF catalogues "
             "installed for the genome; `orfdata setup` installs more.",
    )
    parser.add_argument(
        "--one-per-variant",
        default=False,
        action="store_true",
        help=(
            "Only return highest priority effect overlapping a variant, "
            "otherwise all overlapping transcripts are returned."
        ),
    )
    parser.add_argument(
        "--only-coding",
        default=False,
        action="store_true",
        help="Filter silent and non-coding effects",
    )
    return parser


def attach_effects_parser(parser):
    """Attach both shared variant options and ``effects`` options onto a subparser."""
    add_variant_args(parser)
    return add_effects_args(parser)


arg_parser = make_variants_parser(
    prog="varcode effects",
    description="Annotate variants with overlapping gene names and predicted coding effects",
    epilog="Other subcommands: varcode genes (see `varcode --help`)",
    add_version=True,
)
add_effects_args(arg_parser)


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
    if args_list is None:
        argcomplete.autocomplete(arg_parser)
        args_list = sys.argv[1:]

    args = arg_parser.parse_args(args_list)
    configure_cli_logging()

    from ..effects import Failure
    from ..errors import ReferenceMismatchError

    with exit_on_user_error(arg_parser, errors=(ValueError, OSError, ImportError)):
        if args.output_csv:
            check_output_directory(args.output_csv)
        variants = variant_collection_from_args(args)
        from ..extra_orfs import prepare_extra_orfs
        selection = args.extra_orfs
        if selection and "all" in selection:
            if len(selection) != 1:
                raise ValueError("--extra-orfs all cannot be combined with evidence filters")
            selection = "all"
        extra_orfs = prepare_extra_orfs(selection)
        if extra_orfs is not None:
            for genome in {id(v.genome): v.genome for v in variants}.values():
                extra_orfs.for_genome(genome)
    try:
        effects = variants.effects(
            raise_on_error=False if args.skip_errors else None, extra_orfs=extra_orfs)
    except ReferenceMismatchError as error:
        arg_parser.exit(1, "%s: error: %s Use --skip-errors to retain failed annotations.\n" % (
            arg_parser.prog, error))
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

    effects_dataframe = effects.to_dataframe(include_archive=True)
    if args.skip_errors or failures:
        effects_dataframe["annotation_error"] = [getattr(e, "error", None) for e in effects]
    if args.output_csv:
        with exit_on_user_error(arg_parser, errors=(OSError,)):
            effects_dataframe.to_csv(args.output_csv, index=False)
    else:
        write_cli_dataframe(effects_dataframe, drop_display_columns=("effect_json",))
