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


def add_genes_args(parser):
    """Attach ``varcode genes`` options onto ``parser``."""
    parser.add_argument(
        "--output-csv", "--output", "-o", dest="output_csv",
        help="CSV output: contig, start, ref, alt, gene_name, gene_id, reference columns, plus SV fields",
    )
    return parser


def attach_genes_parser(parser):
    """Attach both shared variant options and ``genes`` options onto a subparser."""
    add_variant_args(parser)
    return add_genes_args(parser)


arg_parser = make_variants_parser(
    prog="varcode genes",
    description="Annotate variants with overlapping gene names",
)
add_genes_args(arg_parser)


def _normalize_genes_dataframe(variants, raise_on_error):
    """Return a CLI genes DataFrame sharing ``contig/start/ref/alt`` normalization with ``effects``."""
    df = variants.to_dataframe(raise_on_error=raise_on_error)
    df = df.rename(columns={"chr": "contig"})
    if len(variants):
        df["start"] = [v.start for v in variants]
        df["ref"] = [v.ref for v in variants]
        df["alt"] = [v.alt for v in variants]
    return df


def main(args_list=None):
    """
    Script which loads variants and annotates them with overlapping genes.

    Example usage:
        varcode genes
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

    with exit_on_user_error(arg_parser):
        if args.output_csv:
            check_output_directory(args.output_csv)
        variants = variant_collection_from_args(args)
    variants_dataframe = _normalize_genes_dataframe(
        variants, raise_on_error=False if args.skip_errors else None
    )
    if args.output_csv:
        with exit_on_user_error(arg_parser, errors=(OSError,)):
            variants_dataframe.to_csv(args.output_csv, index=False)
    else:
        write_cli_dataframe(variants_dataframe)
