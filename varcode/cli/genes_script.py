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

import logging
import sys

from .version_info import print_version_info
from .variant_args import (
    _configure_logging,
    check_output_directory,
    exit_on_user_error,
    make_variants_parser,
    variant_collection_from_args,
)


logger = logging.getLogger(__name__)

arg_parser = make_variants_parser(
    prog="varcode genes",
    description="Annotate variants with overlapping gene names")
arg_parser.add_argument(
    "--output-csv", "--output", "-o", dest="output_csv", help="CSV output: chr, start, ref, alt, gene_name, gene_id, plus SV fields")

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
    _configure_logging()
    print_version_info()
    if args_list is None:
        args_list = sys.argv[1:]
    args = arg_parser.parse_args(args_list)
    with exit_on_user_error(arg_parser):
        if args.output_csv:
            check_output_directory(args.output_csv)
        variants = variant_collection_from_args(args)
    variants_dataframe = variants.to_dataframe(raise_on_error=False if args.skip_errors else None)
    logger.info('\n%s', variants_dataframe)
    if args.output_csv:
        with exit_on_user_error(arg_parser, errors=(OSError,)):
            variants_dataframe.to_csv(args.output_csv, index=False)
