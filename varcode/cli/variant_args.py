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
from argparse import ArgumentParser

from pyensembl import genome_for_reference_name

from ..vcf import load_vcf
from ..maf import load_maf
from ..variant_collection import VariantCollection
from ..variant import Variant


def add_variant_args(arg_parser):
    """
    Extends an ArgumentParser instance with the following commandline arguments:
        --vcf
        --genome
        --maf
        --variant
        --json-variants
    """
    variant_arg_group = arg_parser.add_argument_group(
        title="Variants",
        description="Genomic variant files")

    variant_arg_group.add_argument(
        "--vcf",
        default=[],
        action="append",
        help="Genomic variants in VCF format")

    variant_arg_group.add_argument(
        "--maf",
        default=[],
        action="append",
        help="Genomic variants in TCGA's MAF format",)

    variant_arg_group.add_argument(
        "--variant",
        default=[],
        action="append",
        nargs=4,
        metavar=("CHR", "POS", "REF", "ALT"),
        help=(
            "Individual variant as 4 arguments giving chromsome, position, ref,"
            " and alt. Example: chr1 3848 C G. Use '.' to indicate empty alleles"
            " for insertions or deletions."))

    variant_arg_group.add_argument(
        "--genome",
        type=str,
        help=(
            "What reference assembly your variant coordinates are using. "
            "Examples: 'hg19', 'GRCh38', or 'mm9'. "
            "This argument is ignored for MAF files, since each row includes "
            "the reference. "
            "For VCF files, this is used if specified, and otherwise is guessed from "
            "the header. For variants specfied on the commandline with --variant, "
            "this option is required."))

    variant_arg_group.add_argument(
        "--json-variants",
        default=[],
        action="append",
        help="Path to Varcode.VariantCollection object serialized as a JSON file.")

    return variant_arg_group


def make_variants_parser(**kwargs):
    """
    Parameters
    ----------
    **kwargs : dict
        Passed directly to argparse.ArgumentParser

    Creates argparse.ArgumentParser instance with options needed for loading
    variants from VCF, MAF, or JSON files.
    """
    parser = ArgumentParser(**kwargs)
    add_variant_args(parser)
    return parser

def variant_collection_from_args(args, required=True):
    variant_collections = []

    if args.genome:
        genome = genome_for_reference_name(args.genome)
    else:
        # no genome specified, assume it can be inferred from the file(s)
        # we're loading
        genome = None

    for vcf_path in args.vcf:
        variant_collections.append(load_vcf(vcf_path, genome=genome))
    for maf_path in args.maf:
        variant_collections.append(load_maf(maf_path))

    if args.variant:
        if not genome:
            raise ValueError(
                "--reference-name must be specified when using --variant")

        variants = [
            Variant(
                chromosome,
                start=position,
                ref=ref,
                alt=alt,
                ensembl=genome)
            for (chromosome, position, ref, alt)
            in args.variant
        ]
        variant_collection = VariantCollection(variants)
        variant_collections.append(variant_collection)

    for json_path in args.json_variants:
        with open(json_path, 'r') as f:
            variant_collections.append(
                VariantCollection.from_json(f.read()))

    if required and len(variant_collections) == 0:
        raise ValueError(
            "No variants loaded (use --maf, --vcf, --variant, or --json-variants options)")
    # pylint: disable=no-value-for-parameter
    return VariantCollection.union(*variant_collections)
