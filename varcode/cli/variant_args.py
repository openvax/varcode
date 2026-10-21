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

from argparse import ArgumentParser
from contextlib import contextmanager
import os

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
        "--genome-fasta",
        help="Chromosome DNA FASTA for this reference assembly; enables intronic/SV sequence analysis")

    variant_arg_group.add_argument(
        "--vcf",
        default=[],
        action="append",
        help="Genomic variants in VCF format, including SVs (repeatable)")

    variant_arg_group.add_argument(
        "--include-filtered", action="store_true",
        help="Include VCF records with FILTER values other than PASS or '.'")
    variant_arg_group.add_argument(
        "--skip-errors", action="store_true",
        help="Continue after annotation errors, retaining Failure rows in the output")

    variant_arg_group.add_argument(
        "--maf",
        default=[],
        action="append",
        help="Genomic variants in TCGA's MAF format (repeatable)",)

    variant_arg_group.add_argument(
        "--variant",
        default=[],
        action="append",
        nargs=4,
        metavar=("CHR", "POS", "REF", "ALT"),
        help=(
            "Individual variant as 4 arguments giving chromosome, position, ref,"
            " and alt. Example: chr1 3848 C G. Use '.' to indicate empty alleles"
            " for insertions or deletions. Repeatable."))

    variant_arg_group.add_argument(
        "--genome",
        type=str,
        help=(
            "What reference assembly your variant coordinates are using. "
            "Examples: 'hg19', 'GRCh38', or 'GRCm38'. "
            "Pin a numbered release or annotation date, e.g. 'GRCh38:93', "
            "'GRCh38:2026_04' or 'GRCh37:2013_09'; "
            "otherwise the most recent installed numbered release is used. "
            "This argument is ignored for MAF files, since each row includes "
            "the reference. "
            "For VCF files, this is used if specified, and otherwise is guessed from "
            "the header. For variants specified on the commandline with --variant, "
            "this option is required."))

    variant_arg_group.add_argument(
        "--download-reference-genome-data",
        action="store_true",
        default=False,
        help=(
            ("Automatically download genome reference data required for "
             "annotation using PyEnsembl. Otherwise you must first run "
             "'pyensembl install' for the release/species corresponding "
             "to the genome used in your VCF.")))

    variant_arg_group.add_argument(
        "--json-variants",
        default=[],
        action="append",
        help="Path to a VariantCollection JSON file (repeatable).")

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


@contextmanager
def exit_on_user_error(arg_parser, errors=(ValueError, OSError)):
    """
    Report input and output errors as a one-line ``prog: error: ...`` message
    with exit status 1. Callers keep annotation outside this context so
    internal annotation errors retain their traceback.
    """
    try:
        yield
    except errors as e:
        arg_parser.exit(1, "%s: error: %s\n" % (
            arg_parser.prog, " ".join(str(e).splitlines())))


def check_output_directory(path):
    """
    Fail before any annotation work if an output file can't be written where
    the user asked for it.
    """
    if os.path.isdir(path):
        raise ValueError("output path is a directory: %s" % path)
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        raise ValueError("output directory does not exist: %s" % directory)


# Genomic coordinates are 1-based. PyEnsembl looks up loci in SQLite, whose
# INTEGER is signed 64-bit; larger values raise OverflowError in the query.
MAX_POSITION = 2 ** 63 - 1


def parse_position(position):
    try:
        parsed = int(position)
    except ValueError:
        raise ValueError(
            "--variant position must be an integer, got %r" % (position,)
        ) from None
    if not 1 <= parsed <= MAX_POSITION:
        raise ValueError(
            "--variant position is out of range, got %s (must be between 1 and %d)" % (
                parsed, MAX_POSITION))
    return parsed


def download_and_install_reference_data(variant_collections):
    unique_genomes = {
        variant.ensembl
        for variant_collection in variant_collections
        for variant in variant_collection
    }
    for genome in unique_genomes:
        if not genome.required_local_files_exist():
            genome.download()
            genome.index()


def variant_collection_from_args(args, required=True):
    variant_collections = []
    fasta = getattr(args, "genome_fasta", None)
    if fasta is not None and not os.path.isfile(fasta):
        raise ValueError("--genome-fasta is not a readable FASTA file: %s" % fasta)

    for vcf_path in args.vcf:
        variant_collections.append(
            load_vcf(
                vcf_path, genome=args.genome, parse_structural_variants=True,
                only_passing=not args.include_filtered, warn_on_filtered=True))

    for maf_path in args.maf:
        variant_collections.append(load_maf(maf_path))

    if args.variant:
        if not args.genome:
            raise ValueError(
                "--genome must be specified when using --variant")

        variants = [
            Variant(
                chromosome,
                start=parse_position(position),
                ref=ref,
                alt=alt,
                genome=args.genome,
                convert_ucsc_contig_names=True)
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

    if fasta is not None:
        from ..genome import Genome
        assemblies = {v.reference_name for collection in variant_collections for v in collection}
        if len(assemblies) > 1:
            raise ValueError("--genome-fasta requires inputs from one reference assembly")
        genomes = {}
        wrapped = []
        for collection in variant_collections:
            replacements = {}
            for variant in collection:
                key = id(variant.genome)
                if key not in genomes:
                    genomes[key] = Genome(variant.genome, fasta=fasta)
                replacements[variant] = type(variant).from_dict({
                    **variant.to_dict(), "genome": genomes[key]})
            wrapped.append(VariantCollection(
                [replacements[v] for v in collection],
                sources=collection.sources,
                source_to_metadata_dict={
                    source: {replacements[v]: metadata for v, metadata in entries.items() if v in replacements}
                    for source, entries in collection.source_to_metadata_dict.items()}))
        variant_collections = wrapped

    if args.download_reference_genome_data:
        download_and_install_reference_data(variant_collections)

    # pylint: disable=no-value-for-parameter
    return VariantCollection.union(*variant_collections)
