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


import os
import urllib
import logging
from collections import OrderedDict
from warnings import warn

import pandas
from typechecks import require_string

from .reference import (
    infer_genome,
    ensembl_to_ucsc_reference_names
)
from .variant import Variant, variant_ascending_position_sort_key
from .variant_collection import VariantCollection
from .vcf_parsing import VCFHeader


logger = logging.getLogger(__name__)


def _is_symbolic_allele(alt):
    """Return True if `alt` is a VCF symbolic allele or breakend notation
    that varcode cannot currently represent as a simple ref/alt Variant.

    Examples of alleles that return True:
        <DEL>, <DUP>, <INS>, <INV>, <CN0>, <INS:ME:ALU>  (symbolic)
        G]17:198982], ]17:198982]G, [13:123456[T, T[13:123456[  (breakends)
        *  (spanning deletion placeholder)
    """
    if not alt:
        return False
    if alt.startswith("<"):
        return True
    if "[" in alt or "]" in alt:
        return True
    if alt == "*":
        return True
    return False


def load_vcf(
        path,
        genome=None,
        reference_vcf_key="reference",
        only_passing=True,
        allow_extended_nucleotides=False,
        include_info=True,
        chunk_size=10 ** 5,
        max_variants=None,
        sort_key=variant_ascending_position_sort_key,
        distinct=True,
        normalize_contig_names=True,
        convert_ucsc_contig_names=True,
        parse_structural_variants=False):
    """
    Load reference name and Variant objects from the given VCF filename.

    Currently only local files are supported by this function (no http). If you
    call this on an HTTP URL, it will fall back to `load_vcf`.

    Parameters
    ----------

    path : str
        Path to VCF (*.vcf) or compressed VCF (*.vcf.gz).

    genome : {pyensembl.Genome, reference name, Ensembl version int}, optional
        Optionally pass in a PyEnsembl Genome object, name of reference, or
        PyEnsembl release version to specify the reference associated with a
        VCF (otherwise infer reference from VCF using reference_vcf_key)

    reference_vcf_key : str, optional
        Name of metadata field which contains path to reference FASTA
        file (default = 'reference')

    only_passing : bool, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    allow_extended_nucleotides : bool, default False
        Allow characters other that A,C,T,G in the ref and alt strings.

    include_info : bool, default True
        Whether to parse the INFO and per-sample columns. If you don't need
        these, set to False for faster parsing.

    chunk_size: int, optional
        Number of records to load in memory at once.

    max_variants : int, optional
        If specified, return only the first max_variants variants.

    sort_key : fn
        Function which maps each element to a sorting criterion.
        Set to None to not to sort the variants.

    distinct : bool, default True
        Don't keep repeated variants

    normalize_contig_names : bool, default True
        By default contig names will be normalized by converting integers
        to strings (e.g. 1 -> "1"), and converting any letters after "chr"
        to uppercase (e.g. "chrx" -> "chrX"). If you don't want
        this behavior then pass normalize_contig_names=False.

    convert_ucsc_contig_names : bool, default True
        Convert chromosome names from hg19 (e.g. "chr1") to equivalent names
        for GRCh37 (e.g. "1"). By default this is set to True. If None, it
        also evaluates to True if the genome of the VCF is a UCSC reference.
    """

    require_string(path, "Path or URL to VCF")
    parsed_path = parse_url_or_path(path)

    if parsed_path.scheme and parsed_path.scheme.lower() != "file":
        # pandas.read_table nominally supports HTTP, but it tends to crash on
        # large files and does not support gzip. Switching to the python-based
        # implementation of read_table (with engine="python") helps with some
        # issues but introduces a new set of problems (e.g. the dtype parameter
        # is not accepted). For these reasons, we're currently not attempting
        # to load VCFs over HTTP with pandas directly, and instead download it
        # to a temporary file and open that.

        (filename, headers) = urllib.request.urlretrieve(path)
        try:
            # The downloaded file has no file extension, which confuses pyvcf
            # for gziped files in Python 3. We rename it to have the correct
            # file extension.
            new_filename = "%s.%s" % (
                filename, parsed_path.path.split(".")[-1])
            os.rename(filename, new_filename)
            filename = new_filename
            return load_vcf(
                filename,
                genome=genome,
                reference_vcf_key=reference_vcf_key,
                only_passing=only_passing,
                allow_extended_nucleotides=allow_extended_nucleotides,
                include_info=include_info,
                chunk_size=chunk_size,
                max_variants=max_variants,
                sort_key=sort_key,
                distinct=distinct,
                normalize_contig_names=normalize_contig_names,
                convert_ucsc_contig_names=convert_ucsc_contig_names)
        finally:
            logger.info("Removing temporary file: %s", filename)
            os.unlink(filename)

    # Loading a local file.
    # The file will be opened twice: first to parse the header, then by
    # pandas to read the data. Normalize away any file:// scheme so the
    # opener (and pandas, below) get a plain filesystem path.
    if parsed_path.scheme and parsed_path.scheme.lower() == "file":
        path = parsed_path.path
    header = VCFHeader.from_path(path)

    ####
    # The following code looks a bit crazy because it's motivated by the
    # desired to preserve UCSC reference names even though the Variant
    # objects we're creating will convert them to EnsemblRelease genomes
    # with different reference names.
    #
    # For example, if a VCF is aligned against 'hg19' then we want to create a
    # variant which has 'hg19' as its genome argument, so that serialization
    # back to VCF will put the correct reference genome in the generated
    # header.
    if genome is None:
        if reference_vcf_key not in header.metadata:
            raise ValueError("Unable to infer reference genome for %s" % (path,))
        genome = header.metadata[reference_vcf_key]

    genome, genome_was_ucsc = infer_genome(genome)
    if genome_was_ucsc:
        genome = ensembl_to_ucsc_reference_names[genome.reference_name]

    if convert_ucsc_contig_names is None:
        convert_ucsc_contig_names = genome_was_ucsc

    df_iterator = read_vcf_into_dataframe(
        path,
        include_info=include_info,
        sample_names=header.samples if include_info else None,
        chunk_size=chunk_size)

    if include_info:
        sample_info_parser = header.parse_samples
    else:
        sample_info_parser = None

    variant_kwargs = {
        'genome': genome,
        'allow_extended_nucleotides': allow_extended_nucleotides,
        'normalize_contig_names': normalize_contig_names,
        'convert_ucsc_contig_names': convert_ucsc_contig_names,
    }

    variant_collection_kwargs = {
        'sort_key': sort_key,
        'distinct': distinct
    }

    # TODO: drop chrMT variants from hg19 and warn user about it

    return dataframes_to_variant_collection(
        df_iterator,
        source_path=path,
        info_parser=header.parse_info if include_info else None,
        only_passing=only_passing,
        max_variants=max_variants,
        sample_names=header.samples if include_info else None,
        sample_info_parser=sample_info_parser,
        variant_kwargs=variant_kwargs,
        variant_collection_kwargs=variant_collection_kwargs,
        parse_structural_variants=parse_structural_variants)



def load_vcf_fast(*args, **kwargs):
    """
    Same as load_vcf, keeping this name for backwards compatibility.
    """
    warn(
        "load_vcf_fast is deprecated and has been renamed to load_vcf",
        DeprecationWarning)
    return load_vcf(*args, **kwargs)


def dataframes_to_variant_collection(
        dataframes,
        source_path,
        info_parser=None,
        only_passing=True,
        max_variants=None,
        sample_names=None,
        sample_info_parser=None,
        variant_kwargs={},
        variant_collection_kwargs={},
        parse_structural_variants=False):
    """
    Load a VariantCollection from an iterable of pandas dataframes.

    This takes an iterable of dataframes instead of a single dataframe to avoid
    having to load huge dataframes at once into memory. If you have a single
    dataframe, just pass it in a single-element list.

    Parameters
    ----------
    dataframes
        Iterable of dataframes (e.g. a generator). Expected columns are:
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
        and 'INFO' if `info_parser` is not Null. Columns must be in this
        order.

    source_path : str
        Path of VCF file from which DataFrame chunks were generated.

    info_parser : string -> object, optional
        Callable to parse INFO strings.

    only_passing : boolean, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    max_variants : int, optional
        If specified, return only the first max_variants variants.

    sample_names : list of strings, optional
        Sample names. The final columns of the dataframe should match these.
        If specified, the per-sample info columns will be parsed. You must
        also specify sample_info_parser.

    sample_info_parser : string list * string -> dict, optional
        Callable to parse per-sample info columns.

    variant_kwargs : dict, optional
        Additional keyword paramters to pass to Variant.__init__

    variant_collection_kwargs : dict, optional
        Additional keyword parameters to pass to VariantCollection.__init__.
    """

    expected_columns = (
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] +
        (["INFO"] if info_parser else []))

    if info_parser and sample_names:
        if sample_info_parser is None:
            raise TypeError(
                "Must specify sample_info_parser if specifying sample_names")
        expected_columns.append("FORMAT")
        expected_columns.extend(sample_names)

    variants = []
    metadata = {}
    n_skipped_symbolic = 0
    try:
        for chunk in dataframes:
            assert chunk.columns.tolist() == expected_columns,\
                "dataframe columns (%s) do not match expected columns (%s)" % (
                    chunk.columns, expected_columns)

            for tpl in chunk.itertuples():
                (i, chrom, pos, id_, ref, alts, qual, flter) = tpl[:8]
                if flter == ".":
                    flter = None
                elif flter == "PASS":
                    flter = []
                elif only_passing:
                    continue
                else:
                    flter = flter.split(';')
                if id_ == ".":
                    id_ = None
                qual = float(qual) if qual != "." else None
                alt_num = 0
                info = sample_info = None
                for alt in alts.split(","):
                    if alt != ".":
                        if _is_symbolic_allele(alt):
                            # Spanning-deletion placeholder (``*``) is
                            # always skipped — it's a cross-row
                            # reference, not an allele to annotate.
                            if alt == "*":
                                n_skipped_symbolic += 1
                                alt_num += 1
                                continue
                            if parse_structural_variants:
                                # Parse symbolic / breakend ALT into a
                                # StructuralVariant (PR 8). Falls back
                                # to skipping if the parser can't
                                # recognize the ALT shape.
                                from .sv_allele_parser import (
                                    parse_symbolic_alt)
                                if info_parser is not None and info is None:
                                    info = info_parser(tpl[8])
                                sv = parse_symbolic_alt(
                                    contig=chrom,
                                    start=int(pos),
                                    ref=ref,
                                    alt=alt,
                                    info=info,
                                    genome=variant_kwargs.get("ensembl"),
                                )
                                if sv is not None:
                                    variants.append(sv)
                                    metadata[sv] = {
                                        "id": id_,
                                        "qual": qual,
                                        "filter": flter,
                                        "info": info,
                                        "sample_info": sample_info,
                                    }
                                    alt_num += 1
                                    continue
                            # Flag off or parser rejected the ALT:
                            # preserve the legacy filter-and-warn path.
                            n_skipped_symbolic += 1
                            alt_num += 1
                            continue
                        if info_parser is not None and info is None:
                            info = info_parser(tpl[8])  # INFO column
                            if sample_names:
                                # Sample name -> field -> value dict.
                                sample_info = sample_info_parser(
                                    list(tpl[10:]),  # sample info columns
                                    tpl[9],    # FORMAT column
                                )

                        variant = Variant(
                            chrom,
                            int(pos),  # want a Python int not numpy.int64
                            ref,
                            alt,
                            **variant_kwargs)
                        variants.append(variant)
                        metadata[variant] = {
                            'id': id_,
                            'qual': qual,
                            'filter': flter,
                            'info': info,
                            'sample_info': sample_info,
                            'alt_allele_index': alt_num,
                        }
                        if max_variants and len(variants) > max_variants:
                            raise StopIteration
                    alt_num += 1
    except StopIteration:
        pass

    if n_skipped_symbolic > 0:
        warn(
            "Skipped %d symbolic/breakend allele(s) in %s that varcode "
            "cannot currently represent as simple ref/alt Variants "
            "(e.g. <DEL>, <CN0>, <INS:ME:ALU>, breakends). "
            "Tracked in openvax/varcode#264." % (
                n_skipped_symbolic, source_path)
        )

    return VariantCollection(
        variants=variants,
        source_to_metadata_dict={source_path: metadata},
        **variant_collection_kwargs)


def read_vcf_into_dataframe(
        path,
        include_info=False,
        sample_names=None,
        chunk_size=None):
    """
    Load the data of a VCF into a pandas dataframe. All headers are ignored.

    Parameters
    ----------
    path : str
        Path to local file. HTTP and other protocols are not implemented.

    include_info : boolean, default False
        If true, the INFO field is not parsed, but is included as a string in
        the resulting data frame. If false, the INFO field is omitted.

    sample_names: string list, optional
        Sample names. The final columns of the dataframe should match these.
        If specified (and include_info is also specified), the FORMAT and
        per-sample info columns will be included in the result dataframe.

    chunk_size : int, optional
        If buffering is desired, the number of rows per chunk.

    Returns
    ---------
    If chunk_size is None (the default), a dataframe with the contents of the
    VCF file. Otherwise, an iterable of dataframes, each with chunk_size rows.

    """
    vcf_field_types = OrderedDict()
    vcf_field_types['CHROM'] = str
    vcf_field_types['POS'] = int
    vcf_field_types['ID'] = str
    vcf_field_types['REF'] = str
    vcf_field_types['ALT'] = str
    vcf_field_types['QUAL'] = str
    vcf_field_types['FILTER'] = str
    if include_info:
        vcf_field_types['INFO'] = str
        if sample_names:
            vcf_field_types['FORMAT'] = str
            for name in sample_names:
                vcf_field_types[name] = str

    parsed_path = parse_url_or_path(path)
    if not parsed_path.scheme or parsed_path.scheme.lower() == "file":
        path = parsed_path.path
    else:
        raise NotImplementedError("Only local files are supported.")

    compression = None
    if path.endswith(".gz"):
        compression = "gzip"
    elif path.endswith(".bz2"):
        compression = "bz2"

    reader = pandas.read_csv(
        path,
        sep="\t",
        compression=compression,
        comment="#",
        chunksize=chunk_size,
        dtype=vcf_field_types,
        names=list(vcf_field_types),
        usecols=range(len(vcf_field_types)))
    return reader


def parse_url_or_path(s):
    # urlparse will parse paths with two leading slashes (e.g. "//foo")
    # in a strange way. We collapse these paths to start with just one
    # slash.
    if s.startswith("//"):
        s = "/" + s.lstrip("/")
    return urllib.parse.urlparse(s)
