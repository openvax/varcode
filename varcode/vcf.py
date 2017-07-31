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


from __future__ import absolute_import, print_function, division
import os
import requests
import zlib
import logging
from collections import OrderedDict
from collections import defaultdict
from warnings import warn

from six.moves import urllib
import pandas
from typechecks import require_string
import vcf as pyvcf

from .reference import infer_genome
from .variant import Variant, variant_ascending_position_sort_key
from .variant_collection import VariantCollection


logger = logging.getLogger(__name__)


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
        distinct=True):
    """
    Load reference name and Variant objects from the given VCF filename.

    This is an experimental faster implementation of `load_vcf`. It is
    typically about 2X faster, and with `include_info=False`, about 4X faster.
    If most of the records in the VCF have failed filters (and
    only_passing=True), this function can be orders of magnitude faster than
    `load_vcf`.

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

    only_passing : boolean, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    allow_extended_nucleotides : boolean, default False
        Allow characters other that A,C,T,G in the ref and alt strings.

    include_info : boolean, default True
        Whether to parse the INFO and per-sample columns. If you don't need
        these, set to False for faster parsing.

    chunk_size: int, optional
        Number of records to load in memory at once.

    max_variants : int, optional
        If specified, return only the first max_variants variants.

    sort_key : fn
        Function which maps each element to a sorting criterion.
        Set to None to not to sort the variants.

    distinct : boolean, default True
        Don't keep repeated variants
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
                distinct=distinct)
        finally:
            logger.info("Removing temporary file: %s", filename)
            os.unlink(filename)

    # Loading a local file.
    # The file will be opened twice: first to parse the header with pyvcf, then
    # by pandas to read the data.

    # PyVCF reads the metadata immediately and stops at the first line with
    # data. We can close the file after that.
    handle = PyVCFReaderFromPathOrURL(path)
    handle.close()
    genome = infer_genome_from_vcf(
        genome,
        handle.vcf_reader,
        reference_vcf_key)

    df_iterator = read_vcf_into_dataframe(
        path,
        include_info=include_info,
        sample_names=handle.vcf_reader.samples if include_info else None,
        chunk_size=chunk_size)

    if include_info:
        def sample_info_parser(unparsed_sample_info_strings, format_string):
            """
            Given a format string like "GT:AD:ADP:DP:FS"
            and a list of sample info strings where each entry is like
            "0/1:3,22:T=3,G=22:25:33", return a dict that maps:
            sample name -> field name -> value. Uses pyvcf to parse the fields.
            """
            return pyvcf_calls_to_sample_info_list(
                handle.vcf_reader._parse_samples(
                    unparsed_sample_info_strings, format_string, None))
    else:
        sample_info_parser = None

    return dataframes_to_variant_collection(
        df_iterator,
        source_path=path,
        info_parser=handle.vcf_reader._parse_info if include_info else None,
        only_passing=only_passing,
        max_variants=max_variants,
        sample_names=handle.vcf_reader.samples if include_info else None,
        sample_info_parser=sample_info_parser,
        variant_kwargs={
            'ensembl': genome,
            'allow_extended_nucleotides': allow_extended_nucleotides},
        variant_collection_kwargs={
            'sort_key': sort_key,
            'distinct': distinct})

def variant_collection_to_vcf(variants):
    """Output a VCF file whose records correspond to those in the given `VariantCollection`."""

    # TODO: The variant metadata dictionary (in the `VariantCollection`)
    # contains different data depending on the original input file format (VCF,
    # MAF, CSV).  It's easy to output variants originally in VCF format, but we
    # might want to consider how fields from MAF (for example) map to those in
    # VCF.

    # TODO: The VCF file we output doesn't contain any VCF metadata headers, as
    # the original headers were thrown away when the VCF file was parsed. We
    # may want to keep some of that information and/or infer some of the
    # headers based on the variants themselves. The former is difficult because
    # merge conflicts will inevitably occur; the latter is difficult because
    # the variants themselves don't contain all the information required for
    # these metadata headers (e.g., descriptions).

    # TODO: If we end up needing more functions to "build" VCF record fields
    # (see functions below), we may want to abstract away the individual
    # functions and create a mapping from field to format function.

    # TODO: Currently fails when there is > 1 source for the variant collection.
    # Fix that. Once duplicates are merged, all is well, but we need to access
    # metadata fields in order to merge them, which requires knowing the source.

    def get_metadata_field(key, variant, default='.'):
        """Retrieve field from variant metadata dictionary."""
        val = variants.metadata[variant].get(key)
        if val is None:
            return default
        return val

    def build_filter_field(variant):
        """Build the filter field from the given variant.

        The `filter` field, as stored in the variants metadata dictionary,
        comes in 3 flavors:

            - empty list: 1+ filters were run and none failed
            - non-empty list: 1+ filters were run and 1+ failed
            - `None`: no filters were run

        This function maps each of these internal representations to their
        corresponding VCF representations.
        """
        filter_metadata = get_metadata_field('filter', variant)
        if type(filter_metadata) == list:
            return 'PASS' if filter_metadata == [] else ';'.join(filter_metadata)
        else:
            # TODO: Can the filter field ever be something other than the 3
            # possibilities described in this function's doc comment?
            bad_value_msg = (
                'Filter metadata field took on unexpected value `{}`. Update '
                'code to account for this value.').format(str(filter_metadata))
            assert filter_metadata == '.', bad_value_msg
            return filter_metadata

    def build_info_field(variant):
        """Build the info field from the given variant.

        Format is `<key>=<val>,...;<key>=<val>,...;<key>=<val>,...`.
        """

        def build_info_pair(key, val):
            """Build key/val pair for info field."""
            # Note: Different from `val == True`, which returns True when `val == 1`.
            if val is True:
                return key

            if type(val) == list:
                val = ','.join(map(str, val))
            else:
                val = str(val)
            return '{}={}'.format(key, val)

        info_dict = get_metadata_field('info', variant, default={})
        if not info_dict:
            return '.'

        return ';'.join(build_info_pair(k, v) for k,v in info_dict.iteritems())

    def build_format_field(variant):
        """Build the sample format string from the given variant.

        Each sample column follows this format for the specified variant.
        """
        sample_info = get_metadata_field('sample_info', variant, default={})
        return ':'.join(sample_info.values()[0].keys()) if sample_info else '.'

    def build_sample_fields(variant):
        """Build the sample fields for the given variant."""

        def build_sample_field(sample):
            """Build a specific sample's field (i.e., one sample column)."""
            sample_vals = sample_info[sample].values()
            return ':'.join(build_sample_val(val) for val in sample_vals)

        def build_sample_val(sample_val):
            """Build a specific value for a sample (i.e., one value for one column)."""
            if type(sample_val) is list:
                return ','.join(map(str, sample_val))
            elif sample_val is not None:
                return str(sample_val)
            else:
                return '.'

        sample_info = get_metadata_field('sample_info', variant, default={})
        return list(build_sample_field(sample) for sample in sample_info)

    def build_vcf_record(variant, add_sample_info):
        """Return a list of all the variant's VCF record fields."""
        record = [
            str(variant.original_contig),
            str(variant.original_start),
            get_metadata_field('id', variant),
            variant.original_ref,
            variant.original_alt,
            str(get_metadata_field('qual', variant)),
            build_filter_field(variant),
            build_info_field(variant),
        ]
        if add_sample_info:
            record.append(build_format_field(variant))
            record.extend(build_sample_fields(variant))
        return record

    def merge_duplicate_variants(variants):
        """Merge duplicate variants (according to ID) and return *list* of merged, sorted variants.

        Multiple `Variant`s can have the same VCF id (e.g., those variants which which
        originally had > 1 alternate base), but we can't output a VCF file with multiple
        records having the same id. This function merges those duplicate variants.
        """

        def construct_id2variants():
            """Construct dictionary which maps variant IDs to `Variant`s."""
            id2variants = defaultdict(list)
            for variant in variants:
                # Note: For variants without IDs (e.g., those that came from
                # MAF files), we assign the missing value.
                variant_id = get_metadata_field('id', variant)
                id2variants[variant_id].append(variant)
            return id2variants

        def merge_variant_list(duplicates):
            """Merge duplicate variant list into one."""

            # TODO: Currently assumes that only alternate bases differ, but we may
            # want or need to merge variants that have the same ID but different
            # contigs, positions, or reference bases.

            assert len(duplicates) > 0
            assert all(v.original_contig == duplicates[0].original_contig and
                       v.original_start == duplicates[0].original_start and
                       v.original_ref == duplicates[0].original_ref for v in duplicates)
            import copy
            merged = copy.copy(duplicates[0])
            for variant in duplicates[1:]:
                merged.original_alt += ',' + variant.original_alt
            return merged

        id2variants = construct_id2variants()
        variants_no_id = id2variants.pop('.', [])  # don't want to merge variants w/ no id
        merged_variants = map(merge_variant_list, id2variants.values()) + variants_no_id

        # groups variants by contig; relative ordering of contigs doesn't matter
        return sorted(merged_variants, key=lambda v: (v.original_contig, v.original_start))


    def get_sample_names():
        """Return the sample names for all variants."""

        # TODO: For now, ensures that every variant has the same samples. If they didn't,
        # we'd have to use the missing value for VCF records with no data on a given
        # sample.  In and of itself, that isn't a problem until the format field contains
        # GT (genotype), which is required for all samples if present in the format field.
        #
        # A couple of ways to handle this:
        #
        #   1) Ignore this requirement of VCF files. I'd assume this would still be
        #   compatible with varcode, but may not be with other VCF parsers.
        #   2) Remove the GT field from those records where this rule would be violated.
        #   This is sad because we lose information.
        #
        # It's also important to note that if/when we combine all samples, we'll have to
        # add a `.` (missing value) for those samples in which a variant has no data.
        # Check the VCF spec to make sure this is valid; if not, we may have to write
        # `.:.:.:.` -- one `.` for each field in the format string.
        #
        # Moreover, we'll want to note whether we care about maintaining the relative
        # ordering of the sample names in the original VCF files. This probably isn't
        # important, but the code below does not do this (it effectively alphabetizes the
        # sample names because of the use of `set`).

        sample_names = set()
        for variant in variants:
            sample_info = get_metadata_field('sample_info', variant, default={})
            addl_sample_names = set(sample_info.keys())

            # Ensures all variants have the same samples.
            if sample_names and sample_names != addl_sample_names:
                return []

            sample_names.update(addl_sample_names)

        return list(sample_names)

    headers = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    sample_names = get_sample_names()
    if sample_names:
        headers += ['FORMAT'] + sample_names

    unique_variants_list = merge_duplicate_variants(variants)

    print('##fileformat=VCFv4.2')
    print('\t'.join(headers))
    for variant in unique_variants_list:
        record = build_vcf_record(variant, add_sample_info=bool(sample_names))
        print('\t'.join(record))

def load_vcf_fast(*args, **kwargs):
    """
    Same as load_vcf, keeping this name for backwards compatibility.
    """
    warn(
        "load_vcf_fast is deprecated and has been renamed to load_vcf",
        DeprecationWarning)
    return load_vcf(*args, **kwargs)

def pyvcf_calls_to_sample_info_list(calls):
    """
    Given pyvcf.model._Call instances, return a dict mapping each sample
    name to its per-sample info:
        sample name -> field -> value
    """
    return OrderedDict(
        (call.sample, call.data._asdict()) for call in calls)

def dataframes_to_variant_collection(
        dataframes,
        source_path,
        info_parser=None,
        only_passing=True,
        max_variants=None,
        sample_names=None,
        sample_info_parser=None,
        variant_kwargs={},
        variant_collection_kwargs={}):
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

    reader = pandas.read_table(
        path,
        compression=compression,
        comment="#",
        chunksize=chunk_size,
        dtype=vcf_field_types,
        names=list(vcf_field_types),
        usecols=range(len(vcf_field_types)))
    return reader


class PyVCFReaderFromPathOrURL(object):
    """
    Thin wrapper over a PyVCF Reader object that supports loading over URLs,
    and a close() function (pyvcf somehow doesn't have a close() funciton).

    Attributes
    ----------
    path : string or None
        path that was loaded, if available.

    vcf_reader : pyvcf Reader instance
    """
    def __init__(self, path):
        """
        Construct a new wrapper.

        Parameters
        ----------
        path : string or pyvcf Reader instance
            Path or URL to load, or Reader instance.
        """
        self.path = None  # string path, if available.
        self.vcf_reader = None  # vcf_reader. Will always be set.
        self._to_close = None  # object to call close() on when we're done.

        if isinstance(path, pyvcf.Reader):
            self.vcf_reader = path
        else:
            require_string(path, "Path or URL to VCF")
            self.path = path
            parsed_path = parse_url_or_path(path)
            if not parsed_path.scheme or parsed_path.scheme.lower() == 'file':
                self.vcf_reader = pyvcf.Reader(
                    filename=parsed_path.path,
                    strict_whitespace=True)
            elif parsed_path.scheme.lower() in ("http", "https", "ftp"):
                self._to_close = response = requests.get(path, stream=True)
                response.raise_for_status()  # raise error on 404, etc.
                if path.endswith(".gz"):
                    lines = stream_gzip_decompress_lines(
                        response.iter_content())
                else:
                    lines = response.iter_lines(decode_unicode=True)
                self.vcf_reader = pyvcf.Reader(
                    fsock=lines,
                    compressed=False,
                    strict_whitespace=True)
            else:
                raise ValueError("Unsupported scheme: %s" % parsed_path.scheme)

    def close(self):
        if self._to_close is not None:
            self._to_close.close()


def stream_gzip_decompress_lines(stream):
    """
    Uncompress a gzip stream into lines of text.

    Parameters
    ----------
    Generator of chunks of gzip compressed text.

    Returns
    -------
    Generator of uncompressed lines.
    """
    dec = zlib.decompressobj(zlib.MAX_WBITS | 16)
    previous = ""
    for compressed_chunk in stream:
        chunk = dec.decompress(compressed_chunk).decode()
        if chunk:
            lines = (previous + chunk).split("\n")
            previous = lines.pop()
            for line in lines:
                yield line
    yield previous


def infer_genome_from_vcf(genome, vcf_reader, reference_vcf_key):
    """
    Helper function to make a pyensembl.Genome instance.
    """
    if genome:
        return infer_genome(genome)
    elif reference_vcf_key not in vcf_reader.metadata:
        raise ValueError("Unable to infer reference genome for %s" % (
            vcf_reader.filename,))
    else:
        reference_path = vcf_reader.metadata[reference_vcf_key]
        return infer_genome(reference_path)


def parse_url_or_path(s):
    # urlparse will parse paths with two leading slashes (e.g. "//foo")
    # in a strange way. We collapse these paths to start with just one
    # slash.
    if s.startswith("//"):
        s = "/" + s.lstrip("/")
    return urllib.parse.urlparse(s)
