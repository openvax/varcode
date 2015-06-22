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

# required so that 'import vcf' gets the global PyVCF package,
# rather than our local vcf module
from __future__ import absolute_import

import requests
import zlib
import contextlib
import collections

try:
    from urlparse import urlparse  # Python 2
except ImportError:
    from urllib.parse import urlparse  # Python 3

import pandas
from pyensembl import cached_release
import typechecks
import vcf  # PyVCF

from .reference_name import (
    infer_reference_name,
    ensembl_release_number_for_reference_name
)
from .variant import Variant
from .variant_collection import VariantCollection

def load_vcf(
        path,
        only_passing=True,
        ensembl_version=None,
        reference_name=None,
        reference_vcf_key="reference",
        allow_extended_nucleotides=False,
        include_info=True,
        chunk_size=1e5,
        max_variants=None):
    '''
    Load reference name and Variant objects from the given VCF filename.

    For local files, this uses a fast but somewhat brittle pandas
    implementation. If you encounter issues with it, try the slower but more
    tested `load_vcf_with_pyvcf` function (and file a bug).

    For VCFs accessed over HTTP, this currently falls back to the
    `load_vcf_with_pyvcf` function.

    Parameters
    ----------

    path : str
        Path or URL to VCF (*.vcf) or compressed VCF (*.vcf.gz). Supported URL 
        schemes are "file", "http", "https", and "ftp".

    only_passing : boolean, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    ensembl_version : int, optional
        Which release of Ensembl to use for annotation, by default inferred
        from the reference path. If specified, then `reference_name` and
        `reference_vcf_key` are ignored.

    reference_name : str, optional
        Name of reference genome against which variants from VCF were aligned.
        If specified, then `reference_vcf_key` is ignored.

    reference_vcf_key : str, optional
        Name of metadata field which contains path to reference FASTA
        file (default = 'reference')

    allow_extended_nucleotides : boolean, default False
        Allow characters other that A,C,T,G in the ref and alt strings.

    include_info : boolean, default True
        Whether to parse the info column. If you don't need that column, set to
        False to get a ~4X speed improvement.

    chunk_size: int, optional
        Number of records to load in memory at once.

    max_variants : int, optional
        If specified, return only the first max_variants variants.
    '''

    typechecks.require_string(path, "Path or URL to VCF")
    parsed_path = urlparse(path)

    if parsed_path.scheme.lower() in ("http", "https", "ftp"):
        # Pandas uses a different implementation of read_table for HTTP, which
        # results in issues for us. For now, we fall back to pyvcf to read over
        # HTTP.
        # TODO: use pandas implementation when loading VCFs over HTTP.
        return load_vcf_with_pyvcf(
            path,
            only_passing=only_passing,
            ensembl_version=ensembl_version,
            reference_name=reference_name,
            reference_vcf_key=reference_vcf_key,
            allow_extended_nucleotides=allow_extended_nucleotides,
            max_variants=max_variants)
    else:
        # Local file.
        # We use pyvcf to parse the header, then pandas for the rest.
        wrapper = PyVCFReaderFromPathOrURL(path)
        vcf_reader = wrapper.vcf_reader
        wrapper.close()  # only need metadata, which is read immediately.
        path = parsed_path.path

    ensembl = make_ensembl(
        vcf_reader, ensembl_version, reference_name, reference_vcf_key)

    vcf_field_types = collections.OrderedDict()
    vcf_field_types['CHROM'] = str
    vcf_field_types['POS'] = int
    vcf_field_types['ID'] = str
    vcf_field_types['REF'] = str
    vcf_field_types['ALT'] = str
    vcf_field_types['QUAL'] = str
    vcf_field_types['FILTER'] = str
    if include_info:
        vcf_field_types['INFO'] = str

    reader = pandas.read_table(
        path,
        compression='gzip' if parsed_path.path.endswith('.gz') else None,
        comment="#",
        chunksize=chunk_size,
        dtype=vcf_field_types,
        names=list(vcf_field_types),
        usecols=range(len(vcf_field_types)))

    variants = []
    metadata = {}
    try:
        for chunk in reader:
            if include_info:
                def g():
                    for tpl in chunk.itertuples():
                        info = vcf_reader._parse_info(tpl[-1])
                        yield tpl[:-1] + (info,)
            else:
                def g():
                    for tpl in chunk.itertuples():
                        yield tpl + (None,)
            for (i, chrom, pos, id_, ref, alts, qual, flter, info) in g():
                flter = None if flter == "." else (
                    [] if flter == 'PASS' else flter.split(";"))
                if only_passing and flter:
                    continue
                for alt in alts.split(","):
                    if alt == ".":
                        continue
                    variant = Variant(
                        chrom,
                        int(pos),
                        ref,
                        alt,
                        ensembl=ensembl,
                        allow_extended_nucleotides=allow_extended_nucleotides)
                    variants.append(variant)
                    metadata[variant] = {
                        'id': None if id_ == "." else id_,
                        'qual': None if qual == "." else float(qual),
                        'filter': flter,
                        'info': info,
                    }
                    if max_variants and len(variants) > max_variants:
                        raise StopIteration
    except StopIteration:
        pass

    return VariantCollection(
        variants=variants,
        path=path,
        metadata=metadata)

def load_vcf_with_pyvcf(
        path,
        only_passing=True,
        ensembl_version=None,
        reference_name=None,
        reference_vcf_key="reference",
        allow_extended_nucleotides=False,
        max_variants=None):
    """
    Load reference name and Variant objects from the given VCF filename.

    This uses a PyVCF reader. It is slower than the pandas implementation in
    load_vcf, but handles more unusual VCF files (e.g. records with structural
    variants) and is also more tolerant of minor deviations from the VCF spec
    (e.g. will accept VCFs that are space-delimited, whereas the panda
    implementation requires tabs).

    Parameters
    ----------

    path : str
        Path or URL to VCF (*.vcf) or compressed VCF (*.vcf.gz). Supported URL 
        schemes are "file", "http", "https", and "ftp".

    only_passing : boolean, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    ensembl_version : int, optional
        Which release of Ensembl to use for annotation, by default inferred
        from the reference path. If specified, then `reference_name` and
        `reference_vcf_key` are ignored.

    reference_name : str, optional
        Name of reference genome against which variants from VCF were aligned.
        If specified, then `reference_vcf_key` is ignored.

    reference_vcf_key : str, optional
        Name of metadata field which contains path to reference FASTA
        file (default = 'reference')

    allow_extended_nucleotides : boolean, default False
        Allow characters other that A,C,T,G in the ref and alt strings.

    max_variants : int, optional
        If specified, return only the first max_variants variants.
    """

    typechecks.require_string(path, "Path or URL to VCF")

    wrapper = PyVCFReaderFromPathOrURL(path)
    with contextlib.closing(wrapper):
        vcf_reader = wrapper.vcf_reader

        ensembl = make_ensembl(
            vcf_reader, ensembl_version, reference_name, reference_vcf_key)

        variants = []
        metadata = {}
        try:
            for record in vcf_reader:
                if (not only_passing or
                        not record.FILTER or
                        record.FILTER == "PASS"):
                    for alt in record.ALT:
                        # We ignore "no-call" variants, i.e. those where
                        # X.ALT = [None]
                        if not alt:
                            continue
                        variant = Variant(
                            contig=record.CHROM,
                            start=record.POS,
                            ref=record.REF,
                            alt=alt.sequence,
                            ensembl=ensembl,
                            allow_extended_nucleotides=allow_extended_nucleotides)
                        variants.append(variant)
                        metadata[variant] = {
                            "id": record.ID,
                            "qual": record.QUAL,
                            "filter": record.FILTER,
                            "info": dict(record.INFO),
                        }
                        if max_variants and len(variants) > max_variants:
                            raise StopIteration
        except StopIteration:
            pass
        return VariantCollection(
            variants=variants,
            path=path,
            metadata=metadata)


class PyVCFReaderFromPathOrURL(object):
    '''
    Thin wrapper over a PyVCF Reader object that supports loading over URLs,
    and a close() function (pyvcf somehow doesn't have a close() funciton).
    '''
    def __init__(self, path):
        typechecks.require_string(path, "Path or URL to VCF")

        self.vcf_reader = None
        self.to_close = None

        parsed_path = urlparse(path)
        if not parsed_path.scheme or parsed_path.scheme.lower() == 'file':
            self.vcf_reader = vcf.Reader(filename=parsed_path.path)
        elif parsed_path.scheme.lower() in ("http", "https", "ftp"):
            self.to_close = response = requests.get(path, stream=True)
            response.raise_for_status()  # raise error on 404, etc.
            if path.endswith(".gz"):
                lines = PyVCFReaderFromPathOrURL._stream_gzip_decompress_lines(
                    response.iter_content())
            else:
                lines = response.iter_lines(decode_unicode=True)
            self.vcf_reader = vcf.Reader(fsock=lines, compressed=False)
        else:
            raise ValueError("Unsupported scheme: %s" % parsed_path.scheme)

    def close(self):
        if self.to_close is not None:
            self.to_close.close()

    @staticmethod
    def _stream_gzip_decompress_lines(stream):
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

def make_ensembl(
        vcf_reader,
        ensembl_version,
        reference_name, 
        reference_vcf_key):
    '''
    Helper function to make an ensembl instance.
    '''
    
    if not ensembl_version:
        if reference_name:
            # normalize the reference name in case it's in a weird format
            reference_name = infer_reference_name(reference_name)
        elif reference_vcf_key not in vcf_reader.metadata:
            raise ValueError("Unable to infer reference genome for %s" % (
                vcf_reader.filename,))
        else:
            reference_path = vcf_reader.metadata[reference_vcf_key]
            reference_name = infer_reference_name(reference_path)
        ensembl_version = ensembl_release_number_for_reference_name(
            reference_name)
    return cached_release(ensembl_version)
