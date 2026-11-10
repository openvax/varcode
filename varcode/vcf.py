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
from operator import index
import urllib.parse
import urllib.request
import logging
from collections import OrderedDict
from typing import Literal
import warnings
from warnings import warn

import pandas

from .contigs import _RCRS_LENGTH, _is_grch37
from .pinned_references import (
    find_pinned_name, is_pinned_assembly, names_rcrs, pinned_assembly, pinned_contig_lengths,
    pinned_contig_name, pinned_name)
from .reference import _split_release, infer_genome, infer_reference_name
from .variant import Variant, variant_ascending_position_sort_key
from .variant_collection import VariantCollection
from .vcf_parsing import VCFHeader, _split_struct_fields


logger = logging.getLogger(__name__)


def _is_symbolic_allele(alt):
    """Return True if `alt` is a VCF symbolic allele or breakend notation
    that varcode cannot currently represent as a simple ref/alt Variant.

    Examples of alleles that return True:
        <DEL>, <DUP>, <INS>, <INV>, <CN0>, <INS:ME:ALU>  (symbolic)
        G]17:198982], ]17:198982]G, [13:123456[T, T[13:123456[  (breakends)
        .A, A., .ACGT, ACGT.  (single breakends)
        *  (spanning deletion placeholder)
    """
    if not alt:
        return False
    if alt.startswith("<"):
        return True
    if "[" in alt or "]" in alt:
        return True
    if len(alt) > 1 and (alt[0] == "." or alt[-1] == "."):
        return True
    if alt == "*":
        return True
    return False


def _declares_rcrs_chrm(header):
    """Whether the header's ``##contig`` lines give chrM (or M) only the rCRS
    length, as GENCODE's GRCh37 does; hg19's chrM is 16,571 bp (#710)."""
    lengths = set()
    for line in header.metadata.get("contig", ()):
        line = line.strip()
        if line.startswith("<") and line.endswith(">"):
            fields = _split_struct_fields(line[1:-1])
            if fields.get("ID", "").upper() in ("CHRM", "M"):
                lengths.add(fields.get("length", "").strip())
    return lengths == {str(_RCRS_LENGTH)}


def _names_rcrs_file(genome, reference):
    """Whether the VCF's reference is a file whose chrM is the rCRS: the
    ``genome=`` name if it is a pinned file's, else the header's reference."""
    for name in (genome, reference):
        if isinstance(name, str):
            try:
                found = infer_reference_name(_split_release(name)[0])
            except ValueError:
                continue
            if pinned_name(found) is not None:
                return names_rcrs(found)
    return False


def _check_pinned_contig_lengths(path, header, assembly, chrm_is_rcrs):
    """Raise when ``##contig`` lengths say the VCF isn't on ``assembly``.

    A GRCh38 VCF read as T2T-CHM13v2.0 would put GRCh38 positions on T2T
    transcripts; chr1 is 248,956,422 bp in GRCh38 and 248,387,328 in T2T.
    Contigs the pinned table doesn't name, such as decoys, are not checked.
    """
    lengths = pinned_contig_lengths(assembly, chrm_is_rcrs)
    mismatches = []
    for line in header.metadata.get("contig", ()):
        line = line.strip()
        if not (line.startswith("<") and line.endswith(">")):
            continue
        fields = _split_struct_fields(line[1:-1])
        name = pinned_contig_name(fields.get("ID", ""), assembly, chrm_is_rcrs)
        length = fields.get("length", "").strip()
        if name in lengths and length.isdigit() and int(length) != lengths[name]:
            mismatches.append("%s is %s bp here but %s bp in %s" % (
                fields["ID"], format(int(length), ","), format(lengths[name], ","), assembly))
    if mismatches:
        raise ValueError(
            "%s is not on %s: its ##contig lengths disagree (%s). Annotating it would put another "
            "assembly's positions on %s transcripts; pass the genome its reference matches."
            % (path, assembly, "; ".join(mismatches[:3]) + ("; ..." if len(mismatches) > 3 else ""),
               assembly))


def _check_header_reference(path, reference, assembly):
    """Raise when the header's reference and ``genome=`` put the VCF on two
    assemblies, one of them a genome ORFdata pins, such as T2T-CHM13v2.0.

    A GRCh38 VCF passed with ``genome=`` a T2T genome, or a T2T VCF with
    GRCh38, would annotate one assembly's positions on the other's
    transcripts, and varcode doesn't convert them (#243). A header whose
    reference names no recognized assembly isn't checked (#669).
    """
    if not isinstance(reference, str):
        return
    named = find_pinned_name(reference)
    if named is not None:
        header_assembly = pinned_assembly(named)
    else:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                named = infer_reference_name(reference)
        except ValueError:
            return
        header_assembly = None
        if not is_pinned_assembly(assembly):
            return
    if header_assembly == assembly:
        return
    raise ValueError(
        "%s: its ##reference, %r, names %s, but genome= is %s. One assembly's positions aren't "
        "the other's, and varcode doesn't convert them. Pass the genome the VCF was called on, "
        "or reference_vcf_key=None if the header is wrong." % (path, reference, named, assembly))


def load_vcf(
        path: str | os.PathLike[str],
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
        parse_structural_variants=False,
        genome_fasta=None,
        warn_on_filtered=False) -> VariantCollection:
    """
    Load reference name and Variant objects from the given VCF filename.

    Local files are parsed directly. HTTP/HTTPS URLs are downloaded to a
    temporary file and ``load_vcf`` recurses on the local copy; pandas
    doesn't reliably stream gzipped HTTP responses, so we materialize first.

    Reference-identical ALT alleles are skipped with a counted warning.
    They do not count toward ``max_variants``; retained alleles keep their
    original VCF ALT indices and metadata.

    Parameters
    ----------

    path : str or os.PathLike[str]
        Path to VCF (*.vcf) or compressed VCF (*.vcf.gz).

    genome : {pyensembl.Genome, reference name, Ensembl version int}, optional
        Optionally pass in a PyEnsembl Genome object, name of reference, or
        PyEnsembl release version to specify the reference associated with a
        VCF. Otherwise the reference is inferred from the header line named by
        reference_vcf_key, usually a FASTA path, from an assembly name or
        alias in its file name, or a common build's file name such as
        ``hs37d5``, ``human_g1k_v37`` or the Broad's
        ``Homo_sapiens_assembly38``. A ``genome`` that disagrees with the
        header's reference raises ``ValueError`` when one of the two is a
        genome ORFdata pins, such as T2T-CHM13v2.0: a GRCh38 VCF read as
        T2T, or the reverse.

    reference_vcf_key : str, optional
        Name of metadata field which contains path to reference FASTA
        file (default = 'reference'). With ``genome`` given, None ignores
        the header's reference.

    only_passing : bool, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    warn_on_filtered : bool, default False
        Warn with the number of records skipped by FILTER (before splitting
        multiallelic records). The CLI enables this by default.

    allow_extended_nucleotides : bool, default False
        Allow characters other that A,C,T,G in the ref and alt strings.

    include_info : bool, default True
        Whether to parse the INFO and per-sample columns. If you don't need
        these, set to False for faster parsing.

    chunk_size: int, optional
        Number of records to load in memory at once.

    max_variants : int, optional
        Maximum number of successfully loaded alleles, before sorting or
        deduplication. Applies to small and structural variants. Must be a
        nonnegative integer; zero loads no alleles and None means no limit.

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
        hg19's "chrM" stays "chrM": it is not the rCRS that Ensembl names MT,
        so its variants can't be annotated (see :class:`Variant`). Under
        GRCh37, a header that declares chrM at the rCRS length (16,569 bp, as
        in GENCODE's GRCh37) marks it as the rCRS, so it becomes "MT"; each
        variant then has ``chrm_is_rcrs`` set. Under T2T-CHM13v2.0 the
        reference's file name decides: the ``chm13v2.0_maskedY_rCRS`` set's
        chrM is the rCRS, "NC_012920.1", and any other file's is CHM13's own
        "MT", which has no annotation. ``##contig`` lengths that disagree
        with T2T-CHM13v2.0's, as a GRCh38 VCF's do, raise ``ValueError``.

    genome_fasta : str, path-like, or FASTA object, optional
        Optionally attach a chromosome FASTA to the resolved genome
        before parsing. Equivalent to wrapping with
        ``varcode.Genome(genome, fasta=genome_fasta)`` and passing the
        wrapper; see :class:`varcode.Genome` for the accepted types
        and verification behavior. Without this, features that need
        raw genomic bases (cryptic-exon scoring, indel
        left-alignment, sequence-aware splice prediction) fall back
        to transcript-cDNA coverage only.
    """

    path = os.fspath(path)
    if not isinstance(path, str):
        raise TypeError("Path or URL to VCF must resolve to str, got %s" % type(path).__name__)
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
                genome_fasta=genome_fasta,
                chunk_size=chunk_size,
                max_variants=max_variants,
                sort_key=sort_key,
                distinct=distinct,
                normalize_contig_names=normalize_contig_names,
                convert_ucsc_contig_names=convert_ucsc_contig_names,
                parse_structural_variants=parse_structural_variants,
                warn_on_filtered=warn_on_filtered)
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

    genome_given = genome is not None
    if genome is None:
        if reference_vcf_key not in header.metadata:
            raise ValueError("Unable to infer reference genome for %s" % (path,))
        genome = header.metadata[reference_vcf_key]

    # Preserve the complete UCSC selector, including any explicit release,
    # so variants retain both their reference identity and contig convention.
    original_genome = genome
    genome, genome_was_ucsc = infer_genome(genome)
    assembly = getattr(genome, "assembly", genome.reference_name)
    chrm_is_rcrs = _is_grch37(assembly) and _declares_rcrs_chrm(header)
    if chrm_is_rcrs:
        logger.info("%s declares chrM at the rCRS length; naming it MT", path)
    if genome_given:
        _check_header_reference(path, header.metadata.get(reference_vcf_key), assembly)
    if is_pinned_assembly(assembly):
        # Only the file's name says which mitochondrion a T2T chrM is.
        chrm_is_rcrs = _names_rcrs_file(original_genome, header.metadata.get(reference_vcf_key))
        _check_pinned_contig_lengths(path, header, assembly, chrm_is_rcrs)
    if genome_was_ucsc:
        genome = original_genome

    if convert_ucsc_contig_names is None:
        convert_ucsc_contig_names = genome_was_ucsc

    if genome_fasta is not None:
        # Wrap the resolved genome in varcode.Genome so the FASTA is
        # part of the genome object rather than a side-channel
        # attachment. Late import keeps the cost on callers who opt in.
        from .genome import Genome as VarcodeGenome
        genome = VarcodeGenome(genome, fasta=genome_fasta)

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
    if chrm_is_rcrs:
        variant_kwargs['chrm_is_rcrs'] = True

    variant_collection_kwargs = {
        'sort_key': sort_key,
        'distinct': distinct
    }

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
        parse_structural_variants=parse_structural_variants,
        warn_on_filtered=warn_on_filtered)



def load_vcf_fast(*args, **kwargs):
    """
    Deprecated alias for :func:`load_vcf`; warns with ``DeprecationWarning``.
    """
    warn(
        "load_vcf_fast is deprecated and has been renamed to load_vcf",
        DeprecationWarning, stacklevel=2)
    return load_vcf(*args, **kwargs)


def dataframes_to_variant_collection(
        dataframes,
        source_path,
        info_parser=None,
        only_passing=True,
        max_variants=None,
        sample_names=None,
        sample_info_parser=None,
        variant_kwargs=None,
        variant_collection_kwargs=None,
        parse_structural_variants=False,
        warn_on_filtered=False):
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
        and 'INFO' if `info_parser` is not None. Columns must be in this
        order.

    source_path : str
        Path of VCF file from which DataFrame chunks were generated.

    info_parser : string -> object, optional
        Callable to parse INFO strings.

    only_passing : boolean, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is
        dropped.

    max_variants : int, optional
        Maximum number of successfully loaded alleles, before sorting or
        deduplication. Applies to small and structural variants. Must be a
        nonnegative integer; zero loads no alleles and None means no limit.

    sample_names : list of strings, optional
        Sample names. The final columns of the dataframe should match these.
        If specified, the per-sample info columns will be parsed. You must
        also specify sample_info_parser.

    sample_info_parser : string list * string -> dict, optional
        Callable to parse per-sample info columns.

    variant_kwargs : dict, optional
        Additional keyword parameters to pass to Variant.__init__

    variant_collection_kwargs : dict, optional
        Additional keyword parameters to pass to VariantCollection.__init__.
    """

    if max_variants is not None:
        try:
            max_variants = index(max_variants)
        except TypeError:
            raise ValueError("max_variants must be a nonnegative integer or None") from None
        if max_variants < 0:
            raise ValueError("max_variants must be a nonnegative integer or None")
    if variant_kwargs is None:
        variant_kwargs = {}
    if variant_collection_kwargs is None:
        variant_collection_kwargs = {}

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
    structural_variant_genome = variant_kwargs.get("genome")
    if structural_variant_genome is None:
        structural_variant_genome = variant_kwargs.get("ensembl")
    # Split the skip counters so the user-facing warning can name a fix.
    # ``flag_off``: symbolic ALT seen with ``parse_structural_variants=False``
    #   — the user can opt into SV loading by flipping the flag.
    # ``unparseable``: flag is on but ``parse_symbolic_alt`` couldn't
    #   recognize the ALT shape — typically a malformed or
    #   caller-specific extension; the flag won't help.
    # The spanning-deletion ``*`` placeholder is silent — every consumer
    # drops it and warning every time would just be noise.
    n_skipped_flag_off = 0
    n_skipped_unparseable = 0
    n_skipped_no_change = 0
    n_filtered_records = 0
    limit_reached = max_variants == 0
    for chunk in (() if limit_reached else dataframes):
        if chunk.columns.tolist() != expected_columns:
            raise ValueError(
                "dataframe columns (%s) do not match expected columns (%s)" % (
                    chunk.columns, expected_columns))

        for tpl in chunk.itertuples():
            (i, chrom, pos, id_, ref, alts, qual, flter) = tpl[:8]
            if flter == ".":
                flter = None
            elif flter == "PASS":
                flter = []
            elif only_passing:
                n_filtered_records += 1
                continue
            else:
                flter = flter.split(';')
            if id_ == ".":
                id_ = None
            qual = float(qual) if qual != "." else None
            info = sample_info = None
            for alt_num, alt in enumerate(alts.split(",")):
                # Missing ALT and spanning-deletion placeholders are not
                # independently annotatable alleles and do not count toward
                # the limit. enumerate preserves the original VCF ALT index.
                if alt in (".", "*"):
                    continue
                symbolic = _is_symbolic_allele(alt)
                if symbolic and not parse_structural_variants:
                    n_skipped_flag_off += 1
                    continue
                if info_parser is not None and info is None:
                    info = info_parser(tpl[8])
                    if sample_names:
                        sample_info = sample_info_parser(list(tpl[10:]), tpl[9])

                variant: Variant | None
                if symbolic:
                    from .sv_allele_parser import _confidence_info_for_alt, parse_symbolic_alt
                    variant = parse_symbolic_alt(
                        contig=chrom,
                        start=int(pos),
                        ref=ref,
                        alt=alt,
                        info=_confidence_info_for_alt(info, alt_num, len(alts.split(','))),
                        genome=structural_variant_genome,
                        normalize_contig_names=variant_kwargs.get(
                            "normalize_contig_names", True),
                        convert_ucsc_contig_names=variant_kwargs.get(
                            "convert_ucsc_contig_names"),
                        chrm_is_rcrs=variant_kwargs.get("chrm_is_rcrs", False),
                    )
                    if variant is None:
                        n_skipped_unparseable += 1
                        continue
                else:
                    variant = Variant(chrom, int(pos), ref, alt, **variant_kwargs)
                    # VCF ALT contains non-reference alleles. Compare the
                    # normalized strings so case differences are not edits.
                    # https://samtools.github.io/hts-specs/VCFv4.4.pdf section 1.6
                    if variant.original_ref == variant.original_alt:
                        n_skipped_no_change += 1
                        continue
                variants.append(variant)
                metadata[variant] = {
                    'id': id_,
                    'qual': qual,
                    'filter': flter,
                    'info': info,
                    'sample_info': sample_info,
                    'alt_allele_index': alt_num,
                    'vcf_alt_alleles': tuple(alts.split(',')),
                }
                if max_variants is not None and len(variants) >= max_variants:
                    limit_reached = True
                    break
            if limit_reached:
                break
        if limit_reached:
            break

    if warn_on_filtered and n_filtered_records:
        warn(
            "Skipped %d VCF record(s) with FILTER other than PASS or '.' in %s. "
            "Use --include-filtered on the CLI (only_passing=False in Python) "
            "to include them." % (n_filtered_records, source_path))
    if n_skipped_no_change:
        warn(
            "Skipped %d reference-identical ALT allele(s) in %s (REF == ALT); "
            "these records do not describe a sequence change." % (
                n_skipped_no_change, source_path))
    if n_skipped_flag_off > 0:
        warn(
            "Skipped %d symbolic/breakend allele(s) in %s "
            "(e.g. <DEL>, <CN0>, <INS:ME:ALU>, breakends). "
            "Pass parse_structural_variants=True to load_vcf to "
            "load these as StructuralVariant objects." % (
                n_skipped_flag_off, source_path)
        )
    if n_skipped_unparseable > 0:
        warn(
            "Could not parse %d symbolic/breakend allele(s) in %s — "
            "ALT did not match a recognized symbolic or VCF 4.1 "
            "breakend shape. Tracked in openvax/varcode#264." % (
                n_skipped_unparseable, source_path)
        )

    return VariantCollection(
        variants=variants,
        source_to_metadata_dict={source_path: metadata},
        **variant_collection_kwargs)


def read_vcf_into_dataframe(
        path: str | os.PathLike[str],
        include_info=False,
        sample_names=None,
        chunk_size=None):
    """
    Load the data of a VCF into a pandas dataframe. All headers are ignored.

    Parameters
    ----------
    path : str or os.PathLike[str]
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
    path = os.fspath(path)
    if not isinstance(path, str):
        raise TypeError("Path to VCF must resolve to str, got %s" % type(path).__name__)
    vcf_field_types: OrderedDict[str, type[str] | type[int]] = OrderedDict()
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

    compression: Literal["gzip", "bz2"] | None = None
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


def parse_url_or_path(s: str) -> urllib.parse.ParseResult:
    # urlparse will parse paths with two leading slashes (e.g. "//foo")
    # in a strange way. We collapse these paths to start with just one
    # slash.
    if s.startswith("//"):
        s = "/" + s.lstrip("/")
    return urllib.parse.urlparse(s)
