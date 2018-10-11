# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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
from collections import defaultdict
import sys


def variants_to_vcf(variants, variant_to_metadata, out=sys.stdout):
    """Output a VCF file from a list of Variant records.

    Parameters
    ----------
    variants : iterable
        Variant objects

    variant_to_metadata : dict
        Dictionary mapping each variant in `variants` to a dictionary
        of metadata.
    """

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
    #
    # As a side note, adding headers for certain fields will make them parse
    # correctly -- into an integer instead of a string (the default), for
    # example.

    # TODO: If we maintain headers (see above TODO), what should happen if
    # variant sources use different reference genomes?
    #
    # If we don't maintain headers, what should the default reference genome
    # be? This code chose one fairly arbitrarily.

    # TODO: If we end up needing more functions to "build" VCF record fields
    # (see functions below), we may want to abstract away the individual
    # functions and create a mapping from field to format function.

    def get_metadata_field(key, variant, default='.'):
        """Retrieve field from variant metadata dictionary."""
        val = variant_to_metadata[variant].get(key)
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

        return ';'.join(build_info_pair(k, v) for (k, v) in info_dict.items())

    def build_format_field(variant):
        """Build the sample format string from the given variant.

        Each sample column follows this format for the specified variant.
        """
        sample_info = get_metadata_field('sample_info', variant, default={})
        return ':'.join(list(sample_info.values())[0].keys()) if sample_info else '.'

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

    def merge_duplicate_variants():
        """Merge duplicate variants (according to ID) and return *list* of merged, sorted variants.

        Multiple `Variant`s can have the same VCF id (e.g., those variants which
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
            # contigs, positions, or reference bases. If merging isn't possible (which
            # I think is the case for many situations), it may be easiest to assign
            # the "missing value" to the IDs and move on.

            assert len(duplicates) > 0
            assert all(v.original_contig == duplicates[0].original_contig and
                       v.original_start == duplicates[0].original_start and
                       v.original_ref == duplicates[0].original_ref for v in duplicates)
            import copy
            merged = copy.copy(duplicates[0])
            merged.original_alt = ','.join(duplicate.original_alt for duplicate in duplicates)
            return merged

        id2variants = construct_id2variants()
        variants_no_id = id2variants.pop('.', [])  # don't want to merge variants w/ no id
        merged_variants = list(map(merge_variant_list, id2variants.values())) + variants_no_id

        # groups variants by contig; relative ordering of contigs doesn't matter
        return sorted(
            merged_variants,
            key=lambda v: (str(v.original_contig), str(v.original_start)))

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

    unique_variants_list = merge_duplicate_variants()

    # usually we have just one reference genome for the whole variant collection
    # but if variants from multiple sources have been merged then we might
    # not be able to write out a VCF since the individual variants may be using
    # different coordinate systems
    genome_names = list({v.ensembl.reference_name for v in unique_variants_list})
    if len(genome_names) > 1:
        raise ValueError(
            "Cannot create VCF for variants with multiple reference genomes: %s" % (
                genome_names,))
    genome_name = genome_names[0]
    print('##fileformat=VCFv4.2', file=out)
    print('##reference=%s' % genome_name, file=out)
    print('\t'.join(headers), file=out)
    for variant in unique_variants_list:
        record = build_vcf_record(variant, add_sample_info=bool(sample_names))
        print('\t'.join(record), file=out)
