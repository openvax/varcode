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
        reference_vcf_key="reference"):
    """
    Load reference name and Variant objects from the given VCF filename.
    Drop any entries whose FILTER field is not one of "." or "PASS".

    Parameters
    ----------

    path : str

    only_passing : boolean, optional
        If true, any entries whose FILTER field is not one of "." or "PASS" is dropped.

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

    """

    typechecks.require_string(path, "Path to VCF")

    vcf_reader = vcf.Reader(filename=path)

    if not ensembl_version:
        if reference_name:
            # normalize the reference name in case it's in a weird format
            reference_name = infer_reference_name(reference_name)
        elif reference_vcf_key not in vcf_reader.metadata:
            raise ValueError("Unable to infer reference genome for %s" % (
                path,))
        else:
            reference_path = vcf_reader.metadata[reference_vcf_key]
            reference_name = infer_reference_name(reference_path)
        ensembl_version = ensembl_release_number_for_reference_name(
            reference_name)

    ensembl = cached_release(ensembl_version)

    variants = []
    for record in vcf_reader:
        if not only_passing or not record.FILTER or record.FILTER == "PASS":
            for alt in record.ALT:
                # We ignore "no-call" variants, i.e. those where X.ALT = [None]
                if not alt:
                    continue
                variant = Variant(
                    contig=record.CHROM,
                    start=record.POS,
                    ref=record.REF,
                    alt=alt.sequence,
                    info=record.INFO,
                    ensembl=ensembl)
                variants.append(variant)
    return VariantCollection(
        variants=variants,
        path=path)
