# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Test simple properties of Variant objects, such as their trimming
of shared prefix/suffix strings from ref/alt fields.
"""

import pickle
from pyensembl import ensembl_grch38

from varcode import Variant
from .common import eq_

def test_insertion_shared_prefix():
    variant = Variant(1, start=10, ref="AA", alt="AAT")
    eq_(variant.contig, "1")
    eq_(variant.original_ref, "AA")
    eq_(variant.original_alt, "AAT")
    eq_(variant.original_start, 10)
    # since this variant is just an insertion of a "T", get rid of
    # the prefix context
    eq_(variant.ref, "")
    eq_(variant.alt, "T")
    # the [start,end] interval for an insertion is just the base we're
    # inserting after, which in this case is the 11th position
    eq_(variant.start, 11)
    eq_(variant.end, 11)
    eq_(variant.short_description, "chr1 g.11_12insT")
    assert variant.is_indel
    assert variant.is_insertion
    assert not variant.is_deletion

def test_insertion_no_prefix():
    variant = Variant(1, start=11, ref="", alt="T")
    eq_(variant.contig, "1")
    eq_(variant.original_ref, "")
    eq_(variant.original_alt, "T")
    eq_(variant.original_start, 11)
    eq_(variant.ref, "")
    eq_(variant.alt, "T")
    eq_(variant.start, 11)
    eq_(variant.end, 11)
    eq_(variant.short_description, "chr1 g.11_12insT")
    assert variant.is_indel
    assert variant.is_insertion
    assert not variant.is_deletion

def test_substitution_no_prefix():
    variant = Variant(1, start=11, ref="A", alt="T")
    eq_(variant.contig, "1")
    eq_(variant.original_ref, "A")
    eq_(variant.original_alt, "T")
    eq_(variant.original_start, 11)
    eq_(variant.ref, "A")
    eq_(variant.alt, "T")
    eq_(variant.start, 11)
    eq_(variant.end, 11)
    eq_(variant.short_description, "chr1 g.11A>T")
    assert not variant.is_indel
    assert not variant.is_insertion
    assert not variant.is_deletion

def test_substitution_shared_prefix():
    variant = Variant(1, start=10, ref="AA", alt="AT")
    eq_(variant.contig, "1")
    eq_(variant.original_ref, "AA")
    eq_(variant.original_alt, "AT")
    eq_(variant.original_start, 10)
    eq_(variant.ref, "A")
    eq_(variant.alt, "T")
    eq_(variant.start, 11)
    eq_(variant.end, 11)
    eq_(variant.short_description, "chr1 g.11A>T")
    assert not variant.is_indel
    assert not variant.is_insertion
    assert not variant.is_deletion

def test_deletion_shared_suffix():
    variant = Variant(1, start=10, ref="AAC", alt="C")
    eq_(variant.contig, "1")
    eq_(variant.original_ref, "AAC")
    eq_(variant.original_alt, "C")
    eq_(variant.original_start, 10)
    eq_(variant.ref, "AA")
    eq_(variant.alt, "")
    eq_(variant.start, 10)
    eq_(variant.end, 11)
    eq_(variant.short_description, "chr1 g.10_11delAA")
    assert variant.is_indel
    assert not variant.is_insertion
    assert variant.is_deletion

def test_deletion_no_suffix():
    variant = Variant(1, start=10, ref="AA", alt="")
    eq_(variant.contig, "1")
    eq_(variant.original_ref, "AA")
    eq_(variant.original_alt, "")
    eq_(variant.original_start, 10)
    eq_(variant.ref, "AA")
    eq_(variant.alt, "")
    eq_(variant.start, 10)
    eq_(variant.end, 11)
    eq_(variant.short_description, "chr1 g.10_11delAA")
    assert variant.is_indel
    assert not variant.is_insertion
    assert variant.is_deletion

def test_serialization():
    variants = [
        Variant(
            1, start=10, ref="AA", alt="AAT", genome=ensembl_grch38),
        Variant(10, start=15, ref="A", alt="G"),
        Variant(20, start=150, ref="", alt="G"),
    ]
    for original in variants:
        # This causes the variant's ensembl object to make a SQL connection,
        # which makes the ensembl object non-serializable. By calling this
        # method, we are checking that we don't attempt to directly serialize
        # the ensembl object.
        original.effects()

        # Test pickling.
        serialized = pickle.dumps(original)
        reconstituted = pickle.loads(serialized)
        eq_(original, reconstituted)

        eq_(original.contig, reconstituted.contig)
        eq_(original.ref, reconstituted.ref)
        eq_(original.alt, reconstituted.alt)
        eq_(original.start, reconstituted.start)
        eq_(original.end, reconstituted.end)
        eq_(original.original_ref, reconstituted.original_ref)
        eq_(original.original_alt, reconstituted.original_alt)
        eq_(original.original_start, reconstituted.original_start)

        # Test json.
        serialized = original.to_json()
        reconstituted = Variant.from_json(serialized)
        eq_(original, reconstituted)

def test_deserialization_old_keywords():
    old_variant_representation_json = """
    {
        "ref": "T",
        "contig": "22",
        "start": 23230319,
        "__class__": {
            "__name__": "Variant",
            "__module__": "varcode.variant"
        },
        "normalize_contig_name": true,
        "alt": "G",
        "allow_extended_nucleotides": false,
        "ensembl": {
            "__class__": {
                "__name__": "EnsemblRelease",
                "__module__": "pyensembl.ensembl_release"
            },
            "release": 75,
            "server": "ftp://ftp.ensembl.org",
            "species": {
                "__class__": {
                    "__name__": "Species",
                    "__module__": "pyensembl.species"
                },
                "latin_name": "homo_sapiens"
            }
        }
    }
    """
    variant = Variant.from_json(old_variant_representation_json)
    eq_(variant.contig, "22")
    eq_(variant.ref, "T")
    eq_(variant.alt, "G")
    eq_(variant.reference_name, "GRCh37")
    eq_(variant.normalize_contig_names, True)
    eq_(variant.allow_extended_nucleotides, False)

def test_hg19_chromosome_names():
    # trimming of mithochondrial name
    eq_(Variant("M", 1, "A", "G", "hg19", convert_ucsc_contig_names=True).contig, "MT")
    eq_(Variant("M", 1, "A", "G", "hg19", convert_ucsc_contig_names=False).contig, "M")

    eq_(Variant("chrM", 1, "A", "G", "hg19", convert_ucsc_contig_names=True).contig, "MT")
    eq_(Variant("chrM", 1, "A", "G", "hg19", convert_ucsc_contig_names=False).contig, "chrM")

    # uppercase
    eq_(Variant("chrm", 1, "A", "G", "hg19", convert_ucsc_contig_names=True).contig, "MT")
    eq_(Variant("chrm", 1, "A", "G", "hg19", convert_ucsc_contig_names=False).contig, "chrM")

    # trimming of 'chr' prefix from hg19
    eq_(Variant("chr1", 1, "A", "G", "hg19", convert_ucsc_contig_names=True).contig, "1")
    eq_(Variant("chr1", 1, "A", "G", "hg19", convert_ucsc_contig_names=False).contig, "chr1")

def test_contig_name_normalization():
    eq_(Variant(1, 1, "A", "G", normalize_contig_names=True).contig, "1")
    eq_(Variant(1, 1, "A", "G", normalize_contig_names=False).contig, 1)

    # uppercase
    eq_(Variant(
        "chrm", 1, "A", "G", normalize_contig_names=True, convert_ucsc_contig_names=False).contig, "chrM")
    eq_(Variant(
        "chrm", 1, "A", "G", normalize_contig_names=False, convert_ucsc_contig_names=False).contig, "chrm")


def test_snv_transition_transversion():
    ref_variant = Variant(1, start=100, ref="C", alt="C")
    assert not ref_variant.is_snv

    variant = Variant(1, start=100, ref="C", alt="T")
    assert variant.is_snv
    assert variant.is_transition
    assert not variant.is_transversion

    transversion = Variant(1, start=100, ref="C", alt="A")
    assert transversion.is_snv
    assert not transversion.is_transition
    assert transversion.is_transversion
