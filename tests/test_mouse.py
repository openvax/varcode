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

from .common import eq_

from varcode import load_vcf, load_vcf_fast, Variant
from varcode.effects import Substitution
from pyensembl import Genome, EnsemblRelease
from .data import data_path

MOUSE_ENSEMBL_RELEASE = 95
SERVER = "ftp://ftp.ensembl.org"
MOUSE_GTF_PATH = \
    SERVER + "/pub/release-%d/gtf/mus_musculus/Mus_musculus.GRCm38.%d.gtf.gz" % (
        MOUSE_ENSEMBL_RELEASE, MOUSE_ENSEMBL_RELEASE)
MOUSE_TRANSCRIPT_FASTA_PATH = \
    SERVER + "/pub/release-%d/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
MOUSE_PROTEIN_FASTA_PATH = \
    SERVER + "/pub/release-%d/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz" % (
        MOUSE_ENSEMBL_RELEASE)

MOUSE_VCF = data_path("mouse_vcf_dbsnp_chr1_partial.vcf")

explicit_url_genome = Genome(
    reference_name="GRCm38",
    annotation_name="ensembl",
    annotation_version=MOUSE_ENSEMBL_RELEASE,
    gtf_path_or_url=MOUSE_GTF_PATH,
    transcript_fasta_paths_or_urls=[MOUSE_TRANSCRIPT_FASTA_PATH],
    protein_fasta_paths_or_urls=[MOUSE_PROTEIN_FASTA_PATH])

ensembl_mouse_genome = EnsemblRelease(MOUSE_ENSEMBL_RELEASE, species="mouse")

def test_load_vcf_mouse_with_explicit_urls():
    variants = load_vcf(MOUSE_VCF, genome=explicit_url_genome)
    eq_(len(variants), 217)

def test_load_vcf_mouse_with_ensembl_release():
    variants = load_vcf(MOUSE_VCF, genome=ensembl_mouse_genome)
    eq_(len(variants), 217)

def test_load_vcf_mouse_with_inferred_genome():
    variants = load_vcf(MOUSE_VCF)
    eq_(len(variants), 217)

def test_specific_variant_mouse_with_explicit_urls():
    # Exon #2 at http://useast.ensembl.org/Mus_musculus/Transcript/Exons?
    # db=core;g=ENSMUSG00000017167;r=11:101170523-101190724;t=ENSMUST00000103109
    variant = Variant(
        contig=11,
        start=101177240,
        ref="G",
        alt="T",
        ensembl=explicit_url_genome)
    effects = variant.effects()
    eq_(len(effects), 2)
    substitution_effects = [
        effect
        for effect in effects
        if isinstance(effect, Substitution)
    ]
    eq_(len(substitution_effects), 1)
    substitution_effect = substitution_effects[0]
    # The coding sequence through the sub:
    # ATGATGAGTCTCCGGCTCTTCAGCATCCTGCTCGCCACG
    # GTGGTCTCTGGAGCTTGGGGCTGGGGCTACTACGGTTGC
    # (The final G is the sub: the 77th nucleotide)
    # TGC (C) -> TTC (F)
    # 78 / 3 = 26
    # 0-base = 25
    eq_(substitution_effect.mutant_protein_sequence[25], "F")
    eq_(substitution_effect.original_protein_sequence[25], "C")


def test_specific_variant_mouse_with_ensembl_genome():
    # Exon #2 at http://useast.ensembl.org/Mus_musculus/Transcript/Exons?
    # db=core;g=ENSMUSG00000017167;r=11:101170523-101190724;t=ENSMUST00000103109
    variant = Variant(
        contig=11,
        start=101177240,
        ref="G",
        alt="T",
        ensembl=ensembl_mouse_genome)
    effects = variant.effects()
    eq_(len(effects), 2)
    substitution_effects = [
        effect
        for effect in effects
        if isinstance(effect, Substitution)
    ]
    eq_(len(substitution_effects), 1)
    substitution_effect = substitution_effects[0]
    # The coding sequence through the sub:
    # ATGATGAGTCTCCGGCTCTTCAGCATCCTGCTCGCCACG
    # GTGGTCTCTGGAGCTTGGGGCTGGGGCTACTACGGTTGC
    # (The final G is the sub: the 77th nucleotide)
    # TGC (C) -> TTC (F)
    # 78 / 3 = 26
    # 0-base = 25
    eq_(substitution_effect.mutant_protein_sequence[25], "F")
    eq_(substitution_effect.original_protein_sequence[25], "C")
