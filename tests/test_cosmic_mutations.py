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

from pyensembl import ensembl_grch37 as ensembl
from varcode import Variant
from varcode.effects import (
    Substitution,
    Deletion,
    Insertion,
    FrameShift,
    Silent,
    ExonicSpliceSite,
)

def _get_effect(chrom, pos, dna_ref, dna_alt, transcript_id):
    variant = Variant(chrom, pos, dna_ref, dna_alt, ensembl=ensembl)
    effects = variant.effects()
    transcript_dict = effects.top_priority_effect_per_transcript_id()
    assert transcript_id in transcript_dict, \
        "Expected transcript ID %s for variant %s not found in %s" % (
            transcript_id, variant, transcript_dict)
    effect = transcript_dict[transcript_id]

    # COSMIC seems to ignore exonic splice sites
    if isinstance(effect, ExonicSpliceSite):
        return effect.alternate_effect
    else:
        return effect

def _substitution(chrom, pos, dna_ref, dna_alt, transcript_id, aa_ref, aa_alt):
    effect = _get_effect(chrom, pos, dna_ref, dna_alt, transcript_id)
    assert isinstance(effect, Substitution), \
        "Expected effect to be substitution, got %s" % (effect,)

    assert effect.aa_ref == aa_ref, \
        "Expected aa_ref='%s' : %s but got %s : %s from %s" % (
            aa_ref, type(aa_ref),
            effect.aa_ref, type(effect.aa_ref),
            effect)
    assert effect.aa_alt == aa_alt, \
        "Expected aa_alt='%s' but got %s" % (
            aa_alt, effect)

def _silent(chrom, pos, dna_ref, dna_alt, transcript_id, aa_ref):
    effect = _get_effect(chrom, pos, dna_ref, dna_alt, transcript_id)
    assert isinstance(effect, Silent), \
        "Expected effect to be silent, got %s" % (effect,)
    assert effect.aa_ref == aa_ref, "Expected aa_ref='%s', got '%s'" % (
        aa_ref, effect.aa_ref)

def _deletion(chrom, pos, dna_ref, dna_alt, transcript_id, deleted):
    effect = _get_effect(chrom, pos, dna_ref, dna_alt, transcript_id)
    assert isinstance(effect, Deletion), \
        "Expected deletion, got %s" % (effect,)
    assert effect.aa_ref == deleted, \
        "Expected deletion of '%s' but got deletion of '%s' for %s:%d%s>%s" % (
            deleted, effect.aa_ref, chrom, pos, dna_ref, dna_alt)

def _insertion(chrom, pos, dna_ref, dna_alt, transcript_id, inserted):
    effect = _get_effect(chrom, pos, dna_ref, dna_alt, transcript_id)
    assert isinstance(effect, Insertion), \
        "Expected insertion, got %s" % (effect,)
    assert effect.aa_alt == inserted, \
        "Expected insertion of '%s' but got %s for %s:%d%s>%s" % (
            inserted,
            effect.short_description(),
            chrom,
            pos,
            dna_ref,
            dna_alt)

def _frameshift(
        chrom,
        pos,
        dna_ref,
        dna_alt,
        transcript_id,
        aa_pos,
        aa_ref):
    effect = _get_effect(chrom, pos, dna_ref, dna_alt, transcript_id)
    assert isinstance(effect, FrameShift), \
        "Expected frameshift, got %s" % (effect,)
    effect_aa_pos = effect.aa_mutation_start_offset
    assert effect.aa_ref[0] == aa_ref and effect_aa_pos + 1 == aa_pos, \
        ("Expected frameshift to replace p.%d%s but instead got %s" % (
            aa_pos, aa_ref, effect))

def test_COSM3939556_silent():
    # 22  19222059    COSM3939556 G>T
    # GENE=CLTCL1_ENST00000427926
    # STRAND=-
    # CDS=c.1140C>A
    # AA=p.A380A
    _silent("22", 19222059, "G", "T", "ENST00000427926", "A")

def test_COSM3747785_NBPF10_Q363L():
    # 1   145311839   COSM3747785 A>T
    # GENE=NBPF10_ENST00000369338
    # STRAND=+
    # CDS=c.1088A>T
    # AA=p.Q363L
    _substitution("1", 145311839, "A", "T", "ENST00000369338", "Q", "L")

def test_COSM3368867_SMUG1_Q133L():
    # 12  54576295    COSM3368867 T>A
    # GENE=SMUG1_ENST00000513838
    # STRAND=-
    # CDS=c.398A>T
    # AA=p.Q133L
    _substitution("12", 54576295, "T", "A", "ENST00000513838", "Q", "L")

def test_COSM3508871_FBRS_K224N():
    # 16  30676364    COSM3508871 A>T
    # GENE=FBRS_ENST00000356166
    # STRAND=+
    # CDS=c.1572A>T
    # AA=p.K524N
    _substitution("16", 30676364, "A", "T", "ENST00000356166", "K", "N")

def test_COSM1616161_L1724R():
    # 21  46932218    COSM1616161 T>G
    # GENE=COL18A1_ENST00000359759
    # STRAND=+
    # CDS=c.5171T>G
    # AA=p.L1724R
    _substitution("21", 46932218, "T", "G", "ENST00000359759", "L", "R")

def test_COSM1651074_IL9R_D148Y():
    # X   155234091   COSM1651074 TGG>TCT
    # GENE=IL9R_ENST00000244174
    # STRAND=+
    # CDS=c.441_442GG>CT
    # AA=p.D148Y
    _substitution("X", 155234091, "TGG", "TCT", "ENST00000244174", "D", "Y")

def test_COSM3682816_RBMY1D_V193A():
    # Y   24030663    COSM3682816 A>G
    # GENE=RBMY1D_ENST00000382680
    # STRAND=-
    # CDS=c.578T>C
    # AA=p.V193A
    _substitution("Y", 24030663, "A", "G", "ENST00000382680", "V", "A")

def test_COSM1333672_BCL9_Q1150delQ():
    """
    test_COSM1333672_BCL9_Q1150delQ : in-frame deletion of 3 nucleotides
    """
    # 1   147095923   COSM1333672 ACAG> A
    # GENE=BCL9_ENST00000234739
    # STRAND=+
    # CDS=c.3445_3447delCAG
    # AA=p.Q1150delQ
    _deletion("1", 147095923, "ACAG", "A", "ENST00000234739", "Q")

def test_COSM1190996_FBX011_P57insQQQ():
    """
    test_COSM1190996_FBX011_P57insQQQ : in-frame insertion of 9 nucleotides
    """
    # 2   48132713    COSM1190996 C>CTGCTGCTGC
    # GENE=FBXO11_ENST00000403359
    # STRAND=-
    # CDS=c.146_147insGCAGCAGCA
    # AA=p.Q56_P57insQQQ;CNT=1
    _insertion("2", 48132713, "C", "CTGCTGCTGC", "ENST00000403359", "QQQ")

def test_COSM1732848_CCDC109B_F264fs():
    """
    test_COSM1732848_CCDC109B_F264fs : frame shift from nucleotide deletion
    """
    # 4   110605772   COSM1732848 CT>C
    # GENE=CCDC109B_ENST00000394650
    # STRAND=+
    # CDS=c.787delT
    # AA=p.F264fs*5
    _frameshift(
        "4", 110605772, "CT", "C", "ENST00000394650",
        aa_pos=264,
        aa_ref="F")

def test_COSM87531_SYNE1_E4738fs():
    """
    test_COSM87531_SYNE1_E4738fs : frame shift from nucleotide insertion
    """
    # The given genomic mutation is:
    #    6   152651608   COSM87531   C>CA
    # but through some painful manual checking I realized that
    # the nucleotides here are *not* the correct ones for the
    # forward strand (SYNE1 is on the negative strand) and instead
    # it should be:
    #    6   152651608   COSM87531   C>CT
    # GENE=SYNE1_ENST00000265368
    # STRAND=-
    # CDS=c.14211_14212insA
    # AA=p.E4738fs*34
    _frameshift(
        "6", 152651608, "C", "GT", "ENST00000265368",
        aa_pos=4738,
        aa_ref="E")

def test_COSM27279_CTNNB1_Q4H():
    """
    test_COSM27279_CTNNB1_Q4H : Apply Cosmic mutation COSM27279
    transcript = 'ENST00000405570'
    pos: 41265571,
    ref : A, alt : T
    amino acids = Q -> H  @ pos 4 (mutation = Q4H)
    """
    _substitution("3", 41265571, "A", "T", "ENST00000405570", "Q", "H")
