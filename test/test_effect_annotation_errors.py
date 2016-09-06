# Copyright (c) 2015-2016. Mount Sinai School of Medicine
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

"""
Any time the effect annotation of a variant is erroneous, add it to this
test suite.
"""

from varcode import Variant
from varcode.effects import (
    Silent, PrematureStop, StopLoss, Insertion, Substitution
)

def expect_effect(
        variant,
        transcript_id=None,
        effect_class=None,
        protein_sequence=None,
        **kwargs):
    if transcript_id is None:
        effects = variant.effects()
        effect = effects.top_priority_effect()
    else:
        transcript = variant.ensembl.transcript_by_id(transcript_id)
        effect = variant.effect_on_transcript(transcript)
    if effect_class is not None:
        assert effect.__class__ is effect_class, \
            "Expected effect class %s but got %s" % (
                effect_class.__name__,
                effect.__class__.__name__)
    if protein_sequence is not None:
        assert effect.mutant_protein_sequence == protein_sequence, \
            "Expected protein sequence %s but got %s" % (
                protein_sequence,
                effect.mutant_protein_sequence)
    for field, expected_value in kwargs.items():
        actual_value = getattr(effect, field)
        assert actual_value == expected_value, \
            "Expected %s='%s' but got '%s'" % (
                field, expected_value, actual_value)

def test_issue167_insertion_of_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/167
    """
    The coding sequence at chr1/99772765 has transcript ID ENSMUST00000086738
    and looks like:
    ATG GAT TCT GTA CCA AGA CTG ACC AGC ATT TTG
                          ^ position 99772782
    The protein is ENSMUSP00000083944 and looks like:
        M   D   S   V   P   R   L   T   S   I   L
    ##fileformat=VCFv4.1
    CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    chr1    99772782    .   ACTG    ATGA    5000    .   .   .
    """
    variant = Variant("chr1", 99772782, "A", "ATGA", ensembl="GRCm38")
    expect_effect(
        variant=variant,
        effect_class=PrematureStop,
        protein_sequence="MDSVPR")

def test_issue168_frameshift_creates_silent_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/168
    """
    Using genome GRCm38, over transcript ENSMUST00000086738
    The coding sequence at chr1/99772765 has transcript ID ENSMUST00000086738
    and the end looks like:
    TTC ATC TGA ACT ATT GTG TGG TCA TCT GGT CCT CTT TTT (...)
           ^ stop codon
    Synonymous FrameShift over the stop codon: TGA -> TAG AAC ...

    ##fileformat=VCFv4.1
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    chr1    100484699   .   T   TA  5000    .   .   .
    """
    variant = Variant("chr1", 100484699, "T", "TA", "GRCm38")
    expect_effect(
        variant,
        effect_class=Silent)

def test_issue169_insertion_of_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/169
    """
    The coding sequence at chr1/99772765 looks like:
    ATG GAT TCT GTA CCA AGA CTG ACC AGC ATT TTGx
                          ^ 99772782
    The protein is ENSMUSP00000083944 and looks like:
        M   D   S   V   P   R   L   T   S   I   L
    The variant causes a stop codon to be introduced at position 21 in the
    nucleotide sequence, codon number 7 (with zero based indexing),
    so the protein sequence will look like MDSVPRL.
    Synonymous Substitution + PrematureStop insertion: CTG/TTA both encode L

    # ##fileformat=VCFv4.1
    # #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    # chr1    99772782    .   ACTG    ATTATGA 5000    .   .   .
    """
    # synonymous substitution followed by stop codon
    variant = Variant("chr1", 99772782, "ACTG", "ATTATGA", "GRCm38")
    expect_effect(
        variant,
        effect_class=PrematureStop)

    # non-synonymous substitution followed by stop codon
    variant = Variant("chr1", 99772782, "ACTG", "AACCTGA", "GRCm38")
    expect_effect(
        variant,
        effect_class=PrematureStop)

def test_issue170_stop_loss_does_not_translate_into_3prime_utr():
    # Issue: https://github.com/hammerlab/varcode/issues/170
    """
    Testcase, using genome GRCm38, over transcript ENSMUST00000086738

    The coding sequence at chr1/99772765 has transcript ID ENSMUST00000086738
    and the end looks like:
    TTC ATC TGA ACT ATT GTG TGG TCA TCT GGT CCT CTT TTT TGC AGA GGT TTC CAT CTC
             ^ stop codon
    TTT TTC TTT TCT TTC TTT TAA

    This segment translates as:
    F   I   *   T   I   V   W   S   S   G   P   L   F   C   R   G   F   H   L
    F   F   F   S   F   F   *

    A TGA>TCA stop loss is annotated as:

    StopLoss(
        variant=chr1 g.100484700G>C,
        transcript_name=Cntnap5b-001,
        transcript_id=ENSMUST00000086738,
        effect_description=p.*1293S (stop-loss))
     aa_ref = "*"
     aa_alt = "S"

    aa_alt should be further extended into the UTR resulting in:
        STIVWSSGPLFCRGFHLFFFSFF

    ##fileformat=VCFv4.1
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    chr1    100484699   .   TG  TC  5000    .   .   .
    """
    variant = Variant("chr1", 100484699, "TG", "TC", "GRCm38")
    expect_effect(
        variant,
        transcript_id="ENSMUST00000086738",
        effect_class=StopLoss,
        aa_alt="STIVWSSGPLFCRGFHLFFFSFF")

def test_issue_171_insertion_before_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/171
    """
    Test using genome GRCm38, on transcript ENSMUST00000086738
    Insertion before the stop codon: TGA -> TCC TGA
    This is annotated as:
        StopLoss(
            variant=chr1 g.100484699_100484700insCCT,
            transcript_name=Cntnap5b-001,
            transcript_id=ENSMUST00000086738,
            effect_description=p.*1293S (stop-loss))
        * aa_ref ="*"
        * aa_alt = "S"
    This should just be a clean insertion.

    ##fileformat=VCFv4.1
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    chr1    100484699   .   T   TCCT    5000    .   .   .
    """
    variant = Variant("chr1", 100484699, "T", "TCCT", "GRCm38")
    expect_effect(
        variant,
        effect_class=Insertion)

def test_issue172_insertion_after_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/172
    """
    Insertion immediately after the stop codon
        TGA ACT ATT -> TAG CAG CAG ACT ATT...

    ##fileformat=VCFv4.1
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    chr1    100484701   .   A   ACAGCAG 5000    .   .   .
    """
    variant = Variant("chr1", 100484701, "A", "ACAGCAG", ensembl="GRCm38")
    expect_effect(
        variant=variant,
        effect_class=Silent)

def test_issue174_wrong_aa_ref_for_insertion_of_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/174
    """
    chr1 99772782 . A ATGA 5000 . . .
    PrematureStop (OK)
    * aa_mutation_start_offset  = 6 (OK)
    * aa_mutation_end_offset = 6 (OK)
    * aa_ref = "L" (FAIL, just insertion, ref should be empty)
    * aa_alt = ""
    """
    variant = Variant("chr1", 99772782, "A", "ATGA", "GRCm38")
    expect_effect(
        variant=variant,
        effect_class=PrematureStop,
        aa_ref="",
        aa_alt="")

def test_issue175_wrong_end_offset_for_insertion_with_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/175
    """
    chr1 99772782 . A ACCCTGA 5000 .
    # Annotated as
    PrematureStop (OK)
    * aa_mutation_start_offset  = 6 (OK)
    * aa_mutation_end_offset = 7 (FAIL, should be 6)
    * aa_ref = "L" (FAIL, just insertion, ref should be empty)
    * aa_alt = "P"
    """
    variant = Variant("chr1", 99772782, "A", "ACCCTGA", "GRCm38")
    expect_effect(
        variant=variant,
        aa_ref="",
        aa_alt="",
        aa_mutation_start_offset=6,
        aa_mutation_end_offset=6)

def test_issue176_substitution_before_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/176
    """
    ##fileformat=VCFv4.1
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    # SNP before the stop codon of transcript ENSMUST00000086738
    # ATC TGA -> ACC TGA {I -> T}
    chr1 100484697 . T C 5000 . . .
    This is annotated as a StopLoss. It's a Substitution. VEP gets it right.
    """
    variant = Variant("chr1", 100484697, "T", "C", "GRCm38")
    expect_effect(
        variant=variant,
        transcript_id="ENSMUST00000086738",
        effect_class=Substitution,
        aa_ref="I",
        aa_alt="T")