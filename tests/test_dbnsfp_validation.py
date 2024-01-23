# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

from pyensembl import ensembl_grch37
from varcode import Variant
from varcode.effects import (
    ExonicSpliceSite,
    Substitution,
    TranscriptMutationEffect
)
import pandas as pd

from .data import data_path

def validate_transcript_mutation(
        ensembl_transcript_id,
        chrom,
        dna_position,
        dna_ref,
        dna_alt,
        aa_pos,
        aa_alt):
    variant = Variant(chrom, dna_position, dna_ref, dna_alt, ensembl_grch37)
    effects = variant.effects()
    transcript_id_dict = {
        effect.transcript.id: effect
        for effect in effects
        if isinstance(effect, TranscriptMutationEffect)
    }
    assert ensembl_transcript_id in transcript_id_dict, \
        "%s not found in %s" % (ensembl_transcript_id, transcript_id_dict)
    effect = transcript_id_dict[ensembl_transcript_id]

    if isinstance(effect, ExonicSpliceSite):
        # exonic splice site mutations carry with them an alternate effect
        # which is what we check against dbNSFP (since that database seemed
        # to ignore exonic splicing mutations)
        effect = effect.alternate_effect

    assert isinstance(effect, Substitution), \
        "Expected substitution (aa_pos=%d, aa_alt=%s) but got %s" % (
            aa_pos, aa_alt, effect)
    effect_aa_pos = effect.aa_mutation_start_offset
    effect_aa_alt = effect.mutant_protein_sequence[effect_aa_pos]
    assert (
        effect_aa_pos + 1 == aa_pos and
        effect_aa_alt == aa_alt), \
        "Mutant amino acid %s not found at %d for chr%s:%s %s>%s : %s" % (
            aa_alt,
            aa_pos,
            chrom,
            dna_position,
            dna_ref,
            dna_alt,
            effect)

def test_dbnsfp_validation_set():
    # check that amino acid substitution gives
    # same answer as subset of dbNSFP entries (using Ensembl 75)

    # columns for validation dataset:
    # - aa_pos : base-1 position within protein
    # - dna_alt : non-reference DNA nucleotide
    # - chrom : choromosome
    # - ensembl_transcript : transcript ID
    # - dna_position : base-1 position within chromosome
    # - dna_ref : reference DNA nucleotide

    # pylint: disable=no-member
    # pylint gets confused by read_csv
    validation_set = pd.read_csv(data_path('dbnsfp_validation_set.csv'))
    for _, row in validation_set.iterrows():
        args = (
            row['ensembl_transcript'],
            row['chrom'],
            row['dna_position'],
            row['dna_ref'],
            row['dna_alt'],
            row['aa_pos'],
            row['aa_alt']
        )
        # making this a generator so every row shows up as its
        # owns test in nose
        yield (validate_transcript_mutation,) + args

if __name__ == '__main__':
    for test_tuple in test_dbnsfp_validation_set():
        f = test_tuple[0]
        args = test_tuple[1:]
        f(*args)
