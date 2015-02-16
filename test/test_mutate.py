# Copyright (c) 2014. Mount Sinai School of Medicine
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

from immuno import mutate

from Bio.Seq import Seq

def test_snp_mutation():
    seq = Seq("AACCTT")
    mutated = mutate.mutate(seq, 1, 'A', 'G')
    assert(mutated[0] == seq[0])
    assert(mutated[1] == 'G')
    assert(mutated[2] == seq[2])

def test_insertion_mutation():
    seq = Seq("AACT")
    mutated = mutate.mutate(seq, 1, 'A', 'AG')
    assert(len(mutated) == len(seq) + 1)
    assert(mutated[0] == seq[0])
    assert(mutated[1] == 'A')
    assert(mutated[2] == 'G')
    assert(mutated[3] == 'C')

def test_deletion_mutation():
    seq = Seq("AACT")
    mutated = mutate.mutate(seq, 1, 'ACT', 'T')
    assert(len(mutated) == 2)
    assert(mutated[0] == 'A')
    assert(mutated[1] == 'T')

def test_mutate_protein_from_transcript_snp():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")
    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("AATGCTATTCGTAGT").translate()

    region = mutate.mutate_protein_from_transcript(
        seq, 1, 'C', 'A', padding = 8)
    print(str(region))
    assert(region.seq[0] == 'N')
    assert(str(mutated_seq) == str(region.seq))
    assert(len(region.seq) == 5)
    assert(region.annot == 'T1N')

def test_mutate_protein_from_transcript_indel():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")

    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("AAATGCTATTCGTAGT").translate()
    print (str(mutated_seq))

    region = mutate.mutate_protein_from_transcript(
        seq, 1, 'C', 'AA', padding = 8)
    print(str(region))
    assert(region.seq[0] == 'K')
    assert(region.seq[1] == 'C')
    assert(region.seq[2] == 'Y')
    assert(str(region.seq) == str(mutated_seq[:-1]))
    assert(region.annot == 'T1fs')

def test_mutate_protein_prefix_stop_codon():
    seq = Seq("TAGGCTATTCGTAGT")
    prot_seq = Seq("*AIRS")

    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("TAGGATATTCGTAGT").translate()
    print (str(mutated_seq))

    region = mutate.mutate_protein_from_transcript(
        seq, 4, 'C', 'A', padding = 8)
    print(str(region))
    assert(region.seq[0] == 'D')
    assert(region.seq[1] == 'I')
    assert(str(region.seq) == str(mutated_seq[1:]))
    assert(region.annot == 'A2D')


def test_mutate_protein_from_transcript_snp_coordinates():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")
    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("AATGCTATTCGTAGT").translate()

    region = \
        mutate.mutate_protein_from_transcript(
                seq, 1, 'C', 'A', padding = 8)
    print(str(region))
    assert region.seq[0] == 'N'
    assert str(mutated_seq) == str(region.seq)
    assert len(region.seq) == 5

    assert region.mutation_start == 0
    assert region.n_removed == 1
    assert region.n_inserted == 1
    assert region.annot == 'T1N'

def test_mutate_protein_from_transcript_indel_coordinates():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")

    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("AAATGCTATTCGTAGT").translate()

    region = \
        mutate.mutate_protein_from_transcript(
            seq, 1, 'C', 'AA', padding = 8)
    print(str(region))
    assert region.seq[0] == 'K'
    assert region.seq[1] == 'C'
    assert region.seq[2] == 'Y'
    assert str(region.seq) == str(mutated_seq[:-1])

    assert region.mutation_start == 0
    assert region.n_removed == len(str(prot_seq))
    assert region.n_inserted == len(region.seq)
    assert region.annot == 'T1fs'

def test_del_annotation():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")

    assert str(seq.translate()) == str(prot_seq)

    region = \
        mutate.mutate_protein_from_transcript(
            seq, 3, "GCTATT", "", padding = 8)
    print(str(region))
    assert region.annot == 'AI2del', region

def test_stop_codon_annotation():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")
    assert str(seq.translate()) == str(prot_seq)

    stop_seq =  Seq("ACTTAGCCCATTCGTAGT")
    assert str(stop_seq.translate()) == str(Seq("T*PIRS"))

    # change the 4th-6th chars to stop codon TAG
    region = \
        mutate.mutate_protein_from_transcript(
            seq, 3, "GCT", "TAGCCC", padding = 8)
    print(str(region))
    assert region.seq == 'T', region.seq
    assert region.annot == 'A2*', region


def test_stop_codon_after_subst_annotation():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")
    assert str(seq.translate()) == str(prot_seq)

    stop_seq =  Seq("ACTCCCTAGATTCGTAGT")
    assert str(stop_seq.translate()) == str(Seq("TP*IRS"))

    # change the 4th-6th chars to stop codon TAG
    region = \
        mutate.mutate_protein_from_transcript(
            seq, 3, "GCT", "CCCTAG", padding = 8)
    print(str(region))
    assert region.n_removed == 1
    assert region.n_inserted == 1
    assert region.seq == 'TP', region.seq
    assert region.annot == 'A2P*', region


