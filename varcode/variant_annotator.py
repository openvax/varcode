from variant import Variant

import pyensembl
from pyensembl.biotypes import is_coding_biotype

class VariantAnnotator(object):
    def __init__(self, ensembl_release):
        self.ensembl = pyensembl.EnsemblRelease(ensembl_release)


    def variant_gene_ids(self, contig, pos, number_modified_bases=1):
        """
        Parameters
        ----------

        contig : str
            Chromosome or contig name

        pos : int
            Position in the chromosome

        number_modified_bases : int
            How many reference bases were changed or deleted?
        """
        return self.ensembl.gene_ids_at_locus(
            contig, pos, pos + number_modified_bases)


    def variant_transcript_ids(self, contig, pos, number_modified_bases=1):
        """
        Parameters
        ----------

        contig : str
            Chromosome or contig name

        pos : int
            Position in the chromosome

        number_modified_bases : int
            How many reference bases were changed or deleted?
        """
        return self.ensembl.transcript_ids_at_locus(
            contig, pos, pos + number_modified_bases)

    def make_intergenic(self, variant):
        return Annot(
            variant=variant,
            variant_type='intergenic',
            genes=[],
            transcripts={},
            coding_effects={})


    def make_intronic(self, variant, genes, transcripts):
        return Annot(
            variant=variant,
            variant_type='intronic',
            genes=genes,
            transcripts=transcripts,
            coding_effects={})

    def overlaps_any_exon(self, transcript, contig, start, end):
        return any(
            exon.overlaps(contig=contig, start=start, end=end)
            for exon in transcript.exons)

    def group_by(self, records, field_name):
        groups = {}
        for record in records:
            value = getattr(record, field_name)
            if value in groups:
                groups[value].append(record)
            else:
                groups[value] = [record]
        return groups

    def describe_variant(self, contig, pos, ref, alt):
        variant = Variant(contig=contig, pos=pos, ref=ref, alt=alt)
        end_pos = variant.end_pos

        overlapping_genes = self.ensembl.genes_at_locus(
            contig, pos, variant.end_pos)

        if len(overlapping_genes) == 0:
            return self.make_intergenic(variant)

        overlapping_transcripts = self.ensembl.transcripts_at_locus(
                contig, pos, end_pos)

        assert len(overlapping_transcripts) > 0, \
            "No transcripts found for mutation %s:%d %s>%s" % (
                contig, pos, ref, alt)

        exonic = any(
            self.overlaps_any_exon(
                transcript, contig, start=pos, end=end_pos)
            for transcript in
            overlapping_transcripts
        )

        # group transcripts by their gene ID
        overlapping_transcript_groups = self.group_by(
            overlapping_transcripts, field_name='gene_id')

        if not exonic:
            return self.make_intronic(
                variant=variant,
                genes=overlapping_genes,
                transcripts=overlapping_transcript_groups)

        protein_variants = {}
        for transcript in overlapping_transcripts:
            if is_coding_biotype(transcript.biotype) and transcript.complete:
                seq = transcript.coding_sequence
                variant_start_offset = transcript.spliced_offset(pos)
                variant_end_offset = transcript.spliced_offset(end_pos)
                if variant_start_offset > variant_end_offset:
                    assert transcript.strand == "-"
                    original_cdna = seq[variant_end_offset:variant_start_offset+1:-1]
                    print type(original_cdna)
                    original_dna = original_cdna.complement
                else:
                    assert transcript.strand == "+"
                    original_dna = seq[variant_start_offset:variant_end_offset+1]
                assert original_dna == variant.ref, (variant, original_dna, pos, end_pos)
                original_aa = "V"
                aa_position = 600
                new_aa = "E"
                variant_string = "%s%d%s" % (original_aa, aa_position, new_aa)
                protein_variants[transcript.id] = variant_string

        if len(protein_variants) > 0:
            variant_type = "coding"
        else:
            variant_type = "coding-without-complete-transcripts"

        return Annot(
            variant=variant,
            variant_type=variant_type,
            genes=overlapping_genes,
            transcripts=overlapping_transcripts,
            coding_effects=protein_variants)
