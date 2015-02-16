
from immuno.hla_file import read_hla_file
from immuno.load_file import load_file, expand_transcripts

def test_load_tcga_paad():
    maf_filename = 'data/PAAD.maf'
    hla_filename = 'data/PAAD.hla'
    alleles = read_hla_file(hla_filename)
    transcripts_df, vcf_df, variant_report = load_file(maf_filename)
    assert len(vcf_df) > 0
    assert len(transcripts_df) > 0

def test_load_tcga_skcm():
    maf_filename = 'data/SKCM.maf'
    hla_filename = 'data/SKCM.hla'
    alleles = read_hla_file(hla_filename)
    transcripts_df, vcf_df, variant_report = load_file(maf_filename)
    assert len(vcf_df) > 0
    assert len(transcripts_df) > 0

if __name__ == '__main__':
    from dsltools import testing_helpers
    testing_helpers.run_local_tests()
