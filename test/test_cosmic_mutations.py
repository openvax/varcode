from varcode import VariantAnnotator

annot = VariantAnnotator(75)

def test_COSM3747785_NBPF10_Q363L():
    # 1   145311839   COSM3747785 A>T
    # GENE=NBPF10_ENST00000369338
    # STRAND=+
    # CDS=c.1088A>T
    # AA=p.Q363L
    result = annot.describe_variant("1", 145311839, "A", "T")

def test_COSM1333672_BCL9_Q1150delQ():
    """
    test_COSM1333672_BCL9_Q1150delQ : in-frame deletion of 3 nucleotides
    """
    # 1   147095923   COSM1333672 ACAG> A
    # GENE=BCL9_ENST00000234739
    # STRAND=+
    # CDS=c.3445_3447delCAG
    # AA=p.Q1150delQ
    result = annot.describe_variant("1", 147095923, "ACAG", "T")

def test_COSM1190996_FBX011_P57insQQQ():
    """
    test_COSM1190996_FBX011_P57insQQQ : in-frame insertion of 9 nucleotides
    """
    # 2   48132713    COSM1190996 C>CTGCTGCTGC
    # GENE=FBXO11_ENST00000403359
    # STRAND=-
    # CDS=c.146_147insGCAGCAGCA
    # AA=p.Q56_P57insQQQ;CNT=1
    result = annot.describe_variant("1", 48132713, "C", "CTGCTGCTGC")

def test_COSM1732848_CCDC109B_F264fs():
    """
    test_COSM1732848_CCDC109B_F264fs : frame shift from nucleotide deletion
    """
    # 4   110605772   COSM1732848 CT>C
    # GENE=CCDC109B_ENST00000394650
    # STRAND=+;CDS=c.787delT;AA=p.F264fs*5;CNT=1
    result = annot.describe_variant("4", 110605772, "CT", "C")

def test_COSM87531_SYNE1_E4738fs():
    """
    test_COSM87531_SYNE1_E4738fs : frame shift from nucleotide insertion
    """
    # 6   152651608   COSM87531   C>CA
    # GENE=SYNE1_ENST00000265368
    # STRAND=-
    # CDS=c.14211_14212insA
    # AA=p.E4738fs*34
    result = annot.describe_variant("6", 152651608, "C", "CA")


def test_COSM3368867_SMUG1_Q133L():
    # 12  54576295    COSM3368867 T>A
    # GENE=SMUG1_ENST00000513838
    # STRAND=-
    # CDS=c.398A>T
    # AA=p.Q133L
    result = annot.describe_variant("12", 54576295, "T", "A")

def test_COSM3508871_FBRS_K224N():
    # 16  30676364    COSM3508871 A>T
    # GENE=FBRS_ENST00000356166
    # STRAND=+
    # CDS=c.1572A>T
    # AA=p.K524N
    result = annot.describe_variant("16", 30676364, "A", "T")

def test_COSM4140493_silent():
    """
    test_COSM4140493_silent : silent coding mutation in LRP3
    """
    # 19  33696447    COSM4140493 G>A
    # GENE=LRP3_ENST00000590278
    # STRAND=+
    # CDS=c.393G>A
    # AA=p.R131R;CNT=1
    result = annot.describe_variant("19", 33696447, "G", "A")

def test_COSM1616161_L1724R():
    # 21  46932218    COSM1616161 T>G
    # GENE=COL18A1_ENST00000359759
    # STRAND=+
    # CDS=c.5171T>G
    # AA=p.L1724R
    result = annot.describe_variant("21", 46932218, "T", "G")


def test_COSM3939556_silent():
    # 22  19222059    COSM3939556 G>T
    # GENE=CLTCL1_ENST00000427926
    # STRAND=-
    # CDS=c.1140C>A
    # AA=p.A380A
    result = annot.describe_variant("22", 19222059, "G", "T")

def test_COSM1651074_IL9R_D148Y():
    # X   155234091   COSM1651074 TGG>TCT
    # GENE=IL9R_ENST00000244174
    # STRAND=+
    # CDS=c.441_442GG>CT
    # AA=p.D148Y
    result = annot.describe_variant("X", 155234091, "TGG", "TCT")

def test_COSM3682816_RBMY1D_V193A():
    # Y   24030663    COSM3682816 A>G
    # GENE=RBMY1D_ENST00000382680
    # STRAND=-
    # CDS=c.578T>C
    # AA=p.V193A
    result = annot.describe_variant("Y", 24030663, "A", "G")


