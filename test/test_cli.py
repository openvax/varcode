from varcode.cli.variants import main as varcode_variants_main
from .data import ov_wustle_variants, db_snp_variants

from tempfile import NamedTemporaryFile
import pandas as pd

def test_varcode_variants_script():
    """
    Load a variant collection with combines the ovarian cancer test VCF
    and a small number of variants from dbSNP
    """
    commandline_args = ["--genome", "grch37"]
    commandline_args.extend(["--maf", ov_wustle_variants.path])
    for variant in db_snp_variants:
        commandline_args.append("--variant")
        commandline_args.append(str(variant.contig))
        commandline_args.append(str(variant.start))
        commandline_args.append(str(variant.original_ref))
        commandline_args.append(str(variant.original_alt))
    with NamedTemporaryFile(mode="r+", delete=True) as f:
        commandline_args.extend(["--output-csv", f.name])
        varcode_variants_main(commandline_args)
        f.flush()
        combined_variants = pd.read_csv(f.name)
        assert len(combined_variants) == (len(ov_wustle_variants) + len(db_snp_variants))
