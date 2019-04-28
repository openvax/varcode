from nose.tools import eq_
from varcode.cli.effects_script import main as run_script
from varcode import Variant
from tempfile import NamedTemporaryFile
import pandas as pd


def test_varcode_effects_script_kras_g12d_top_effect():
    """
    Load a variant collection with combines the ovarian cancer test VCF
    and a small number of variants from dbSNP
    """
    kras_g12d_variant = Variant(
        12,
        25398284,
        "C",
        "T",
        "GRCh37")
    commandline_args = ["--genome", "grch37", "--only-coding", "--one-per-variant"]
    commandline_args.append("--variant")
    commandline_args.append(str(kras_g12d_variant.contig))
    commandline_args.append(str(kras_g12d_variant.start))
    commandline_args.append(str(kras_g12d_variant.original_ref))
    commandline_args.append(str(kras_g12d_variant.original_alt))
    with NamedTemporaryFile(mode="r+", delete=True) as f:
        commandline_args.extend(["--output-csv", f.name])
        run_script(commandline_args)
        f.flush()
        df = pd.read_csv(f.name)
    eq_(len(df), 1)
    eq_(df.loc[0].gene_name, "KRAS")
    eq_(df.iloc[0].effect, "p.G12D")

