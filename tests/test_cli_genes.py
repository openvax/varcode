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

from varcode.cli.genes_script import main as run_script
from .data import ov_wustle_variants, db_snp_variants

from tempfile import NamedTemporaryFile
import pandas as pd
import pytest


def test_varcode_effects_script():
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
        run_script(commandline_args)
        f.flush()
        combined_variants = pd.read_csv(f.name)
        assert len(combined_variants) == (len(ov_wustle_variants) + len(db_snp_variants))


def test_varcode_genes_script_missing_maf_exits_cleanly(capsys):
    with pytest.raises(SystemExit) as exc_info:
        run_script(["--maf", "/nonexistent/x.maf"])
    assert exc_info.value.code == 1
    stderr = capsys.readouterr().err
    assert "Traceback" not in stderr
    assert len(stderr.splitlines()) == 1
    assert ": error: " in stderr
    assert "/nonexistent/x.maf" in stderr


def test_varcode_command_subcommands(capsys):
    from varcode.cli.main import main
    assert main(["--version"]) == 0
    assert "varcode" in capsys.readouterr().out
    assert main(["--help"]) == 0
    assert "genes" in capsys.readouterr().out
    with pytest.raises(SystemExit) as exc_info:
        main(["genes", "--maf", "/nonexistent/x.maf"])
    assert exc_info.value.code == 1


def test_subcommand_usage_names_the_subcommand(capsys):
    # Bare `varcode <options>` means `varcode effects`, so a usage line
    # without the subcommand would send `genes` options to `effects`.
    from varcode.cli.main import main
    for name in ("effects", "genes", "check-samples"):
        with pytest.raises(SystemExit) as exc_info:
            main([name, "--help"])
        assert exc_info.value.code == 0
        assert ("usage: varcode %s " % name) in capsys.readouterr().out


def test_output_is_synonym_for_output_csv():
    from varcode.cli.genes_script import arg_parser
    assert arg_parser.parse_args(["-o", "x.csv"]).output_csv == "x.csv"
    assert arg_parser.parse_args(["--output", "x.csv"]).output_csv == "x.csv"


def test_run_subcommand_rejects_unknown_name():
    from varcode.cli.main import run_subcommand
    with pytest.raises(ValueError, match="Unknown varcode subcommand 'nope'"):
        run_subcommand("nope", [])
