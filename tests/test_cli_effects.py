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

from tempfile import NamedTemporaryFile
import pandas as pd
import pytest

from varcode.cli.effects_script import main as run_script
from varcode import Variant

from .common import eq_
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


@pytest.mark.parametrize(
    "commandline_args, message",
    [
        ([], "No variants loaded"),
        (["--vcf", "/nonexistent/x.vcf"], "No such file or directory"),
        (["--variant", "12", "25398284", "C", "T"], "--genome must be specified"),
        (
            ["--genome", "grch37", "--variant", "12", "abc", "C", "T"],
            "--variant position must be an integer, got 'abc'",
        ),
        (
            ["--genome", "grch37", "--variant", "12", "9" * 20, "C", "T"],
            "--variant position is out of range",
        ),
        (
            ["--genome", "grch37", "--variant", "12", "25398284", "C", "T",
             "--output-csv", "/nonexistent/dir/out.csv"],
            "output directory does not exist: /nonexistent/dir",
        ),
    ],
)
def test_varcode_effects_script_user_errors_exit_cleanly(commandline_args, message, capsys):
    with pytest.raises(SystemExit) as exc_info:
        run_script(commandline_args)
    eq_(exc_info.value.code, 1)
    stderr = capsys.readouterr().err
    assert "Traceback" not in stderr
    assert message in stderr


def test_varcode_effects_script_output_csv_directory_fails_before_annotating(
        tmp_path, capsys):
    """
    An --output-csv path which is itself a directory should be rejected by the
    same pre-flight check as a missing directory, rather than after the whole
    annotation run.
    """
    commandline_args = [
        "--genome", "grch37",
        "--variant", "12", "25398284", "C", "T",
        "--output-csv", str(tmp_path),
    ]
    with pytest.raises(SystemExit) as exc_info:
        run_script(commandline_args)
    eq_(exc_info.value.code, 1)
    stderr = capsys.readouterr().err
    assert "Traceback" not in stderr
    assert "output path is a directory: %s" % tmp_path in stderr
