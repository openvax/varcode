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

