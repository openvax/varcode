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

"""Regression tests for #302 — `import varcode` must not drag in PyVCF3
(and transitively rpy2/R) until someone actually calls load_vcf*.

Each test spawns a fresh subprocess so sys.modules starts empty.
"""

import subprocess
import sys
import textwrap


def _run(script):
    result = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(script)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"subprocess failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    return result.stdout.strip()


def test_plain_import_does_not_load_vcf_module():
    # Bare `import varcode` must not pull in varcode.vcf (and therefore
    # must not drag PyVCF3 / rpy2 into sys.modules).
    out = _run(
        """
        import sys
        import varcode  # noqa: F401
        print("varcode.vcf" in sys.modules)
        print("vcf" in sys.modules)
        """
    )
    assert out.splitlines() == ["False", "False"]


def test_attribute_access_triggers_lazy_load():
    # Accessing varcode.load_vcf resolves via __getattr__ and pulls in
    # the vcf module on demand. Also covers `from varcode import load_vcf`.
    out = _run(
        """
        import sys
        import varcode
        assert callable(varcode.load_vcf)
        assert callable(varcode.load_vcf_fast)
        print("varcode.vcf" in sys.modules)
        """
    )
    assert out == "True"


def test_from_import_triggers_lazy_load():
    out = _run(
        """
        import sys
        from varcode import load_vcf, load_vcf_fast  # noqa: F401
        print("varcode.vcf" in sys.modules)
        """
    )
    assert out == "True"


def test_unknown_attribute_raises_attribute_error():
    out = _run(
        """
        import varcode
        try:
            varcode.definitely_not_a_real_attribute
        except AttributeError as e:
            print("ok:", str(e))
        """
    )
    assert out.startswith("ok: module 'varcode' has no attribute")
