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

"""Shared pytest fixtures for the varcode test suite.

The ``dual_annotator`` fixture parametrizes every test that uses it
over both ``"fast"`` and ``"protein_diff"`` annotators. Tests that
call ``variant.effects()`` inside a ``use_annotator(annotator_name)``
scope exercise both code paths automatically.

The ``annotator_scope`` autouse fixture sets the default annotator
for the entire test function based on the ``--annotator`` CLI option
(default: ``"fast"``). This lets CI run the full suite under
``protein_diff`` with ``pytest --annotator=protein_diff`` to catch
parity regressions across ALL tests, not just the explicit parity
harness.
"""

import os

import pytest

import varcode
import varcode.reference as _reference
from pyensembl import EnsemblRelease as _EnsemblRelease
from pyensembl import cached_release as _cached_release


# --- Cap bare-"GRCh38" resolution at the newest INSTALLED release ------------
#
# A bare reference name like "GRCh38" (e.g. ``Variant("7", ..., "GRCh38")`` or
# ``from_csv(genome="GRCh38")``) resolves, via pyensembl, to the latest release
# pyensembl KNOWS about -- currently 115 -- which is not necessarily one that's
# been downloaded. CI installs GRCh38 only up to release 95 (the
# openvax/ensembl-data mirror tops out there), so those tests would try to
# auto-fetch release 115 from Ensembl's FTP at *runtime*. That's slow and
# flaky, and is the source of intermittent, whole-job
# "GTF database needs to be created, run: pyensembl install --release 115"
# failures scattered across the suite whenever the FTP is unreachable.
#
# During the test session -- and only during the test session; the library is
# untouched -- transparently cap GRCh38 string resolution at the newest release
# whose GTF index actually exists on disk. Effect annotations are stable across
# recent GRCh38 releases (verified: the full suite passes identically under 95
# and 115), so this removes the network dependency without changing any
# expected value. GRCh38 begins at Ensembl release 76, so the walk floors there
# and never crosses into a GRCh37 release.
_GRCH38_FIRST_RELEASE = 76


def _release_index_exists(release_number):
    """True if the sqlite GTF index for a human release is already built."""
    try:
        return os.path.exists(_EnsemblRelease(release_number).db.local_db_path)
    except Exception:
        return False


_original_reference_resolver = _reference.get_genome_for_ensembl_reference_name


def _resolve_grch38_to_installed_release(reference_name):
    genome = _original_reference_resolver(reference_name)
    if getattr(genome, "reference_name", None) != "GRCh38":
        return genome
    resolved_release = getattr(genome, "release", 0)
    if _release_index_exists(resolved_release):
        return genome
    for candidate in range(resolved_release - 1, _GRCH38_FIRST_RELEASE - 1, -1):
        if _release_index_exists(candidate):
            return _cached_release(candidate)
    # Nothing installed: leave the original genome so pyensembl's usual
    # "please install release N" error still surfaces with a clear message.
    return genome


_reference.get_genome_for_ensembl_reference_name = _resolve_grch38_to_installed_release


def pytest_addoption(parser):
    parser.addoption(
        "--annotator",
        action="store",
        default=None,
        help=(
            "Run the full test suite under a specific annotator "
            "(e.g. --annotator=protein_diff). Default: no override "
            "(uses whatever each test sets, which is fast unless "
            "the test explicitly picks something else)."
        ),
    )


@pytest.fixture(autouse=True)
def annotator_scope(request):
    """When ``--annotator=<name>`` is passed on the CLI, temporarily
    set it as the default for every test. Otherwise no-op.
    """
    name = request.config.getoption("--annotator")
    if name is not None:
        with varcode.use_annotator(name):
            yield
    else:
        yield


@pytest.fixture(params=["fast", "protein_diff"])
def dual_annotator(request):
    """Parametrize a test over both annotators. Use this on tests
    that exercise ``variant.effects()`` to get automatic dual-
    annotator coverage::

        def test_something(dual_annotator):
            with varcode.use_annotator(dual_annotator):
                effects = variant.effects()
                ...
    """
    return request.param
