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
(default: ``"fast"``). CI uses ``--annotator=protein_diff -m parity``
for the shared point-edit domain. Partial experimental annotators are not
expected to satisfy tests for unsupported inputs (for example SVs).
"""

import pytest

import varcode
import varcode.reference as _reference
from pyensembl import cached_release as _cached_release


# Match the suite's shared GRCh38 fixtures in tests/data.py. Choosing a test
# reference must not inspect .db (which can download missing files) or depend
# on unrelated installations in the user's cache (#493). CI and CONTRIBUTING
# explicitly install release 81; tests needing another release select it.
_original_reference_resolver = _reference.get_genome_for_ensembl_reference_name


def _resolve_test_reference(reference_name):
    if reference_name == "GRCh38":
        return _cached_release(81)
    return _original_reference_resolver(reference_name)


_reference.get_genome_for_ensembl_reference_name = _resolve_test_reference


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
