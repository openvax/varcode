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

from pathlib import Path
import sys
import time

import pytest
from pyensembl import cached_release

from varcode.util import random_variants

from .common import run_with_cpu_limit

# Pin to an installed release. genome_for_reference_name("GRCh38") (the
# random_variants default) resolves to pyensembl's latest release, which is
# not necessarily downloaded (CI installs 75/81/95 only); release 81 matches
# the rest of the suite.
ensembl_grch38 = cached_release(81)

def _time_variant_annotation(variant_collection):
    """Annotate a collection and time it.

    Returns
    -------
    tuple of float
        CPU seconds used by this process, and wall-clock seconds elapsed.
    """
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    effects = variant_collection.effects()
    cpu_t = time.process_time() - start_cpu
    wall_t = time.perf_counter() - start_wall
    assert len(effects.groupby_variant()) == len(variant_collection)
    return cpu_t, wall_t


def _check_effect_timing(
        n_variants=100,
        random_seed=0,
        n_warmup_variants=5):
    """Check the default pipeline's throughput after a reproducible warmup."""
    warmup_collection = random_variants(
        n_warmup_variants,
        random_seed=1,
        genome=ensembl_grch38)
    warmup_collection.effects()

    variant_collection = random_variants(
        n_variants,
        random_seed=random_seed,
        genome=ensembl_grch38)
    cpu_t, wall_t = _time_variant_annotation(variant_collection)
    print("CPU: %0.4f s, wall-clock: %0.4f s for %d variants" % (
        cpu_t, wall_t, n_variants))
    assert cpu_t / n_variants < 0.1, \
        "Should use less than 100ms of CPU per variant!"


@pytest.mark.no_cover
def test_effect_timing():
    """Measure ordinary execution independently of earlier suite state.

    Public-cohort tests retain large reference and annotation caches. A fresh
    interpreter keeps their allocation/GC history and coverage instrumentation
    out of this throughput measurement. Imports and warmup remain outside the
    timer; annotation still uses the default pipeline and the original budget.

    The budget is on CPU time. Warm annotation is CPU-bound, so CPU time
    tracks varcode's own cost, while wall-clock time also counts every other
    process on the machine, including this suite's xdist workers (#679). The
    guard on the whole child, imports and warmup included, is on CPU time for
    the same reason (#702).
    """
    result = run_with_cpu_limit(
        [sys.executable, "-m", "tests.test_timings"],
        cpu_seconds=120,
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True)
    print(result.stdout)
    assert result.returncode == 0, result.stdout + result.stderr


def test_cpu_limit_stops_a_spinning_child_but_not_a_waiting_one():
    with pytest.raises(AssertionError, match="more than 1 s of CPU"):
        run_with_cpu_limit([sys.executable, "-c", "while True: pass"], cpu_seconds=1)
    # Waiting, like a child starved of CPU on a loaded machine, doesn't count against the limit.
    waiting = run_with_cpu_limit([sys.executable, "-c", "import time; time.sleep(2)"], cpu_seconds=1)
    assert waiting.returncode == 0


if __name__ == "__main__":
    _check_effect_timing()
