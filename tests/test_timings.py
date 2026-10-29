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

import time

import pytest
from pyensembl import cached_release

from varcode.util import random_variants

# Pin to an installed release. genome_for_reference_name("GRCh38") (the
# random_variants default) resolves to pyensembl's latest release, which is
# not necessarily downloaded (CI installs 75/81/95 only); release 81 matches
# the rest of the suite.
ensembl_grch38 = cached_release(81)

def _time_variant_annotation(variant_collection):
    start_t = time.perf_counter()
    effects = variant_collection.effects()
    end_t = time.perf_counter()
    assert len(effects.groupby_variant()) == len(variant_collection)
    elapsed_t = end_t - start_t
    return elapsed_t


@pytest.mark.no_cover
def test_effect_timing(
        n_variants=100,
        random_seed=0,
        n_warmup_variants=5):
    """Measure default throughput without coverage's tracing overhead.

    Functional tests cover annotation separately. Keep the same runtime
    budget for the ordinary default pipeline, using a reproducible warmup.
    """
    warmup_collection = random_variants(
        n_warmup_variants,
        random_seed=1,
        ensembl=ensembl_grch38)
    warmup_collection.effects()

    variant_collection = random_variants(
        n_variants,
        random_seed=random_seed,
        ensembl=ensembl_grch38)
    elapsed_t = _time_variant_annotation(variant_collection)
    print("Elapsed: %0.4f for %d variants" % (elapsed_t, n_variants))
    assert elapsed_t / n_variants < 0.1, \
        "Should be faster than 100ms / variant!"

if __name__ == "__main__":
    test_effect_timing()
