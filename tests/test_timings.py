# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
import time

from varcode.util import random_variants

def _time_variant_annotation(variant_collection):
    start_t = time.time()
    effects = variant_collection.effects()
    end_t = time.time()
    assert len(effects.groupby_variant()) == len(variant_collection)
    elapsed_t = end_t - start_t
    return elapsed_t


def test_effect_timing(
        n_variants=100,
        random_seed=0,
        n_warmup_variants=5):
    warmup_collection = random_variants(
        n_warmup_variants,
        random_seed=None)
    warmup_collection.effects()

    variant_collection = random_variants(
        n_variants,
        random_seed=random_seed)
    elapsed_t = _time_variant_annotation(variant_collection)
    print("Elapsed: %0.4f for %d variants" % (elapsed_t, n_variants))
    assert elapsed_t / n_variants < 0.1, \
        "Should be faster than 100ms / variant!"

if __name__ == "__main__":
    test_effect_timing()
