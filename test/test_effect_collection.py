# Copyright (c) 2016. Mount Sinai School of Medicine
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

"""
Test properties of EffectCollection
"""

from nose.tools import eq_

from .data import tcga_ov_variants
tcga_ov_effects = tcga_ov_variants.effects()

def test_to_dataframe():
    df = tcga_ov_effects.to_dataframe()
    eq_(len(tcga_ov_effects), len(df))
