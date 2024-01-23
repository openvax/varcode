# Copyright (c) 2018. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from varcode import Variant
from varcode.effects import ExonicSpliceSite, PrematureStop


def test_STAT1_stop_gain_at_exon_boundary():
    # top priority effect for this variant should be PrematureStop,
    # even though it's also ExonicSpliceSite
    stat1_variant = Variant("2", "191872291", "G", "A", "GRCh37")
    effects = stat1_variant.effects()
    print(effects)
    assert any([e.__class__ is ExonicSpliceSite for e in effects])
    top_effect = effects.top_priority_effect()
    print(top_effect)
    assert top_effect.__class__ is PrematureStop
