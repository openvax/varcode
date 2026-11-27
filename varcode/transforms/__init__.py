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


"""``VariantCollection`` -> ``VariantCollection`` transforms.

A transform is a pure function: given a ``VariantCollection`` and
optional auxiliary inputs (reference, phase resolver, ...), it returns
a new ``VariantCollection``. Cardinality may be preserved, reduced,
or increased; the docstring of each transform states which.

Every variant in the output carries a ``source_variants`` attribute
recording its provenance:

* ``()`` — pass-through; the variant was not derived from any other.
* ``(v,)`` — derived from a single source variant
  (e.g. left-aligned position).
* ``(v1, v2, ...)`` — combined from multiple source variants
  (e.g. paired breakends merged into one combined SV row).

Transforms compose by application and are organized by domain:

* :mod:`varcode.transforms.breakends` — :func:`pair_breakends` (reduces)
* :mod:`varcode.transforms.indels` — :func:`left_align_indels` (preserves)
* :mod:`varcode.transforms.cis_snvs` — :func:`combine_cis_snvs` (reduces)
"""

from .breakends import pair_breakends
from .cis_snvs import combine_cis_snvs
from .indels import left_align_indels

__all__ = ["pair_breakends", "left_align_indels", "combine_cis_snvs"]
