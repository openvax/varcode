# Copyright (c) 2014. Mount Sinai School of Medicine
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

from varcode import Variant, VariantCollection

def test_drop_duplicates():
    v1 = Variant("1", 3000, "A", "G")
    v1_copy = Variant("1", 3000, "A", "G")
    v2 = Variant("2", 10, "G", "T")
    collection_with_duplicates = VariantCollection(
        variants=[v1, v1, v1_copy, v2],
        reference_name="hg19")
    assert len(collection_with_duplicates) == 4
    collection_without_duplicates = collection_with_duplicates.drop_duplicates()
    assert len(collection_without_duplicates) == 2
