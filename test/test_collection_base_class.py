# Copyright (c) 2015. Mount Sinai School of Medicine
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

"""
Test methods shared by both VariantCollection and EffectCollection
"""
from nose.tools import eq_, assert_not_equal
from varcode import Collection

def test_collection_len():
    collection = Collection([1, 2, 3])
    assert len(collection) == 3
    collection = Collection([])
    assert len(collection) == 0

def test_collection_eq():
    elements = ["a", "b", "c"]
    collection = Collection(elements)
    eq_(collection, collection, "collection should equal itself")
    eq_(collection, Collection(elements))


def test_collection_neq():
    elements = ["a", "b", "c"]
    c1 = Collection(elements)
    c2 = Collection(elements[:2])
    assert_not_equal(c1, c2)

def test_filter():
    collection = Collection([1, 2, 3, 4])
    filtered = collection.filter(lambda x: x < 3)
    expected = Collection([1, 2])
    eq_(filtered, expected)

def test_groupby():
    collection = Collection(["alpha", "able", "beta", "backroom"])
    groups = collection.groupby(lambda x: x[0])
    eq_(len(groups), 2)
    eq_(set(groups.keys()), {"a", "b"})

