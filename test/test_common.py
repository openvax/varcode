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

from nose.tools import eq_
import varcode

def test_memoize():
    class State(object):
        def __init__(self):
            self.x = 0

        def incr(self):
            self.x += 1

    state1 = State()
    # call incr twice and expect state to increment twice
    state1.incr()
    state1.incr()
    eq_(state1.x, 2)

    state2 = State()
    memoized = varcode.common.memoize(state2.incr)
    # call twice but should only increase once
    memoized()
    memoized()
    eq_(state2.x, 1)

def test_groupby_field():
    class Record(object):
        def __init__(self, x, y):
            self.x = x
            self.y = y

        def __eq__(self, other):
            return self.x == other.x and self.y == other.y

        def __str__(self):
            return "Record(%s, %s)" % (self.x, self.y)

        def __repr__(self):
            return str(self)

    r1_2 = Record(1, 2)
    r10_20 = Record(10, 20)
    r1_3 = Record(1, 3)
    data = [r1_2, r10_20, r1_3]
    grouped_dict = varcode.common.groupby_field(data, 'x')
    eq_(tuple(sorted(grouped_dict.keys())), (1, 10))
    eq_(grouped_dict[1], [r1_2, r1_3])
    eq_(grouped_dict[10], [r10_20])
