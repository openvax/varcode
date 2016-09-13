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

from __future__ import print_function, division, absolute_import
from collections import defaultdict

from functools import wraps

def groupby_field(records, field_name):
    """
    Given a list of objects, group them into a dictionary by
    the unique values of a given field name.
    """

    # create an empty list for every new key
    groups = defaultdict(list)
    for record in records:
        value = getattr(record, field_name)
        if value is not None:
            groups[value].append(record)
    return dict(groups)

def memoize(fn):
    """Simple memoization decorator for functions and methods,
    assumes that all arguments to the function can be hashed and
    compared.
    """
    memoized_values = {}

    @wraps(fn)
    def wrapped_fn(*args, **kwargs):
        key = (args, tuple(sorted(kwargs.items())))
        try:
            return memoized_values[key]
        except KeyError:
            memoized_values[key] = fn(*args, **kwargs)
            return memoized_values[key]

    return wrapped_fn
