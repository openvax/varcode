# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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


def apply_groupby(records, fn, skip_none=False):
    """
    Given a list of objects, group them into a dictionary by
    applying fn to each one and using returned values as a dictionary
    key.

    Parameters
    ----------
    records : list

    fn : function

    skip_none : bool
        If False, then None can be a key in the returned dictionary,
        otherwise records whose key value is None get skipped.

    Returns dict.
    """

    # create an empty list for every new key
    groups = defaultdict(list)
    for record in records:
        value = fn(record)
        if value is not None or not skip_none:
            groups[value].append(record)
    return dict(groups)


def groupby_field(records, field_name, skip_none=True):
    """
    Given a list of objects, group them into a dictionary by
    the unique values of a given field name.
    """
    return apply_groupby(
        records,
        lambda obj: getattr(obj, field_name),
        skip_none=skip_none)


def memoize(fn):
    """
    Simple memoization decorator for functions and methods,
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
