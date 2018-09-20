# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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

from collections import OrderedDict
from os.path import dirname
from .. import __file__ as package_init_file_path
from .. import __version__


def collect_version_info():
    """
    Collection the version and path of Varcode.

    TODO: add a `dependencies=False` option to also collect this info from
    major Python dependencies such as PyEnsembl
    """
    d = OrderedDict()
    d["Varcode"] = (__version__, dirname(package_init_file_path))
    return d


def print_version_info(dependencies=False):
    for (program, (version, path)) in collect_version_info().items():
        print(program)
        print("  Version: %s" % version)
        print("  Path: %s" % path)
