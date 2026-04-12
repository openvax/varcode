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
Shared helpers for reading and writing ``# key=value`` metadata headers
at the top of CSV files produced by VariantCollection / EffectCollection.

The format is:

    # varcode_version=2.0.0
    # reference_name=GRCh38
    # ...
    col_a,col_b,col_c
    1,2,3

Any line starting with ``#`` before the first non-comment line is
treated as a metadata header of the form ``# key=value``. Lines that
don't fit that shape are ignored (but still skipped on read via
``pandas.read_csv(comment='#')``).
"""

from collections import OrderedDict


HEADER_PREFIX = "#"


def write_metadata_header(file_obj, metadata):
    """Write ``# key=value`` lines for each item in ``metadata``.

    Values that are None are skipped so the header stays tidy when a
    field wasn't populated.
    """
    for key, value in metadata.items():
        if value is None:
            continue
        file_obj.write("%s %s=%s\n" % (HEADER_PREFIX, key, value))


def read_metadata_header(path):
    """Parse leading ``#`` lines of a CSV into a dict.

    Blank / whitespace-only lines at the top are tolerated and skipped;
    parsing only terminates when a non-blank line that doesn't start
    with ``#`` is encountered (the CSV body).

    Returns
    -------
    OrderedDict
        ``key -> value`` pairs in the order they appeared. Values are
        always strings; callers are responsible for coercion.
    """
    metadata = OrderedDict()
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                # Tolerate blank lines between comment block and body.
                continue
            if not stripped.startswith(HEADER_PREFIX):
                break
            content = stripped[len(HEADER_PREFIX):].strip()
            if "=" in content:
                key, _, value = content.partition("=")
                metadata[key.strip()] = value.strip()
    return metadata
