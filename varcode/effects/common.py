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


import contextlib
import os
import sys
import warnings

_PACKAGE_ROOT = os.path.realpath(os.path.dirname(os.path.dirname(__file__))) + os.sep
_CONTEXTLIB_FILE = os.path.realpath(contextlib.__file__)


def warn_deprecated(message: str) -> None:
    """Emit a ``DeprecationWarning`` attributed to the first caller outside varcode."""
    stacklevel = 2
    try:
        frame = sys._getframe(1)
        level = 2
        while frame is not None:
            filename = os.path.realpath(frame.f_code.co_filename)
            if not filename.startswith(_PACKAGE_ROOT) and filename != _CONTEXTLIB_FILE:
                stacklevel = level
                break
            frame = frame.f_back
            level += 1
    except ValueError:
        pass
    warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)


def validate_deletion_notation(deletion_notation):
    if deletion_notation not in ("sequence", "range"):
        raise ValueError("deletion_notation must be 'sequence' or 'range', got %r" % (deletion_notation,))
