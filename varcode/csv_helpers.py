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

# VariantCollection historically emits a "chr" column and
# EffectCollection emits "contig" — the reader should accept either so
# CSVs are interchangeable across the two collection types. Tracked in
# openvax/varcode#274.
CONTIG_COLUMN_ALIASES = ("contig", "chr")


def resolve_contig_column(columns):
    """Return the first of :data:`CONTIG_COLUMN_ALIASES` that appears in
    ``columns``, or ``None`` if none do.
    """
    for name in CONTIG_COLUMN_ALIASES:
        if name in columns:
            return name
    return None


def _parse_major_minor(version_string):
    """Return (major, minor) integers from a semver-ish version string.

    Tolerates trailing ``.patch``, ``+build`` and ``-pre`` suffixes;
    returns ``None`` when the first two dot-separated components aren't
    both integers.
    """
    if not version_string:
        return None
    # Strip build/pre-release metadata.
    for sep in ("+", "-"):
        if sep in version_string:
            version_string = version_string.split(sep, 1)[0]
    parts = version_string.split(".")
    if len(parts) < 2:
        return None
    try:
        return (int(parts[0]), int(parts[1]))
    except ValueError:
        return None


def warn_on_version_drift(header_metadata, current_version, source_path):
    """Emit ``warnings.warn`` when the CSV's ``varcode_version`` header
    reports a major-version mismatch against ``current_version``.

    Minor and patch drift are silent (semver guarantees compatibility
    within a major line). Major drift is load-bearing because
    annotation logic can change — the effects reconstructed on read
    may differ from the ones that were written.

    Also checks ``annotator`` and ``annotator_version`` when both are
    present in the header: warns when the CSV was written with a
    different annotator than the current default, or when the
    annotator version has a major mismatch (#271 stage 3b).
    """
    import warnings

    serialized = header_metadata.get("varcode_version")
    if serialized:
        serialized_mm = _parse_major_minor(serialized)
        current_mm = _parse_major_minor(current_version)
        if serialized_mm is not None and current_mm is not None:
            if serialized_mm[0] != current_mm[0]:
                warnings.warn(
                    "CSV at %s was written by varcode %s but you are reading it "
                    "with varcode %s. Because from_csv re-runs annotation on "
                    "read, results may differ across major versions. See "
                    "openvax/varcode#275 for context." % (
                        source_path, serialized, current_version))

    serialized_annotator = header_metadata.get("annotator")
    if serialized_annotator:
        # Lazy import to avoid a load-time cycle with varcode.annotators.
        from .annotators.registry import get_default_annotator
        current_annotator = getattr(get_default_annotator(), "name", None)
        if current_annotator and serialized_annotator != current_annotator:
            warnings.warn(
                "CSV at %s was written by annotator %r but the current "
                "default annotator is %r. from_csv re-runs annotation on "
                "read, so reconstructed effects reflect %r, not the "
                "annotator that originally produced the CSV. Use "
                "varcode.use_annotator(%r) around from_csv if you need "
                "the original annotator's output." % (
                    source_path, serialized_annotator, current_annotator,
                    current_annotator, serialized_annotator))

    serialized_annotator_version = header_metadata.get("annotator_version")
    if serialized_annotator_version:
        from .annotators.registry import get_default_annotator
        current_annotator_version = getattr(
            get_default_annotator(), "version", None)
        if current_annotator_version:
            serialized_mm = _parse_major_minor(serialized_annotator_version)
            current_mm = _parse_major_minor(current_annotator_version)
            if (serialized_mm is not None and current_mm is not None
                    and serialized_mm[0] != current_mm[0]):
                warnings.warn(
                    "CSV at %s was written by annotator version %s; "
                    "current annotator version is %s (major mismatch). "
                    "Reconstructed effects may differ from the originals." % (
                        source_path,
                        serialized_annotator_version,
                        current_annotator_version))


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
