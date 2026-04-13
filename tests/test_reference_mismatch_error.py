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
Regression tests for https://github.com/openvax/varcode/issues/215
(and the duplicate symptom in #246).

When a variant's ref allele doesn't match the reference genome at the
variant's position — typically because the variant was called against
a different build, the ref field was populated with the patient's
germline allele, or there's strand confusion — varcode should raise a
dedicated ``ReferenceMismatchError`` with an actionable message, not
a generic ``ValueError`` or ``AssertionError``.
"""

import pytest

import varcode
from varcode import Variant, ReferenceMismatchError


def _construct_mismatching_variant():
    """Construct a variant whose ref doesn't match the GRCh38 genome.

    Uses CFTR exon 4 (chr7:117530899-117531114 on GRCh38, + strand)
    where the real genome has specific bases. We claim ref='Z'... well,
    varcode rejects unknown nucleotides, so instead we use a valid
    base that doesn't match the genome at that position.
    """
    # chr7:117531114 on GRCh38 is G (last base of CFTR exon 4). Claim
    # the variant has ref='T' (which is wrong). This will fail the
    # transcript-vs-variant ref check.
    return Variant("7", 117531114, "T", "A", "GRCh38")


def test_ref_mismatch_raises_reference_mismatch_error():
    variant = _construct_mismatching_variant()
    with pytest.raises(ReferenceMismatchError):
        variant.effects()


def test_reference_mismatch_error_is_value_error_subclass():
    # Keep the existing contract: callers that catch ValueError still
    # see this. (predict_variant_effect_on_transcript_or_failure relies
    # on this for the Failure-effect fallback path.)
    assert issubclass(ReferenceMismatchError, ValueError)


def test_reference_mismatch_error_message_is_actionable():
    variant = _construct_mismatching_variant()
    try:
        variant.effects()
    except ReferenceMismatchError as e:
        msg = str(e)
        # Names the variant so the user can find it.
        assert "117531114" in msg
        # Shows both the expected (genome) and observed (variant) bases.
        assert "'T'" in msg   # variant's claimed ref
        assert "'G'" in msg   # actual genome base
        # Suggests the most common causes.
        assert "genome build" in msg or "germline" in msg or "strand" in msg
        # Points at the escape hatch.
        assert "raise_on_error=False" in msg
    else:
        raise AssertionError("Expected ReferenceMismatchError")


def test_reference_mismatch_error_carries_structured_fields():
    variant = _construct_mismatching_variant()
    try:
        variant.effects()
    except ReferenceMismatchError as e:
        assert e.variant == variant
        assert e.transcript is not None
        assert e.expected_ref == "G"
        assert e.observed_ref == "T"
    else:
        raise AssertionError("Expected ReferenceMismatchError")


def test_ref_mismatch_with_raise_on_error_false_returns_failure():
    # When the user opts into error suppression, the mismatch should
    # collapse into a Failure effect (the existing contract).
    from varcode.effects import Failure
    variant = _construct_mismatching_variant()
    effects = variant.effects(raise_on_error=False)
    assert any(isinstance(e, Failure) for e in effects), \
        "Expected at least one Failure effect when raise_on_error=False"


def test_reference_mismatch_error_exposed_at_package_level():
    # Users should be able to catch varcode.ReferenceMismatchError
    # without importing from a submodule.
    assert varcode.ReferenceMismatchError is ReferenceMismatchError
