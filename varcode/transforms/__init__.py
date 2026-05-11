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

"""``VariantCollection`` -> ``VariantCollection`` transforms.

A transform is a pure function: given a ``VariantCollection`` and
optional auxiliary inputs (reference, phase resolver, ...), it returns
a new ``VariantCollection``. Cardinality may be preserved, reduced,
or increased; the docstring of each transform states which.

Every variant in the output carries a ``source_variants`` attribute
recording its provenance:

* ``()`` — pass-through; the variant was not derived from any other.
* ``(v,)`` — derived from a single source variant
  (e.g. left-aligned position).
* ``(v1, v2, ...)`` — combined from multiple source variants
  (e.g. paired breakends merged into one combined SV row).

Transforms compose by application. The intermediate indices a
transform builds (MATEID -> variant, phase set -> variants, ...) live
inside the function and are dropped on return.

Shipped transforms:

* :func:`pair_breakends` — merges MATEID-paired BND rows into a
  single combined ``StructuralVariant``. (reduces)

Planned (separate PRs):

* ``combine_cis_snvs`` — combine adjacent in-codon SNVs into MNVs
  given phase. (reduces)
* ``left_align_indels`` — canonical-position indel normalization.
  (preserves)
* ``split_multiallelic`` — split a multi-ALT row into one variant
  per ALT, making today's implicit parse-time behavior explicit.
  (increases)
* ``decompose_mnvs`` — inverse of ``combine_cis_snvs``: split MNVs
  into constituent SNVs for cross-tool comparison. (increases)
"""


def pair_breakends(vc):
    """Merge MATEID-paired BND rows into a single combined
    ``StructuralVariant`` per rearrangement event. (reduces)

    For each pair of :class:`~varcode.StructuralVariant` rows where
    row A's ``MATEID`` references row B's VCF ID (and vice versa),
    emit one combined row carrying both endpoints. The combined
    variant's ``source_variants`` attribute holds the two originals.

    Non-BND variants, single-row TRA, and single-ended BNDs (no
    ``MATEID``) pass through unchanged with ``source_variants=()``.

    Pairing rules:

    * Primary key: ``MATEID`` field on each variant's ``info`` dict
      against the VCF row ID stored in the collection's source
      metadata.
    * Alias: ``PARID`` (used by older GRIDSS) is treated as ``MATEID``.
    * Symmetric: row A's ``MATEID`` must equal row B's ID and row B's
      ``MATEID`` must equal row A's ID. Asymmetric references are
      logged and left unpaired.
    * If a ``MATEID`` points to an ID not present in this collection
      (filtered out, chunked load), a warning is emitted and the
      variant passes through unpaired.
    * If three or more rows share a ``MATEID`` group, the whole group
      is left unpaired with a warning (caller bug; pairing is
      ambiguous).

    Metadata merge for paired rows:

    * Genotype: both halves must agree on the per-sample ``GT``. The
      combined variant inherits the shared genotype. Disagreement
      raises (caller bug).
    * ``alt_assembly``: if exactly one half carries an assembled
      sequence, the combined variant inherits it. If both differ,
      the lexicographically earlier source row wins (deterministic
      across runs).
    * Other INFO fields: preserved from the lexicographically earlier
      source row; the other half's INFO is available via
      ``combined.source_variants``.

    Parameters
    ----------
    vc : VariantCollection
        Input collection. May contain a mix of structural and
        non-structural variants.

    Returns
    -------
    VariantCollection
        A new collection. SV rows that were halves of paired BNDs
        are replaced by combined rows; everything else (including
        SNVs/indels/MNVs and unpaired SVs) passes through.
    """
    raise NotImplementedError(
        "varcode.transforms.pair_breakends — implementation in progress, "
        "tracked in openvax/varcode#364")
