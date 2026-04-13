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
Per-sample genotype representation — see openvax/varcode#267.

A ``Genotype`` captures one sample's call at one variant locus: the
alleles observed, whether the call was phased, and any supporting
FORMAT fields (AD, DP, GQ, PS) that were present in the VCF. Zygosity
is computed relative to a specific alt allele index, which matters for
multi-allelic sites where a sample may carry a different alt from the
one being queried.

This module is intentionally free of circular imports so it can be
used anywhere in varcode without pulling in the rest of the package.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple


class Zygosity(Enum):
    """Zygosity of a sample's genotype relative to a specific alt allele.

    ``ABSENT`` is distinct from ``MISSING``: ABSENT means the call
    exists but doesn't include the alt in question (e.g. the sample
    is ref-ref, or carries a *different* alt at a multi-allelic
    site). MISSING means the call itself is ``./.`` or the sample
    wasn't called.
    """
    ABSENT = "absent"
    HETEROZYGOUS = "het"
    HOMOZYGOUS = "hom"
    MISSING = "missing"


def parse_gt_string(gt_str):
    """Parse a VCF ``GT`` string into ``(alleles, phased)``.

    Handles:
    * Standard diploid: ``"0/1"``, ``"1/1"``, ``"1/2"``
    * Phased diploid: ``"0|1"``, ``"1|0"``
    * Haploid (chrY, chrM): ``"1"``, ``"0"``
    * Polyploid: ``"0/1/1"`` → alleles = (0, 1, 1)
    * Missing: ``"./."``, ``"."``, ``""``, ``None``
    * Partial missing: ``"./1"`` → alleles = (None, 1)

    Returns a tuple ``(alleles, phased)`` where ``alleles`` is a tuple
    of ``Optional[int]`` (``None`` for missing) and ``phased`` is
    ``True`` iff the string used the ``|`` delimiter.
    """
    if not gt_str or gt_str == ".":
        return ((None,), False)
    if "|" in gt_str:
        phased = True
        parts = gt_str.split("|")
    else:
        phased = False
        parts = gt_str.split("/")
    alleles = tuple(
        None if p == "." else int(p)
        for p in parts
    )
    return alleles, phased


@dataclass(frozen=True)
class Genotype:
    """One sample's genotype at one variant locus.

    The ``alleles`` tuple encodes the observed alleles using VCF GT
    semantics: ``0`` is the reference allele, ``1`` is the first ALT
    listed on the VCF row, ``2`` is the second, and so on. ``None``
    indicates a no-call on that haplotype.

    For varcode's variant-level API, note that ``Variant.alt`` is a
    *specific* alt (a multi-allelic VCF row is split into one Variant
    per alt). When querying zygosity relative to a Variant, use the
    variant's ``alt_allele_index`` from the collection's metadata and
    add 1 to get the GT-encoded index, then call
    :meth:`zygosity_for_alt` or :meth:`carries_alt`.
    """

    raw_gt: str
    alleles: Tuple[Optional[int], ...]
    phased: bool = False
    phase_set: Optional[int] = None
    allele_depths: Optional[Tuple[int, ...]] = None
    total_depth: Optional[int] = None
    genotype_quality: Optional[int] = None

    # ---- construction ------------------------------------------------

    @classmethod
    def from_sample_info(cls, sample_info):
        """Build a Genotype from pyvcf's ``call.data._asdict()`` output.

        Handles the keys varcode normally sees: ``GT``, ``AD``, ``DP``,
        ``GQ``, ``PS``. Missing keys default to ``None``.
        """
        if sample_info is None:
            return cls(raw_gt="./.", alleles=(None, None), phased=False)
        gt_str = sample_info.get("GT")
        if gt_str is None:
            gt_str = "./."
        alleles, phased = parse_gt_string(gt_str)
        ad = sample_info.get("AD")
        return cls(
            raw_gt=gt_str,
            alleles=alleles,
            phased=phased,
            phase_set=sample_info.get("PS"),
            allele_depths=tuple(ad) if ad is not None else None,
            total_depth=sample_info.get("DP"),
            genotype_quality=sample_info.get("GQ"),
        )

    # ---- general predicates (alt-agnostic) --------------------------

    @property
    def is_called(self) -> bool:
        """True if at least one allele is non-None."""
        return any(a is not None for a in self.alleles)

    @property
    def is_missing(self) -> bool:
        return not self.is_called

    @property
    def ploidy(self) -> int:
        """Number of alleles in the call (including missing)."""
        return len(self.alleles)

    @property
    def is_haploid(self) -> bool:
        return self.ploidy == 1

    # ---- alt-relative predicates ------------------------------------

    def carries_alt(self, alt_index: int) -> bool:
        """True if this sample's genotype contains the given alt.

        ``alt_index`` uses VCF GT encoding: ``1`` is the first alt on
        the row, ``2`` is the second, etc. (i.e. one more than
        ``alt_allele_index`` from the VariantCollection metadata).
        """
        return any(
            a == alt_index
            for a in self.alleles
            if a is not None
        )

    def copies_of_alt(self, alt_index: int) -> int:
        """Number of haplotypes carrying the given alt."""
        return sum(
            1 for a in self.alleles
            if a is not None and a == alt_index
        )

    def zygosity_for_alt(self, alt_index: int) -> Zygosity:
        """Classify the sample's zygosity relative to one alt allele.

        Multi-allelic aware: ``GT=1/2`` queried for alt ``1`` returns
        ``HETEROZYGOUS`` (one copy of this alt, one of a different
        alt); queried for alt ``3`` it returns ``ABSENT``.
        """
        called = [a for a in self.alleles if a is not None]
        if len(called) == 0:
            return Zygosity.MISSING
        n_copies = sum(1 for a in called if a == alt_index)
        if n_copies == 0:
            return Zygosity.ABSENT
        if n_copies == len(called):
            return Zygosity.HOMOZYGOUS
        return Zygosity.HETEROZYGOUS

    # ---- depth helpers ----------------------------------------------

    def depth_for_alt(self, alt_index: int) -> Optional[int]:
        """Per-allele read depth for a given alt, from the ``AD`` field.

        ``AD`` is indexed with ref at position 0 and alt #1 at
        position 1, etc., so ``alt_index`` should use GT encoding
        (``1`` = first alt).
        """
        if self.allele_depths is None:
            return None
        if alt_index >= len(self.allele_depths):
            return None
        return self.allele_depths[alt_index]
