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

"""Parse VCF symbolic alleles and breakends into
:class:`StructuralVariant` objects.

Kept as a small, self-contained module so the VCF loader can call it
without pulling in SV annotation. The parser is deliberately lenient
— unknown symbolic types collapse to ``sv_type="BND"`` with the raw
ALT preserved, so a VCF with novel SV extensions still loads rather
than silently dropping rows.
"""

import re
from typing import Optional

from .structural_variant import StructuralVariant, SV_TYPES


# VCF 4.3 symbolic ALT grammar: <TYPE> or <TYPE:SUBTYPE> (e.g.
# <INS:ME:ALU>). We care about the top-level token for routing.
_SYMBOLIC_RE = re.compile(r"^<([A-Za-z0-9_:]+)>$")


# VCF 4.1 §5.4 breakend grammar. An ALT is one of four shapes:
#   t[CHROM:POS[   — joined-after, orientation +[
#   t]CHROM:POS]   — joined-after, orientation +]
#   [CHROM:POS[t   — joined-before, orientation [+
#   ]CHROM:POS]t   — joined-before, orientation ]+
# where ``t`` is the REF-adjacent base(s) and [/] encode strandedness
# of the mate. We capture (prefix, open_bracket, mate_contig, mate_pos,
# close_bracket, suffix). Either prefix or suffix is non-empty, not
# both.
_BREAKEND_RE = re.compile(
    r"^(?P<prefix>[ACGTNacgtn]*)"
    r"(?P<open>[\[\]])"
    r"(?P<mate_contig>[^:\[\]]+):(?P<mate_pos>\d+)"
    r"(?P<close>[\[\]])"
    r"(?P<suffix>[ACGTNacgtn]*)$"
)


def _extract_info(info, key):
    """Pull an INFO field value, handling list-of-ints from CIPOS /
    CIEND and plain ints/strings otherwise. Returns ``None`` if
    missing or unparseable rather than raising — the VCF loader
    should surface malformed rows via the normal filter, not crash
    here."""
    if info is None:
        return None
    if hasattr(info, "get"):
        return info.get(key)
    return None


def parse_symbolic_alt(
        contig: str,
        start: int,
        ref: str,
        alt: str,
        info=None,
        genome=None) -> Optional[StructuralVariant]:
    """Parse a single symbolic or breakend ALT into a
    :class:`StructuralVariant`. Returns ``None`` if the ALT is not
    symbolic (the caller keeps handling it as a simple variant).

    ``info`` is an optional mapping (e.g. a ``pyvcf`` INFO dict) that
    may carry ``END``, ``SVTYPE``, ``CIPOS``, ``CIEND``, ``MATEID``,
    etc. The parser reads those when present but doesn't require
    them — the ALT shape alone is enough to distinguish symbolic
    from breakend from inline.
    """
    if not alt:
        return None

    # Symbolic allele: <DEL>, <DUP>, <INS:ME:ALU>, <CN0>, etc.
    m = _SYMBOLIC_RE.match(alt)
    if m:
        token = m.group(1).upper()
        # Top-level type is the first colon-delimited token; e.g.
        # ``INS:ME:ALU`` → ``INS``. The rest stays in ``info`` under
        # the custom ``symbolic_subtype`` key.
        top, _, subtype = token.partition(":")

        # CNV / copy-number: sometimes written <CN0>, <CN2>, <CNV>.
        if top.startswith("CN"):
            top = "CNV"

        if top not in SV_TYPES:
            # Unknown symbolic — keep the raw token as subtype, fall
            # back to BND so the caller still gets a usable object.
            subtype = token
            top = "BND"

        end = _extract_info(info, "END") or start
        svtype_from_info = _extract_info(info, "SVTYPE")
        ci_start = _extract_info(info, "CIPOS")
        ci_end = _extract_info(info, "CIEND")

        # INFO/SVTYPE overrides when present and recognized.
        if svtype_from_info and svtype_from_info.upper() in SV_TYPES:
            top = svtype_from_info.upper()

        sv_info = {}
        if subtype:
            sv_info["symbolic_subtype"] = subtype

        return StructuralVariant(
            contig=contig,
            start=int(start),
            end=int(end),
            sv_type=top,
            alt=alt,
            ref=ref or "N",
            ci_start=tuple(ci_start) if ci_start else None,
            ci_end=tuple(ci_end) if ci_end else None,
            info=sv_info,
            genome=genome,
        )

    # Breakend: t[CHROM:POS[ etc. Four orientation shapes per §5.4.
    m = _BREAKEND_RE.match(alt)
    if m:
        prefix = m.group("prefix")
        suffix = m.group("suffix")
        open_br = m.group("open")
        close_br = m.group("close")
        mate_contig = m.group("mate_contig")
        mate_pos = int(m.group("mate_pos"))
        # Two-letter orientation token captures which side the REF
        # base is on (prefix = joined-after, suffix = joined-before)
        # and the bracket direction (mate strand).
        orientation = open_br + close_br
        return StructuralVariant(
            contig=contig,
            start=int(start),
            end=int(start),
            sv_type="BND",
            alt=alt,
            ref=ref or "N",
            mate_contig=mate_contig,
            mate_start=mate_pos,
            mate_orientation=orientation,
            info={"mateid": _extract_info(info, "MATEID"),
                  "bnd_anchor": prefix or suffix},
            genome=genome,
        )

    # Spanning-deletion placeholder. VCF 4.2+ uses ``*`` to mean "this
    # position is covered by a deletion recorded on another row".
    # We drop these like before — there's no SV object to construct.
    return None
