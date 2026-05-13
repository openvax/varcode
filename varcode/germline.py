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

"""Germline-aware annotation: input contract (#268, #289).

Effect prediction against the reference genome is the wrong baseline
for somatic-variant analysis in a real patient. The correct baseline
is the patient's germline. This module ships the *input contract* for
that — the :class:`GermlineContext` object — without yet wiring it
through annotator dispatch (which lands in subsequent slices of #268).

A :class:`GermlineContext` is the patient's germline plus enough
metadata for the effect classifier to handle it correctly:

* the variants themselves (as a :class:`~varcode.VariantCollection`);
* a :class:`Completeness` flag declaring whether the call set is a real
  germline (every position is implicitly ref/ref unless listed) or a
  sparser shape (somatic caller's "normal" column, hotspots-only,
  panel of normals) where absence does not imply ref/ref;
* the reference build, so a mismatch with the somatic VCF surfaces as
  a hard error before any nonsense annotation gets emitted.

The four canonical input shapes from #268 each have a constructor:

============================== ===========================
:meth:`GermlineContext.from_germline_vcf`         Route 1
:meth:`GermlineContext.from_multi_sample_vcf`     Route 2
:meth:`GermlineContext.from_variants`             Direct (tests, custom pipelines)
:meth:`GermlineContext.empty`                     Route 3 explicit fallback
============================== ===========================

Route 5 (``INFO=GERMLINE`` flags on a single VCF) is deferred — it's a
small adapter that lives in the loader once it has a clear consumer.

This module *only* ships the input contract and validation. Effect-side
work (window-based lookup, transcript construction with germline
applied, phase enumeration → multi-outcome packaging, LOH detection,
splice signal recomputation) lands in subsequent slices.
"""
from __future__ import annotations

import bisect
import enum
import warnings
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING,
    Any,
    Iterable,
    List,
    Mapping,
    Optional,
    Tuple,
)

if TYPE_CHECKING:
    from .variant_collection import VariantCollection


class Completeness(enum.Enum):
    """How exhaustive the germline call set is — the load-bearing
    flag that pins what *absence* of a call at a position means.

    The same data structure ("a list of germline variants") can come
    from very different pipelines, and downstream effect prediction
    cannot make the right call without knowing which:

    * If a position is **absent** from a real germline VCF emitted
      by a germline caller that examined the entire normal BAM, the
      patient is ref/ref there. Effect prediction proceeds
      reference-relative at that codon.
    * If a position is **absent** from the ``NORMAL`` column of a
      somatic-caller VCF, it likely means the somatic caller didn't
      emit a row — not that the position is ref/ref. The patient's
      germline state at that codon is unknown. The honest output is
      a possibility set including "unknown germline."
    * If a position is **absent** from a panel-of-normals filter
      list, it definitely doesn't imply ref/ref — the file only
      lists curated hotspots.

    Mis-treating "absent" as "ref/ref" silently produces wrong
    germline-aware effects on somatic variants in long stretches of
    the genome the somatic caller never touched. The flag exists so
    that mistake fails loud (or at least produces an honest
    possibility set) instead of silently corrupting clinical
    annotation.

    Values
    ------

    +-------------------+-----------------------------------+--------------------------+
    | Value             | Typical pipeline of origin        | Absence at a position    |
    +===================+===================================+==========================+
    | :attr:`COMPLETE`  | Germline caller (DeepVariant,     | ⇒ ref/ref                |
    |                   | HaplotypeCaller, Strelka2         |                          |
    |                   | germline) on the normal BAM       |                          |
    +-------------------+-----------------------------------+--------------------------+
    | :attr:`SPARSE`    | ``NORMAL`` column of a somatic    | ⇒ unknown (probably      |
    |                   | tumor-vs-normal VCF (Mutect2,     |   ref/ref but not        |
    |                   | Strelka2 somatic, VarScan2        |   queried). Honest output|
    |                   | somatic)                          |   is a possibility set.  |
    +-------------------+-----------------------------------+--------------------------+
    | :attr:`HOTSPOTS_  | Panel-of-normals filter list,     | ⇒ definitely unknown.    |
    | ONLY`             | ClinVar pathogenic list, single-  |   Strictly weaker        |
    |                   | hotspot allowlists                |   evidence than SPARSE.  |
    +-------------------+-----------------------------------+--------------------------+
    | :attr:`EMPTY`     | "I have no germline data" —       | n/a (no germline-aware   |
    |                   | explicit fallback, used so users  | logic runs; equivalent to|
    |                   | opt into reference-relative       | not passing germline= at |
    |                   | annotation rather than getting it | all)                     |
    |                   | by accident from a missing kwarg  |                          |
    +-------------------+-----------------------------------+--------------------------+

    What downstream slices do with this
    -----------------------------------

    Slice 3 of #268 wires ``germline=`` through annotator dispatch.
    When a somatic variant lands in a transcript window that has no
    germline calls, the annotator reads this flag to decide between:

    * ``COMPLETE`` → patient is ref/ref in this window; emit a
      single reference-relative effect.
    * ``SPARSE`` / ``HOTSPOTS_ONLY`` → patient's germline is
      unknown in this window; emit a possibility set including the
      reference-relative effect plus "germline-unknown" outcomes
      so the user sees the uncertainty.
    * ``EMPTY`` → no germline-aware logic; reference-relative.

    Constructors and defaults
    -------------------------

    :meth:`GermlineContext.from_germline_vcf` defaults to
    ``COMPLETE`` because that's almost always what
    a real germline VCF is.

    :meth:`GermlineContext.from_multi_sample_vcf` requires the
    caller to declare completeness explicitly (no default) — a
    multi-sample VCF could be either, and silently defaulting
    either direction is a correctness bug waiting to happen.

    :meth:`GermlineContext.empty` always sets ``EMPTY``.
    """
    COMPLETE = "complete"
    SPARSE = "sparse"
    HOTSPOTS_ONLY = "hotspots_only"
    EMPTY = "empty"


class GenomeBuildMismatchError(ValueError):
    """Raised when a germline VCF and a somatic VCF were called against
    different reference genome builds (e.g. GRCh37 vs GRCh38). Effect
    coordinates from the two VCFs cannot be meaningfully composed.

    Subclasses :class:`ValueError` so callers that already catch
    ``ValueError`` for ``ReferenceMismatchError`` continue to work.
    Set ``validate_reference=False`` on the call site if the user
    has explicitly lifted over one VCF into the other build and
    knows what they're doing.
    """

    def __init__(self, somatic_reference, germline_reference):
        self.somatic_reference = somatic_reference
        self.germline_reference = germline_reference
        super().__init__(
            "Genome build mismatch: somatic VCF uses %r, germline "
            "context uses %r. Effect coordinates from the two cannot "
            "be composed without lift-over. If you've already lifted "
            "over one VCF into the other build, pass "
            "validate_reference=False to skip this check." % (
                somatic_reference, germline_reference))


@dataclass(frozen=True)
class GermlineContext:
    """The patient's germline, packaged with completeness metadata
    and reference-build info for cross-VCF validation.

    Construct via the ``from_*`` classmethods rather than instantiating
    directly; the constructors apply the input-shape-specific
    validation each route needs.

    Attributes
    ----------
    variants
        The germline variants as a :class:`~varcode.VariantCollection`.
        Empty for :meth:`empty` contexts.
    completeness
        How to interpret absence-of-a-call (see :class:`Completeness`).
    reference_name
        The genome reference these variants were called against —
        ``"GRCh37"``, ``"GRCh38"``, ``"hg19"``, etc. ``None`` when the
        context is empty or the reference is unknown. Used by
        :meth:`validate_against` to fail fast on cross-VCF build
        mismatches.
    metadata
        Open-ended dict for caller-supplied annotations (source
        caller name, sample identifier, normalization tool, etc.).
        Not interpreted by varcode; rides along for downstream
        consumers and serialization.

    Examples
    --------
    Route 1 — full germline call set::

        ctx = GermlineContext.from_germline_vcf("normal.vcf")

    Route 2 — multi-sample VCF, extract a column. The user must
    declare completeness explicitly because absence-from-a-multi-
    sample column rarely means ref/ref::

        ctx = GermlineContext.from_multi_sample_vcf(
            "merged.vcf", sample="NORMAL", completeness=Completeness.SPARSE)

    Direct construction (tests, custom pipelines)::

        ctx = GermlineContext.from_variants(
            germline_variants, completeness=Completeness.COMPLETE,
            reference_name="GRCh38")

    Explicit empty context — opt-in to reference-relative fallback::

        ctx = GermlineContext.empty()
    """

    variants: "VariantCollection"
    completeness: Completeness = Completeness.COMPLETE
    reference_name: Optional[str] = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    # --- Constructors -------------------------------------------------------

    @classmethod
    def from_germline_vcf(
            cls,
            path: str,
            *,
            completeness: Completeness = Completeness.COMPLETE,
            metadata: Optional[Mapping[str, Any]] = None,
            **load_vcf_kwargs) -> "GermlineContext":
        """Load a full germline VCF into a context.

        ``load_vcf_kwargs`` are passed through to
        :func:`varcode.load_vcf` — for example ``genome=`` or
        ``only_passing=False``. The returned context defaults to
        ``Completeness.COMPLETE``; pass ``completeness=`` only if the
        VCF is something other than a real germline call set.
        """
        # Lazy import to avoid pulling vcf.py into the import graph
        # for callers who don't need it.
        from .vcf import load_vcf
        vc = load_vcf(path, **load_vcf_kwargs)
        if len(vc) == 0:
            warnings.warn(
                "Loaded germline VCF %r contains zero variants. This "
                "is almost always a wrong-file error; effect "
                "prediction will silently fall through to "
                "reference-relative on every somatic variant." % path)
        return cls(
            variants=vc,
            completeness=completeness,
            reference_name=cls._reference_name_of(vc),
            metadata=dict(metadata or {}),
        )

    @classmethod
    def from_multi_sample_vcf(
            cls,
            path: str,
            sample: str,
            *,
            completeness: Completeness,
            metadata: Optional[Mapping[str, Any]] = None,
            **load_vcf_kwargs) -> "GermlineContext":
        """Load a multi-sample VCF and extract one sample's calls as
        the germline.

        ``completeness`` is required (no default) — multi-sample VCFs
        from somatic callers (Mutect2's ``NORMAL`` column, e.g.) are
        almost always sparse, but pure-germline multi-sample VCFs
        (1000G, gnomAD batch genotyping) are complete. Forcing the
        caller to declare prevents subtle correctness bugs from
        treating a sparse column as if absence implied ref/ref.

        The sample is filtered post-load. If you need per-sample
        zygosity information, pass ``include_info=True`` (the default)
        and consult ``vc.metadata[variant]["sample_info"][sample]``
        downstream.
        """
        from .vcf import load_vcf
        vc = load_vcf(path, **load_vcf_kwargs)
        if sample not in cls._samples_of(vc):
            from .errors import SampleNotFoundError
            raise SampleNotFoundError(
                "Sample %r not present in %s. Available samples: %s" % (
                    sample, path, sorted(cls._samples_of(vc))))
        # Filter to variants where the named sample has a non-ref call.
        # The base VariantCollection isn't intrinsically sample-aware;
        # this is a conservative subset that the caller can refine.
        filtered = cls._variants_called_in_sample(vc, sample)
        meta = dict(metadata or {})
        meta.setdefault("source_path", path)
        meta.setdefault("sample", sample)
        return cls(
            variants=filtered,
            completeness=completeness,
            reference_name=cls._reference_name_of(vc),
            metadata=meta,
        )

    @classmethod
    def from_variants(
            cls,
            variants,
            *,
            completeness: Completeness = Completeness.COMPLETE,
            reference_name: Optional[str] = None,
            metadata: Optional[Mapping[str, Any]] = None) -> "GermlineContext":
        """Construct from an already-built :class:`VariantCollection`
        or any iterable of :class:`Variant` objects.

        Useful for tests, hand-built pipelines, and downstream tools
        that already have variants in memory and don't need to re-parse
        a VCF. ``reference_name`` should be passed explicitly when not
        carried by the variants themselves; otherwise cross-VCF
        validation will be a no-op.
        """
        from .variant_collection import VariantCollection
        if isinstance(variants, VariantCollection):
            vc = variants
        else:
            vc = VariantCollection(list(variants))
        if reference_name is None:
            reference_name = cls._reference_name_of(vc)
        return cls(
            variants=vc,
            completeness=completeness,
            reference_name=reference_name,
            metadata=dict(metadata or {}),
        )

    @classmethod
    def empty(cls) -> "GermlineContext":
        """Explicit no-germline context. Use this in pipelines where
        ``germline=`` is structurally required but the caller has no
        germline data — it documents intent better than passing
        ``None``, and downstream code can rely on the kwarg always
        being a :class:`GermlineContext`.

        Effect prediction with an empty context falls through to
        reference-relative annotation (no patient transcript
        construction), with no warnings — the caller has explicitly
        opted in to the fallback.
        """
        from .variant_collection import VariantCollection
        return cls(
            variants=VariantCollection([]),
            completeness=Completeness.EMPTY,
            reference_name=None,
            metadata={},
        )

    # --- Public API ---------------------------------------------------------

    def __bool__(self) -> bool:
        """Truthy when there's something to apply. ``EMPTY`` contexts
        are falsy so ``if germline_context:`` reads idiomatically."""
        return self.completeness is not Completeness.EMPTY and len(self.variants) > 0

    def __len__(self) -> int:
        return len(self.variants)

    def __iter__(self):
        return iter(self.variants)

    def validate_against(
            self,
            somatic,
            *,
            validate_reference: bool = True) -> None:
        """Cross-validate this context with a somatic
        :class:`VariantCollection`. Hard error on reference-build
        mismatch unless ``validate_reference=False``; warn on
        suspicious shapes (empty germline, sparse coverage with no
        overlap with somatic, etc.).

        Called automatically by :meth:`Variant.effects` /
        :meth:`VariantCollection.effects` when a context is supplied;
        callers running validation manually can do so up front to fail
        fast before annotation.
        """
        if not isinstance(self, GermlineContext):
            raise TypeError(
                "validate_against expects self to be a GermlineContext")
        somatic_ref = self._reference_name_of(somatic)
        if validate_reference and self.reference_name and somatic_ref:
            if self.reference_name != somatic_ref:
                raise GenomeBuildMismatchError(
                    somatic_reference=somatic_ref,
                    germline_reference=self.reference_name)
        if (self.completeness is not Completeness.EMPTY
                and len(self.variants) == 0):
            warnings.warn(
                "GermlineContext is non-empty by completeness flag (%s) "
                "but holds zero variants — likely a wrong-file or "
                "filter-too-aggressive error." % self.completeness.value)

    def variants_in_window(
            self,
            contig: str,
            start: int,
            end: int) -> Tuple:
        """Germline variants overlapping ``[start, end]`` on
        ``contig`` (inclusive on both ends).

        Used by the window-based lookup machinery (slice 2 of #268).
        Lazy interval index is built on first call and cached on the
        instance — subsequent calls are O(log N) per contig.

        Returns a tuple (immutable) so callers can safely cache the
        result without worrying about the underlying index mutating.
        """
        index = self._index()
        starts, variants = index.get(contig, ((), ()))
        if not starts:
            return ()
        # Variants are indexed by start. Find candidates whose
        # start <= end, then filter by their actual end span. Most
        # germline variants are SNVs / short indels so the candidate
        # window is small.
        cutoff = bisect.bisect_right(starts, end)
        result = []
        for i in range(cutoff):
            v = variants[i]
            v_end = getattr(v, "end", None) or v.start
            if v_end >= start:
                result.append(v)
        return tuple(result)

    # --- Internals ----------------------------------------------------------

    def _index(self) -> Mapping[str, Tuple[List[int], List]]:
        """Build (or return cached) a per-contig sorted index of
        germline variants by start position.

        Bisect-based lookup is plenty for typical germline sizes (a
        few million variants); an interval tree would be faster only
        for very wide windows, which we don't use.

        Cached on the instance via ``object.__setattr__`` because the
        dataclass is frozen.
        """
        cached = getattr(self, "_index_cache", None)
        if cached is not None:
            return cached
        by_contig: dict = {}
        for v in self.variants:
            by_contig.setdefault(v.contig, []).append(v)
        sorted_index = {
            contig: ([v.start for v in sorted(vs, key=lambda x: x.start)],
                     sorted(vs, key=lambda x: x.start))
            for contig, vs in by_contig.items()
        }
        object.__setattr__(self, "_index_cache", sorted_index)
        return sorted_index

    @staticmethod
    def _reference_name_of(vc) -> Optional[str]:
        """Best-effort: pull the reference name off a
        VariantCollection. Different load paths populate this
        differently; we try a few attribute names and bail to None
        if nothing matches."""
        for attr in ("reference_name", "genome_reference_name"):
            value = getattr(vc, attr, None)
            if value:
                return value
        # VariantCollection.reference_names() returns a set; if all
        # variants agree on one reference, use that.
        try:
            names = vc.reference_names()
        except Exception:
            return None
        if len(names) == 1:
            return next(iter(names))
        return None

    @staticmethod
    def _samples_of(vc) -> Iterable[str]:
        """Sample names declared on a VariantCollection — used by the
        multi-sample-VCF constructor to validate the requested
        sample exists before extracting it."""
        # VariantCollection.metadata is a {variant: {...}} dict;
        # sample names come from sample_info subdicts. Using a set
        # keeps the lookup O(unique samples) rather than O(variants).
        names: set = set()
        meta = getattr(vc, "metadata", None)
        if not meta:
            return names
        for per_variant in meta.values():
            sample_info = per_variant.get("sample_info") if isinstance(
                per_variant, Mapping) else None
            if sample_info:
                names.update(sample_info.keys())
        return names

    @staticmethod
    def _variants_called_in_sample(vc, sample):
        """Subset ``vc`` to variants where ``sample`` has a non-ref
        genotype.

        For slice 1 we keep this conservative: a variant is "called
        in the sample" if the GT field is present and not all-ref
        (``./.`` / ``0/0`` / ``0|0`` are ref-or-missing). Real
        sample-aware extraction (handling ``.|0`` mixed phase,
        somatic-vs-germline GT semantics, etc.) lives in slice 8 of
        the umbrella — this is the minimal sufficient implementation
        to ship the input contract.
        """
        from .variant_collection import VariantCollection
        keep = []
        meta = getattr(vc, "metadata", {})
        # The collection's source_to_metadata_dict shape varies; pull
        # per-variant metadata by lookup so we don't depend on a
        # particular VariantCollection internal layout.
        for v in vc:
            per_variant = meta.get(v) if isinstance(meta, Mapping) else None
            sample_info = (
                per_variant.get("sample_info") if isinstance(
                    per_variant, Mapping) else None)
            if not sample_info:
                continue
            cell = sample_info.get(sample)
            if not cell:
                continue
            gt = cell.get("GT") if isinstance(cell, Mapping) else None
            if gt and gt not in ("./.", ".|.", "0/0", "0|0", "."):
                keep.append(v)
        return VariantCollection(keep)


# =====================================================================
# Germline-aware effect prediction
#
# This is where the input contract from above flows through to actual
# effect classification. The shape from the design pivot in #268:
# germline isn't a wrapper around an effect, it's a transcript modifier
# applied before classification. ``predict_germline_aware_effect`` is
# the single entry point — ``predict_variant_effects`` calls it
# whenever a non-empty ``GermlineContext`` is supplied; otherwise it
# stays out of the way.
#
# Algorithm sketch (per-(somatic, transcript) call):
#
#   1. Look up germline variants whose coordinates overlap the somatic's
#      "window" on the transcript (default: the same codon for coding
#      variants — pluggable via ``window_fn``).
#   2. If the window is empty, the patient transcript is identical to
#      the reference at this locus. Annotate normally. (When the
#      context is SPARSE / HOTSPOTS_ONLY, flag the resulting effect with
#      ``germline_unknown=True`` so consumers see the uncertainty.)
#   3. Detect LOH — ``somatic`` shares (position, alt) with a germline
#      het call. Sets ``effect.is_loh = True`` and proceeds.
#   4. Resolve phase between somatic and each germline-in-window via
#      ``phase_resolver``. Hemizygous variants (chrX/Y/M) → automatic
#      cis. PS-tagged variants → use the resolver's answer. Unknown
#      phase → enumerate up to 2^n hypotheses (cap configurable).
#   5. For each hypothesis, build the patient haplotype the somatic
#      landed on (just the cis germline variants on that haplotype),
#      then build the post-somatic haplotype (cis germline + somatic).
#      Classify the diff between them using the protein-diff
#      classifier — same machinery the protein_diff annotator uses.
#   6. Single hypothesis → return the single classified effect.
#      Multiple hypotheses → wrap in ``PhaseCandidateSet`` (a
#      ``MultiOutcomeEffect`` subclass; not a wrapper around a
#      reference-relative effect).
#
# Out of scope for v1 (refinements that don't change this API):
#   * Splice-signal recomputation when germline already disrupts a
#     splice site varcode would otherwise classify against (#363).
#     Falls out partially because the patient transcript carries the
#     germline edits, but downgrade-by-already-broken needs targeted
#     work.
#   * Coverage-track-aware "unknown" regions (for SPARSE contexts
#     with explicit BedGraph coverage).
# Mitochondrial codon table (NCBI table 2) is selected automatically
# from the transcript's contig — see
# ``varcode.effects.codon_tables.codon_table_for_transcript``.
# =====================================================================


def apply_germline_to_transcript(transcript, germline_ctx, somatic_variant=None):
    """Apply germline variants from ``germline_ctx`` to ``transcript``,
    returning the patient's baseline :class:`MutantTranscript`.

    Lower-level entry point for callers that want the patient
    transcript directly without going through full effect prediction.
    Used internally by :func:`predict_germline_aware_effect`; exposed
    publicly for downstream tools (Isovar, Exacto) that want to
    compute a custom analysis on the patient haplotype.

    Behaviour:

    * If ``germline_ctx`` is empty, returns ``None``.
    * If ``somatic_variant`` is provided, restricts germline to the
      somatic's window (per :func:`default_germline_window`); else
      applies all germline variants overlapping any exon of the
      transcript.
    * If germline edits conflict (overlapping cDNA ranges) or land
      outside the CDS, returns ``None`` and the caller falls back.

    The returned object is the same shape that
    :func:`varcode.mutant_transcript.apply_variants_to_transcript`
    produces: a :class:`MutantTranscript` carrying the germline
    edits with ``mutant_protein_sequence`` populated when the edits
    land after the CDS start.
    """
    if not germline_ctx:
        return None
    from .mutant_transcript import apply_variants_to_transcript
    if somatic_variant is not None:
        contig, start, end = default_germline_window(
            somatic_variant, transcript)
        germline_in_window = germline_ctx.variants_in_window(
            contig, start, end)
    else:
        # Whole-transcript window: cover every exon. Useful for
        # callers building a single patient transcript independent of
        # any specific somatic — e.g., to translate the patient
        # protein for cohort-level analysis.
        germline_in_window = []
        try:
            for exon in transcript.exons:
                germline_in_window.extend(
                    germline_ctx.variants_in_window(
                        exon.contig, exon.start, exon.end))
        except Exception:
            return None
    if not germline_in_window:
        return None
    return apply_variants_to_transcript(
        list(germline_in_window), transcript)


def default_germline_window(somatic_variant, transcript) -> Tuple[str, int, int]:
    """Default window for looking up germline variants relevant to a
    somatic variant on a transcript.

    Returns ``(contig, start, end)`` covering the codon containing the
    somatic variant — three reference bases on each side of
    ``somatic_variant.start``. This is the window from #268's table
    for in-exon coding variants.

    Larger windows (splice signal region for splice-adjacent variants,
    same exon for frameshift candidates) are useful refinements but
    don't change the API. Callers that need them pass a custom
    ``window_fn`` to :func:`predict_germline_aware_effect`.

    Splice-adjacent: when the somatic is within 6bp of an exon-intron
    boundary, expand to a 12bp window centered on the boundary so
    germline edits to the donor / acceptor signal show up in the
    lookup. This catches the "germline broke the splice site" case
    without forcing the caller to wire up a separate window function.
    """
    contig = somatic_variant.contig
    pos = somatic_variant.start
    # Default: the codon containing the somatic. Three bases on either
    # side — slightly wider than strictly necessary so overlapping
    # frame-aware codon membership is conservative.
    start = pos - 3
    end = (getattr(somatic_variant, "end", None) or pos) + 3
    # Splice-adjacent expansion: if any of this transcript's exon
    # boundaries is within 6bp of the somatic, widen to capture the
    # canonical splice signal region (MAG | GURAGU and YAG | R, ~12bp).
    try:
        for exon in transcript.exons:
            for boundary in (exon.start, exon.end):
                if abs(boundary - pos) <= 6:
                    start = min(start, boundary - 6)
                    end = max(end, boundary + 6)
    except Exception:
        # Hand-built transcripts in tests may not have exons; the
        # default codon window is a fine fallback.
        pass
    return contig, max(1, start), end


def detect_loh(somatic_variant, germline_in_window) -> bool:
    """True when ``somatic_variant`` is identical at (position, alt)
    to a germline variant in the window.

    LOH is the most common "looks somatic but isn't really" case —
    the patient was germline het at this position, and the tumor lost
    the reference allele, so the variant call says "alt" in tumor and
    "het" in normal but the alt itself is the germline allele. We
    flag the resulting effect with ``is_loh=True`` so consumers can
    distinguish a true somatic mutation from a zygosity change.

    Only same-position-and-alt counts. A position where germline and
    somatic disagree on alt is a different mutation, not LOH.
    """
    for g in germline_in_window:
        if (g.contig == somatic_variant.contig
                and g.start == somatic_variant.start
                and g.ref == somatic_variant.ref
                and g.alt == somatic_variant.alt):
            return True
    return False


@dataclass(frozen=True)
class PhaseHypothesis:
    """One way the somatic variant might be phased relative to the
    germline variants in its window.

    ``cis`` are germline variants on the same haplotype as the
    somatic; ``trans`` are on the other haplotype. The somatic's
    effect at this locus is determined entirely by the cis set —
    edits in trans live on the haplotype the somatic *didn't* hit and
    don't change the codon the somatic edits.

    ``haplotype`` is a stable opaque label (``"A"``, ``"B"``, etc.)
    that callers use to align outcomes across axes — e.g. an RNA
    resolver returning ``evidence={"haplotype": "B"}`` aligns with the
    hypothesis whose ``haplotype == "B"``.

    ``phase_state`` mirrors the table in #268: ``"phased"``,
    ``"implicit"`` (hemizygous), or ``"unknown"`` (enumerated).
    """
    cis: Tuple = ()
    trans: Tuple = ()
    haplotype: str = "unknown"
    phase_state: str = "unknown"


def enumerate_phase_hypotheses(
        somatic_variant,
        germline_in_window,
        phase_resolver=None,
        max_hypotheses: int = 8) -> Tuple[PhaseHypothesis, ...]:
    """Enumerate plausible phase configurations of ``somatic_variant``
    relative to ``germline_in_window``.

    Three regimes:

    * **Hemizygous chromosome** (chrX/Y/M, male X) — single haplotype;
      all germline-in-window is implicitly cis. One hypothesis.
    * **Resolver answers for every pair** (``phase_resolver.in_cis``
      returns True/False for each ``(somatic, germline_v)``) — a
      single deterministic hypothesis with cis/trans assigned per
      the resolver. ``phase_state="phased"``.
    * **Phase unknown** — enumerate all 2^n cis/trans assignments
      across n germline variants. Cap at ``max_hypotheses``; emit a
      single ``"unknown"`` placeholder when the cap is exceeded
      (consumers see a ``TooManyHypotheses`` evidence flag).

    The cap is configurable so downstream pipelines that tolerate more
    hypotheses (long-read with rich phasing, manual analyses) can
    raise it. Default 8 = up to 3 germline variants in a window
    fully unphased.
    """
    # Hemizygous: chrX (in males), chrY, chrM. Detection is
    # heuristic since varcode doesn't carry sex info — the X
    # chromosome detection conservatively requires the genome to
    # claim hemizygosity. For v1 we treat chrM (mitochondrial) and
    # chrY as definitely hemizygous; chrX is treated as diploid
    # by default (the female case; the male case undergenerates
    # but doesn't misgenerate).
    contig = str(somatic_variant.contig).lstrip("chr").upper()
    if contig in ("M", "MT", "Y"):
        return (PhaseHypothesis(
            cis=tuple(germline_in_window),
            trans=(),
            haplotype="A",
            phase_state="implicit"),)

    # Resolver-based phase: ask the resolver for each germline
    # variant. If the resolver answers for all of them, single
    # deterministic hypothesis.
    if (phase_resolver is not None
            and hasattr(phase_resolver, "in_cis")
            and germline_in_window):
        cis_list = []
        trans_list = []
        all_answered = True
        for g in germline_in_window:
            try:
                answer = phase_resolver.in_cis(somatic_variant, g)
            except Exception:
                all_answered = False
                break
            if answer is True:
                cis_list.append(g)
            elif answer is False:
                trans_list.append(g)
            else:
                all_answered = False
                break
        if all_answered:
            return (PhaseHypothesis(
                cis=tuple(cis_list),
                trans=tuple(trans_list),
                haplotype="A",
                phase_state="phased"),)

    # Phase unknown: enumerate 2^n cis/trans assignments across the
    # n germline variants. Cap to avoid blow-up.
    n = len(germline_in_window)
    if n == 0:
        # No germline in window — single trivial hypothesis (no
        # germline edits). Caller usually short-circuits before
        # reaching this branch, but it's a safe default.
        return (PhaseHypothesis(
            cis=(),
            trans=(),
            haplotype="A",
            phase_state="phased"),)

    if 2 ** n > max_hypotheses:
        # Bail out cleanly — emit a single "too many" hypothesis with
        # all germline marked cis (the conservative case where
        # somatic effect is most strongly germline-modified).
        return (PhaseHypothesis(
            cis=tuple(germline_in_window),
            trans=(),
            haplotype="unknown",
            phase_state="too_many_hypotheses"),)

    hypotheses: List[PhaseHypothesis] = []
    germline_tuple = tuple(germline_in_window)
    for mask in range(2 ** n):
        cis = []
        trans = []
        for i, g in enumerate(germline_tuple):
            if mask & (1 << i):
                cis.append(g)
            else:
                trans.append(g)
        # Label haplotype "A" for the all-cis case, "B" for all-trans,
        # "mixed" otherwise. These are opaque tags consumers use for
        # cross-axis matching with RNA evidence.
        if not trans:
            hap = "A"
        elif not cis:
            hap = "B"
        else:
            hap = "A_mixed_%d" % mask
        hypotheses.append(PhaseHypothesis(
            cis=tuple(cis),
            trans=tuple(trans),
            haplotype=hap,
            phase_state="unknown"))
    return tuple(hypotheses)


def _classify_against_patient_baseline(
        somatic_variant, transcript, hypothesis: PhaseHypothesis):
    """Build the patient-baseline and post-somatic haplotypes for a
    single phase hypothesis, then classify the somatic effect via
    the protein-diff classifier.

    Returns a :class:`MutationEffect` instance. The effect's class
    (Missense, FrameShift, etc.) reflects the diff against the
    patient's baseline — not the reference's — when ``hypothesis.cis``
    is non-empty.

    For the all-trans case (``hypothesis.cis`` empty), the patient
    baseline is identical to the reference, so the diff degenerates
    to today's reference-relative classification. That's the
    structural back-compat: when no germline lands on the somatic's
    haplotype, this branch produces what the protein-diff annotator
    would produce on its own.
    """
    from .effects.classify import classify_from_protein_diff
    from .effects.effect_classes import (
        IncompleteTranscript,
        NoncodingTranscript,
    )
    from .mutant_transcript import apply_variants_to_transcript

    if not transcript.is_protein_coding:
        return NoncodingTranscript(somatic_variant, transcript)
    if not transcript.complete:
        return IncompleteTranscript(somatic_variant, transcript)

    cis_germline = list(hypothesis.cis)

    # Patient baseline: just the cis germline applied. Empty cis →
    # baseline is the reference transcript (no edits).
    if cis_germline:
        baseline = apply_variants_to_transcript(cis_germline, transcript)
        if baseline is None:
            # The germline edits conflict (overlapping cDNA ranges) or
            # land outside the CDS — fall back to reference-relative
            # classification of the somatic alone. The caller still
            # sees a per-hypothesis result; it just degenerates.
            from .effects import predict_variant_effect_on_transcript
            return predict_variant_effect_on_transcript(
                somatic_variant, transcript)
        baseline_protein = baseline.mutant_protein_sequence
        baseline_length_delta = baseline.total_length_delta
    else:
        # No cis germline → patient haplotype == reference at this
        # locus. Use the reference protein directly.
        baseline_protein = str(transcript.protein_sequence) if (
            getattr(transcript, "protein_sequence", None)) else None
        baseline_length_delta = 0

    # Post-somatic: cis germline + somatic on the same haplotype.
    joint_variants = cis_germline + [somatic_variant]
    joint = apply_variants_to_transcript(joint_variants, transcript)
    if joint is None or joint.mutant_protein_sequence is None:
        # Joint build failed (conflicting edits, or somatic lands
        # before CDS start) — fall back to reference-relative
        # classification of the somatic alone.
        from .effects import predict_variant_effect_on_transcript
        return predict_variant_effect_on_transcript(
            somatic_variant, transcript)

    if baseline_protein is None:
        # Edge: transcript has no reference protein (incomplete or
        # non-coding got past our gate); fall back.
        from .effects import predict_variant_effect_on_transcript
        return predict_variant_effect_on_transcript(
            somatic_variant, transcript)

    # Length delta of the somatic alone, after stripping the germline-
    # only baseline shift. classify_from_protein_diff expects the diff
    # delta, not the joint delta.
    somatic_length_delta = joint.total_length_delta - baseline_length_delta

    return classify_from_protein_diff(
        variant=somatic_variant,
        transcript=transcript,
        ref_protein=baseline_protein,
        mut_protein=joint.mutant_protein_sequence,
        length_delta=somatic_length_delta,
        mutant_transcript=joint)


def predict_germline_aware_effect(
        somatic_variant,
        transcript,
        germline_ctx: GermlineContext,
        annotator,
        phase_resolver=None,
        window_fn=default_germline_window,
        max_hypotheses: int = 8):
    """Predict the effect of ``somatic_variant`` on ``transcript``
    against the patient's germline-applied transcript.

    Single entry point for germline-aware effect prediction.
    :func:`varcode.effects.predict_variant_effects` calls this whenever
    a non-empty :class:`GermlineContext` is supplied; otherwise it
    bypasses the germline path entirely and the existing annotator
    dispatch produces today's reference-relative output unchanged.

    Behaviour by case:

    * **No germline in the somatic's window** — patient transcript ≡
      reference transcript at this locus; delegate to ``annotator``
      directly. SPARSE / HOTSPOTS_ONLY contexts mark the result with
      ``effect.germline_unknown = True`` so consumers can see the
      uncertainty.
    * **Germline in window, phase known** (resolver answers, or
      hemizygous, or all-cis-by-zygosity) — single patient haplotype;
      classify against it via :func:`_classify_against_patient_baseline`.
    * **Germline in window, phase unknown** — enumerate hypotheses
      (capped via ``max_hypotheses``), classify each, wrap in
      :class:`~varcode.effects.effect_classes.PhaseCandidateSet`.

    LOH (``somatic`` matches germline at position+alt with het zygosity)
    sets ``effect.is_loh = True`` regardless of which branch ran.

    ``window_fn`` is the pluggable window selector — defaults to
    :func:`default_germline_window` (codon-level, with splice-signal
    expansion when the somatic is splice-adjacent). Callers that
    need different windows pass their own.
    """
    from pyensembl import Transcript
    if not isinstance(transcript, Transcript):
        # Mirrors annotator entry: SVs and other non-Transcript
        # consumers don't go through germline-aware prediction.
        return annotator.annotate_on_transcript(somatic_variant, transcript)

    contig, start, end = window_fn(somatic_variant, transcript)
    germline_in_window = germline_ctx.variants_in_window(contig, start, end)
    is_loh = detect_loh(somatic_variant, germline_in_window)

    if not germline_in_window:
        effect = annotator.annotate_on_transcript(
            somatic_variant, transcript)
        if germline_ctx.completeness in (
                Completeness.SPARSE, Completeness.HOTSPOTS_ONLY):
            effect.germline_unknown = True
        if is_loh:
            effect.is_loh = True
        return effect

    hypotheses = enumerate_phase_hypotheses(
        somatic_variant,
        germline_in_window,
        phase_resolver=phase_resolver,
        max_hypotheses=max_hypotheses)

    if len(hypotheses) == 1:
        effect = _classify_against_patient_baseline(
            somatic_variant, transcript, hypotheses[0])
        # Stash the phase metadata on the effect so consumers /
        # serializers can recover it. Single-hypothesis effects don't
        # need a possibility set, but the evidence is still useful.
        effect.germline_phase_state = hypotheses[0].phase_state
        effect.germline_variants_in_window = tuple(germline_in_window)
        if is_loh:
            effect.is_loh = True
        return effect

    # Multiple hypotheses → possibility set.
    candidates = tuple(
        _classify_against_patient_baseline(
            somatic_variant, transcript, h)
        for h in hypotheses)
    from .effects.effect_classes import PhaseCandidateSet
    effect = PhaseCandidateSet(
        variant=somatic_variant,
        transcript=transcript,
        candidates=candidates,
        hypotheses=hypotheses,
        germline_variants=tuple(germline_in_window))
    if is_loh:
        effect.is_loh = True
    return effect
