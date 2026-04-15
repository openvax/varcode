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

"""SequenceDiffEffectAnnotator — classify effects from protein diff
instead of offset arithmetic (openvax/varcode#271 stage 3d, #309).

**WORK IN PROGRESS — not registered in the annotator registry,
not functional.** This module currently ships a class skeleton,
a classifier signature with a documented decision table, and
TODOs. See the PR description and #309 for the rollout plan.

See :doc:`/effect_annotation` for the user-facing guide covering
how legacy, sequence-diff, splice outcomes, and the SV roadmap
fit together. This module docstring captures implementation
detail relevant to reviewers of the classifier itself.

Algorithm
---------

1. **Gate**: non-coding / incomplete transcripts and splice-adjacent
   variants go straight to :class:`LegacyEffectAnnotator`. Sequence-
   diff owns protein-level diff classification; splice classification
   stays offset-based (#305 will rewrite that module on top of
   :class:`MutantTranscript` without sequence-diff's involvement).

2. **Fast path**: trivial single-codon SNVs dispatch to the shared
   :func:`try_fast_path_snv` helper (stage 3c). Both annotators hit
   the same code here so they can't diverge on the common case.

3. **Slow path**: materialize a :class:`MutantTranscript` via
   :func:`apply_variant_to_transcript`, translate its cDNA with the
   transcript's codon table (mt-aware via #294), and compare the
   resulting protein to the reference. The minimal edit comes from
   :func:`trim_shared_flanking_strings`; classification is a
   finite decision table (see :func:`classify_from_protein_diff`).

Generalization to structural variants (forward-looking)
-------------------------------------------------------

The current :class:`MutantTranscript` carries a single
``reference_transcript``, which expresses point variants and indels
cleanly but doesn't fit SV-class events:

* **Gene fusions** compose segments from two reference transcripts.
  See #252.
* **Translocations** can join a transcript to an intergenic region,
  producing many candidate ORFs resolvable only with RNA evidence
  (long-read transcript assembly in particular). See #257, #259.
* **Large inversions / duplications / tandem repeats** can produce
  multiple plausible mutant transcripts per DNA event.

The cross-cutting pattern is **one DNA event → possibility set of
plausible mutant proteins**, narrowed by RNA evidence (#259) — the
same shape as splice outcomes (#262/#305). The planned extensions:

1. ``MutantTranscript`` grows a ``reference_segments`` alternative
   to ``reference_transcript``, each segment a ``(transcript,
   cdna_range)``. Fusion = two segments from different transcripts;
   translocation-to-intergenic = one transcript segment plus a
   genomic-interval segment.
2. An SV-annotator returns ``List[MutantTranscript]`` (or a
   :class:`MultiOutcomeEffect` wrapper per the #299 generalization)
   when outcomes are ambiguous. RNA evidence narrows the set at
   a later stage.
3. The :func:`classify_from_protein_diff` decision table below is
   unchanged for SVs — once the mutant protein is materialized,
   classification is the same. The difference is upstream: how
   you get there.

None of this SV work is in scope for 3d (#309). Calling it out
here so the point-variant design doesn't close off the SV path.

Ownership boundary with topiary's ProteinFragment
-------------------------------------------------

Topiary's :class:`ProteinFragment` is the downstream universal
substrate for "mutant protein with source-type tag + target
intervals" (openvax/topiary#121). The split should be:

* :class:`MutantTranscript` — transcript-level reasoning. Lives in
  varcode. Owns splicing decisions, UTR readthrough, codon-table
  selection, edit provenance.
* :class:`ProteinFragment` — prediction-boundary substrate. Lives
  in topiary. Owns MHC-binding input framing, peptide-window
  slicing, source-type tagging.

The conversion from MutantTranscript to ProteinFragment is 1:1:
``fragment.sequence = mutant_transcript.mutant_protein_sequence``,
``fragment.source_type = <variant_class>``, ``fragment.fragment_id
= <transcript_id>:<variant_hash>``, and provenance fields
(``annotator_name``, evidence dict) thread through as fragment
metadata. Downstream tools should use MutantTranscript for
transcript-level reasoning and convert to ProteinFragment at the
prediction boundary — ad-hoc conversions drift across tools.

Possibility-set prior for DNA-only callers
------------------------------------------

When a single DNA event produces a possibility set of mutant
proteins (splice outcomes, fusion candidates, translocation
ORFs), and RNA evidence isn't available to narrow it, callers
still need to pick a top outcome for downstream pipelines that
consume one protein per variant (neoantigen prediction, vaccine
peptide selection). The :class:`MultiOutcomeEffect` wrapper
should carry a **prior over outcomes** even when that prior is
crude:

* Splice outcomes: the plausibility tables in
  ``varcode.splice_outcomes`` already serve this role.
* Fusions: canonical-transcript-preferred; known recurrent
  junctions (BCR-ABL, EML4-ALK, TMPRSS2-ERG) get
  possibility-set size 1.
* Translocations producing many ORFs: rank by
  (in-frame, ORF length, canonical-transcript involvement).

The API convention: ``.most_likely`` always returns the
highest-prior candidate, ``.candidates`` returns all of them
sorted prior-descending. Callers that want "best guess from DNA
alone" read ``.most_likely``; callers with RNA evidence filter
``.candidates`` by observation support. Same shape as the
existing :class:`SpliceOutcomeSet` — extending it uniformly
across variant classes is exactly what #299 tracks.

Recurrent-junction registry (for fusions)
-----------------------------------------

Known recurrent fusion junctions should short-circuit to a
single MutantTranscript rather than enumerate candidates. The
registry contract needs resolution before fusion support lands:

* Does varcode ship the registry, or does each fusion annotator
  (Exacto's plugin, Isovar's plugin) consult its own?
* If varcode ships it, what's the data source — COSMIC fusion
  entries, a curated list, user-supplied?
* If annotators consult their own, what's the common lookup
  interface so registries can coexist?

None of this needs resolution for 3d's point-variant classifier.
Flagging so the design doesn't lock in a particular answer too
early. Likely lives in its own tracking issue once #252 fusion
support starts.
"""

from ..version import __version__ as _varcode_version

# Intentional: `apply_variant_to_transcript` will be imported here in
# 3d's classifier. Keeping the import path documented in the module
# docstring only for now to avoid a ruff F401 on the WIP scaffolding.


class SequenceDiffEffectAnnotator:
    """Classify effects by diffing translated mutant protein against
    the reference protein. **WIP — skeleton only.** See module
    docstring.
    """

    name = "sequence_diff"
    version = _varcode_version
    supports = frozenset({"snv", "indel", "mnv"})
    """Variant kinds handled directly. Splice effects delegate to
    legacy; SV-class variants (``sv``, ``fusion``, ``translocation``,
    ``inversion``, ``duplication``) get added once the
    ``reference_segments`` generalization lands — see module
    docstring."""

    def __repr__(self):
        # Include the registered name in the repr so users who see
        # "<SequenceDiffEffectAnnotator ...>" in an error message can
        # immediately map it back to the annotator= kwarg string.
        return "SequenceDiffEffectAnnotator(name=%r, version=%r)" % (
            self.name, self.version)

    def annotate_on_transcript(self, variant, transcript):
        """Dispatch a single (variant, transcript) pair.

        **Not implemented.** Stage 3d will fill this in against the
        validation corpus described in #309. Raises
        :class:`NotImplementedError` so accidental registration
        catches you early.
        """
        raise NotImplementedError(
            "SequenceDiffEffectAnnotator is scaffolding only in %s. "
            "Use LegacyEffectAnnotator until #309 lands the "
            "classifier." % _varcode_version)


def classify_from_protein_diff(mutant_transcript, variant, transcript):
    """Classify a :class:`MutationEffect` from a reference/mutant
    protein pair.

    **Not implemented.** Decision table the implementation will follow
    (from #309), applied after
    :func:`trim_shared_flanking_strings` has reduced
    ``(reference_protein, mutant_protein)`` to the minimal edit
    ``(ref_delta, alt_delta)`` with shared ``prefix`` and ``suffix``:

    +----------------------------------------------+-----------------------+
    | Condition                                    | Effect                |
    +==============================================+=======================+
    | ``ref == mut``                               | Silent                |
    +----------------------------------------------+-----------------------+
    | first codon is a non-ATG alternate start     | AlternateStartCodon   |
    | and proteins otherwise match                 | (non-diff special     |
    |                                              | case, see #304)       |
    +----------------------------------------------+-----------------------+
    | ``aa_offset == 0`` and                       | StartLoss             |
    | ``not mut.startswith('M')``                  |                       |
    +----------------------------------------------+-----------------------+
    | ``ref_delta`` ends at ``len(ref)`` and       | StopLoss (requires    |
    | ``len(mut) > len(ref)``                      | 3'UTR readthrough —   |
    |                                              | extend                |
    |                                              | apply_variant_to_     |
    |                                              | transcript)           |
    +----------------------------------------------+-----------------------+
    | ``total_length_delta % 3 != 0`` and          | FrameShiftTruncation  |
    | ``len(alt_delta) == 0``                      |                       |
    +----------------------------------------------+-----------------------+
    | ``total_length_delta % 3 != 0``              | FrameShift            |
    +----------------------------------------------+-----------------------+
    | ``len(mut) < len(ref)`` with in-frame edit   | PrematureStop         |
    +----------------------------------------------+-----------------------+
    | ``len(ref_delta) == len(alt_delta) == 1``    | Substitution          |
    +----------------------------------------------+-----------------------+
    | ``len(ref_delta) == 0``                      | Insertion             |
    +----------------------------------------------+-----------------------+
    | ``len(alt_delta) == 0``                      | Deletion              |
    +----------------------------------------------+-----------------------+
    | otherwise                                    | ComplexSubstitution   |
    +----------------------------------------------+-----------------------+
    """
    raise NotImplementedError(
        "classify_from_protein_diff is WIP scaffolding — see #309.")


def _variant_is_splice_adjacent(variant, transcript):
    """True when ``variant`` sits within the canonical splice window
    (last 3 bases of an exon or first 3–6 bases of an intron) on
    ``transcript``.

    **Not implemented.** Sequence-diff delegates splice-adjacent
    variants to the legacy classifier; this gate decides when. The
    implementation can reuse the distance-to-exon helpers in
    :mod:`varcode.effects.effect_prediction`.

    Dual-dispatch decision (from vaxrank feedback on PR #310)
    ---------------------------------------------------------

    A splice-adjacent variant in the **exon body** (e.g. a missense
    at the last exonic base) is biologically *both* a splice effect
    AND a protein-level coding effect. Legacy's
    :class:`ExonicSpliceSite` carries an ``alternate_effect`` field
    for exactly this reason — it captures "what happens if splicing
    proceeds normally." Options for sequence-diff:

    **(a) Dual-dispatch**: when the variant is exon-body *and* in
    the splice window, run legacy for the splice classification
    AND return sequence-diff's protein-level classification as
    the ``alternate_effect`` on the resulting
    :class:`ExonicSpliceSite`. Output is an :class:`ExonicSpliceSite`
    whose ``alternate_effect`` comes from sequence-diff's diff
    logic instead of legacy's offset arithmetic — the two
    classifications agree on the common case because they go
    through the same shared fast path (stage 3c), and disagree
    only on the cases sequence-diff provably handles better
    (boundary codons, frameshift realignment, etc.).

    **(b) Legacy-only**: exon-body splice-adjacent variants
    delegate entirely to legacy. Protein-level consequences come
    from legacy's offset-based logic, even for the cases where
    sequence-diff's diff would be more accurate. Neoantigen
    pipelines lose the sequence-diff-corrected protein consequence
    on splice-adjacent missense.

    **Recommended: (a)**, dual-dispatch for exon-body variants.
    Pure-intronic splice variants (inside the intron, e.g. +1/+2)
    stay legacy-only since there's no protein-level diff to
    compute — the variant doesn't sit in translated sequence.

    This is the contract the classifier will honor in 3d; calling
    it out explicitly here (rather than in a commit message) so
    the splice-adjacency behaviour is reviewable alongside the
    rest of the design.
    """
    raise NotImplementedError(
        "_variant_is_splice_adjacent is WIP scaffolding — see #309.")
