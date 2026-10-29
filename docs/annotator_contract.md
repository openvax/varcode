# Writing an annotator

An annotator predicts an effect for a variant on a transcript. Register a custom
implementation when you want to compare your model through Varcode's public
annotation interface. For the built-in implementations, see
[Annotation pipeline and compatibility](experimental_annotators.md).

## Required method

An annotator needs a `name` and an `annotate_on_transcript` method. It may also
provide a `version` for provenance.

Return a `MutationEffect` for a supported input, or Python's `NotImplemented`
singleton for an input your implementation cannot handle. No capability list
is required.

For example, this adapter takes a callable that returns an effect or `None`:

```python
from collections.abc import Callable
from types import NotImplementedType

from pyensembl import Transcript
from varcode import MutationEffect, Variant


class MyAnnotator:
    name = "my_model"

    def __init__(self, model: Callable[[Variant, Transcript], MutationEffect | None]):
        self.model = model

    def annotate_on_transcript(
            self, variant: Variant, transcript: Transcript,
    ) -> MutationEffect | NotImplementedType:
        prediction = self.model(variant, transcript)
        if prediction is None:
            return NotImplemented
        return prediction
```

Register an instance, then select it explicitly:

```python
import varcode

annotator = MyAnnotator(my_model)
varcode.register_annotator(annotator)
effects = variant.effects(annotator=annotator.name)
```

Here `my_model` is your prediction callable and `variant` is a loaded variant.

The annotations let a type checker verify the required EffectAnnotator
protocol without making your class inherit from it. `register_annotator()`
returns the same concrete instance. For concurrent work, pass `annotator=`
per call; `use_annotator()` temporarily changes the process-wide default.

## Unsupported inputs and errors

Public APIs convert `NotImplemented` into
`Unresolved(mechanism="unsupported_annotation", reason=...)`, preserving the
selected annotator's provenance. They do not silently fall back to the default.
For example, selecting `protein_diff` for an SV returns `Unresolved`, including
inside `use_annotator("protein_diff")`.

Returning `None` directly is a plugin error. Exceptions retain the ordinary
`raise_on_error` behavior; do not turn unexpected failures into unsupported inputs.

## Germline context

An optional `annotate_with_context(variant, transcript, germline_ctx,
phase_resolver=None)` method has the same return contract. Without it, nonempty
germline context produces `Unresolved` rather than being ignored or sent to
another implementation. Empty context uses `annotate_on_transcript`.

`effects()`, `effect_on_transcript()`, and
`predict_variant_effect_on_transcript()` use the same selection rules, including
the scoped default.

## Joint haplotypes

To predict combined effects, implement the optional method:

```python
def annotate_haplotype(self, variants, transcript,
                       germline_ctx=None, phase_resolver=None):
    # variants is a tuple of two or more known-cis variants on this transcript.
    return NotImplemented
```

The collection groups variants using the phase resolver, then calls the selected
annotator for each group. Return a `MutationEffect` or `NotImplemented`, just as
for individual annotation. Missing hooks and refusals produce
`Unresolved(mechanism="unsupported_haplotype")`; `None` is a plugin error.
Exceptions follow `raise_on_error`, with handled errors retained as `Failure`.
There is no automatic fallback to another annotator.

Joint results coexist with individual effects and carry `variants`,
`phase_source`, `annotator`, and `annotator_version`. Handle all supplied context
or decline the group; ignoring germline changes would produce the wrong baseline.
JSON round trips retain this group context and any attached observed RNA models.
The phase resolver may also supply observed mutant transcripts, whose source and
coverage must remain distinguishable from inferred sequences.

The default uses its multi-variant engine and returns an ordinary classified
effect with alternatives in `candidates`. Germline lookup covers the transcript
span. Supplied RNA models refine the same candidates, preserving partial or
unmatched models as unresolved comparisons. No observation sidecar is needed.
The explicit `fast` and `protein_diff` selections retain the legacy point-edit
builder and `HaplotypeEffect` result type.

## Default implementation

The default centralizes genomic layout editing, paired splice enumeration,
transcript realization and classification. Existing helpers in
`varcode.effects.structural` own typed junctions, BNDs, inversions, nonlocal
rearrangements and supplied transcript assemblies. Unsupported inputs
remain explicit; there is no automatic fallback between annotators.

Structural transcript models can retain the builder label `structural_variant`
in their provenance, while the effect collection records the selected annotator.
A stored builder label is not a selectable annotator name.

## Testing an integration

Test supported inputs, explicit refusal, context-dependent refusal, and genuine
errors separately. Check direct transcript calls as well as collections, and
preserve provenance and unresolved results when serializing.

Compare implementations only on inputs they both support. Passing a parity
test for point edits does not establish SV support. When sending results to
other tools, retain candidate evidence, transcript identity, and sequence
completeness; a protein string alone does not establish a complete expressed
protein.
