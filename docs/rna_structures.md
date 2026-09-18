# Importing observed RNA structures

Attach observed RNA structures to [structural variants](structural_variants.md)
without treating an RNA junction as a DNA breakpoint. The examples below assume
that the DNA variants and relevant reference transcripts are already loaded.

Varcode represents supplied structures and predicts their coding consequences;
RNA reconstruction and read-based reconciliation belong to Isovar or another
RNA producer. See [the library responsibilities](library_roles.md) for the split
and the remaining end-to-end SV workflow.

## Import Exacto transcripts

Use an existing DNA `StructuralVariant` as the anchor, so an RNA splice junction
is not mistaken for a genomic breakpoint. Multiple observed RNA models can be
attached to that same variant:

```python
from varcode import load_exacto_fusions

rna = load_exacto_fusions(
    "transcript_structures.tsv", "integrated_variants.tsv",
    variants_by_id={"42": structural_variant},
)
effects = variants.effects(rna_resolver=rna)
for candidate in rna.candidates:
    sequence = candidate.effect.mutant_transcript.cdna_sequence
    model_id = candidate.evidence["transcript_model_id"]
```

`variants_by_id` maps **Exacto DNA call IDs**, not assumed VCF record IDs, to
already loaded variants. Only these IDs are selected from the integration file.
The model key includes both `transcript_model_id` and the reference transcript
ID set. Different models/reference interpretations and their DNA links are kept;
duplicate integration rows do not duplicate a candidate.

The adapter follows Exacto's `index` ordering and inclusive `read_start` /
`read_end` base-row coordinates. Sequences already run 5' to 3' even for reverse
strand genes: they are never reverse-complemented a second time. Event rows
preserve splice structure but do not invent missing bases. The complete original
structure and integration rows are in `candidate.evidence` and the mutant
transcript's evidence. Row counts are **not** supporting-read counts.

This is a deliberately limited adapter for linear, two-locus, SV-linked models
with an annotated sense 5' anchor. A missing or antisense 3' annotation gives
`TranslocationToIntergenic`, whose existing meaning includes non-sense joins.
Unknown transcript IDs and incomplete/unsupported structures raise errors.
This does not implement all Exacto DNA/RNA variant formats, circular RNA,
or multi-gene paths. Those remain in the broader
[#260](https://github.com/openvax/varcode/issues/260) roadmap.

## Other RNA producers

For other producers, use the same existing effect classes directly:

```python
from varcode import RNAEvidence, make_fusion_outcome

candidate = make_fusion_outcome(
    structural_variant, five_prime_transcript,
    partner_transcript=three_prime_transcript,  # only when sense-oriented
    sequence=observed_sequence,               # already 5' to 3'
    transcript_model_id="assembled-model-7",
    source="rna_assembler", read_count=6,
    extra_evidence={"sequence_status": "junction_fragment"},
)
effects = variants.effects(rna_resolver=RNAEvidence([candidate]))
```

Do not supply a sense partner merely because a gene name appears near a
breakpoint. Leave `partner_transcript=None` for unresolved/intergenic/antisense
partners and record the observed loci/orientation in `extra_evidence`.
The importer preserves the supplied sequence, including insertions, without
appending reference exons. Read counts do not become probabilities or change
default candidate ranking. DNA predictions and RNA observations coexist.

Protein sequence is `None` by default. A caller can pass `cds_start=<0-based
offset>` to `make_fusion_outcome`, or a `cds_starts` mapping keyed by
`("model-id", ("sorted-reference-id-1", "sorted-reference-id-2"))` to the Exacto
loader. This requires a valid unambiguous start-to-stop ORF in the supplied
sequence. It produces a **sequence-predicted protein**, not evidence of
translation, full-length RNA, tumor-cell identity, or a mature two-gene fusion.

## Import Exacto's protein predictions

When Exacto has already selected reading frames, supply its native
primary-structures table instead of choosing `cds_starts` yourself:

```python
rna = load_exacto_fusions(
    "transcript_structures.tsv", "integrated_variants.tsv",
    variants_by_id=variants_by_id,
    primary_structures_path="primary_structures.tsv",
)
for candidate in rna.candidates:
    print(candidate.evidence.get("exacto_peptide_id"))
    print(candidate.effect.mutant_protein_sequence)
    print(candidate.evidence.get("protein_completeness"))
```

Each `(model, reference-transcript group, peptide_id)` remains separate. The
loader validates native `primary_structure_index`, `amino_acid_index`, and
`codon_index`, checks each nucleotide against its `transcript_structure_index`
and inclusive read coordinates, and checks the amino acid against its codon.
The amino acid repeated on all three base rows is emitted once. Terminal `*`
is retained in provenance but omitted from `mutant_protein_sequence`.

Partial peptides are available, with `protein_completeness` set to
`partial_start`, `partial_end`, or `partial_both`; start-to-stop predictions
are labeled `start_to_stop`. Trailing incomplete codons cannot claim an amino
acid. Missing peptide rows leave that model's protein `None`. Per-base variant
IDs, frameshift state and all original fields remain in `exacto_primary_structure`.
Exacto's standard genetic-code choice is recorded as `protein_translation_table=1`.

These sequences are **predictions from RNA**, not evidence of translation.
A complete start-to-stop ORF need not cross the rearrangement or establish a
full-length fusion transcript. Downstream users must check completeness and
source coordinates rather than treating every protein string as a complete
expressed fusion protein.

Known limitation: SV sequence-change flags currently compare imported partial
proteins with the full reference ([#462](https://github.com/openvax/varcode/issues/462)).
Until that is fixed, a positive flag alone does not establish a change in a
partial observation; check its mapped observed region and completeness.

The small osteosarc regression fixtures deliberately exercise the negative
case: GABBR1 joins sequence upstream of SLC29A1, OTUD7A joins an antisense FMN1
intron, and the KLF15-side reads are intronic. Their RNA junctions are retained
without manufacturing coding fusions or protein sequences.

Format reference: [Exacto's structure translation source, pinned revision
307c086](https://github.com/pirl-unc/exacto/blob/307c08670d5e706734bddf393bcebc84db497f9f/exacto/exacto-translator/src/algorithms/translation.rs).

## Related reference

- [Transcript models](transcript_models.md): partial structures and sequence access.
- [RNA and transcript API](api_rna.md): import parameters and evidence objects.
- [Phasing](phasing.md): evidence linking variants on the same allele.
