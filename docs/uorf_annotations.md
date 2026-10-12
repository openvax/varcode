# Upstream ORF variants and reference MS evidence

Install `pip install 'varcode[uorf]'` to use the shared offline catalogue from
[oncoref](https://github.com/pirl-unc/oncoref). It contains the 3,771 reference
models originally released in Isovar, including 1,219 with reported MS evidence.
Oncoref owns the empirical data, source/license notices and builder; Varcode
owns variant consequences; Isovar owns RNA reconstruction. No dependency from
Varcode to Isovar is introduced. Supplied compatible catalogue/model objects
can be used without the optional dependency.

```python
from pyensembl import cached_release
from varcode import Variant

variant = Variant('7', 66205486, 'C', 'G', genome=cached_release(95))
canonical = variant.effects()  # Existing transcript/CDS annotations.
for result in variant.orf_effects():
    print(result.orf_id, result.effect)
    annotation = result.to_dict()  # JSON-compatible, including every outcome.
```

This TPST1 +4 substitution is P2A in uORF `c7riboseqorf56`. It also changes
initiation context. The canonical TPST1 annotation remains 5′ UTR. The reference
native peptide `KIILSKISV`, residues 31–39, remains sequence-compatible with this
product; that reference observation does not establish translation or HLA
presentation of the patient's mutant protein.

Each reference model is validated against the chosen transcript's assembly,
identity, entire ordered spliced coding path, stop and translated protein.
Transcript versions/releases and catalogue provenance are exported. A newer
transcript with the same stable ID cannot inherit evidence unless the complete
path/protein checks pass. Mismatches produce explicit unresolved query results;
invalid genomic REF raises `ReferenceMismatchError`, including VCF padding.
Alleles always use forward-genomic orientation.

## Consequences and alternatives

`annotate_uorf_variants(variants, record, transcript, genome_build=...)` accepts
an explicitly known-cis group. It uses the existing transcript-model engine
for nucleotide edits, local structural edits, phase and splice alternatives.
Use `germline_variants`, `homozygous_germline`, `phase_resolver` and
`sequence_provider` as in `predict_transcript_model_effect`. A genomic sequence
provider is required where intronic bases are needed. Unknown sequence stays
unresolved. Hypothesis limits never silently truncate the alternatives.

`result.effect` is the ordinary Varcode consequence in this ORF's frame.
`result.to_dict()` retains the model, reference/alternate cDNA and proteins,
initiation contexts, stop/boundary status, every realized hypothesis and
reference peptide overlays. Native residues are mapped by three-base codon
ancestry and sequence, so matching letters in a new frame cannot borrow native
MS support. Reference peptide observations retain quality and source-mapping
ambiguity. RNA support, reference MS, mutant translation and mutant HLA
presentation remain separate.

Loss of a recognized start does not prove no initiation: alternative and
near-cognate initiation are unassessed. Stop-loss and frameshift predictions
use the available downstream transcript; a sequence boundary is not an
observed stop. Initiation context is a sequence prior, not a measured
translation efficiency. Downstream main-CDS regulation is unassessed.

## Sid RNA products

`annotate_uorf_product(record, transcript, product, genome_build=...)` accepts a
`RealizedTranscriptProduct` with one `(contig, position, strand)` or `None`
entry per nucleotide in `product.evidence['base_origins']`. Coordinates are
one-based. A uniquely retained complete reference start is required; repeated
or unavailable starts remain unresolved. This uses the realized RNA sequence,
including its actual fusion/splice path, without adding canonical sequence
that the reads do not contain.

Isovar's `annotate_sv_uorf_evidence(result, transcripts)` adapts its zero-based
RNA paths to this interface and preserves candidate-specific full-interval
support in its JSON/TSV exports. With the packaged original Sid PacBio reads,
TPST1–CRCP has 30-, 61-, and 72-aa alternatives with respectively 15, 1, and 1
full-interval fragment witnesses. All lose the native peptide at residues
31–39. None receives MS support for its altered tail. Shared signal ancestry,
missing qualities and unresolved RNA polarity stay in the source evidence.

The catalogue queries reference models only. An empty result does not rule
out a newly created start or another translated ORF. De novo start discovery,
complete contemporary shotgun curation, and general DNA/RNA hypothesis
composition remain tracked in [#583](https://github.com/openvax/varcode/issues/583)
and [Isovar #466](https://github.com/openvax/isovar/issues/466).

Primary scientific references: [Mudge et al. 2022](https://doi.org/10.1038/s41587-022-01369-0),
[Deutsch et al. 2026](https://doi.org/10.1038/s41586-026-10459-x), and
[Kozak's short-leader experiments](https://pubmed.ncbi.nlm.nih.gov/1820208/).
The catalogue's pinned source manifest defines the actual evidence coverage.
