# Bounded discovery of newly created start codons

Add an opt-in product consumer, `discover_created_orfs(effect, policy=...)`, and
a lower-level paired-product API. Reuse the transcript-model outcomes rather
than reconstructing phase or splicing. Scan only start triplets changed in
sequence or adjacency, comparing mapped base ancestry with the paired baseline.
Coordinate shifts alone must not create a new start. Ambiguous ancestry and
missing products remain explicit; no sequence alignment guesses missing origins.

The policy declares accepted start codons (ATG by default), maximum scanned RNA,
maximum translated codons, and maximum candidates. Limits and partial products
are returned, including an empty result's assessment scope. Genetic codes come
from the selected transcript/reference. Kozak context is descriptive evidence.

Each candidate is a normal `MutationEffect` with host transcript/gene identity,
baseline and mutant RNA accessors, and a conditional mutant peptide. A focal
`discovered_orf` field describes the mutant coding interval, ancestry, genetic
code and termination; it is not a fabricated reference ORF. The baseline protein
for this new initiation hypothesis remains unknown, rather than borrowing the
canonical protein or asserting no translation. Preserve source hypotheses,
discovery policy and reference identity through JSON/pickle/CSV and Proteome.

Initial support covers paired mapped products, including supplied noncoding RNA
and realized coding-transcript products. Legacy results without paired RNA,
unmapped structural products and unknown/partial ancestry report unavailable
discovery. This does not infer expression, downstream CDS repression, initiation
probabilities, recoding, or all alternative starts.

Validate with the published CDKN2A c.-34G>T start gain (Liu et al., 1999,
doi:10.1038/5082), independently translated pinned RNA, plus controlled tests
for both strands, exon-junction starts, indels, phase, mitochondrial code,
ambiguous/missing sequence, bounds, and lossless exports. Supplied ORF annotation
remains covered by the existing common-field and native ORF regressions.
