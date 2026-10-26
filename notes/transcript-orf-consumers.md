# Common transcript ORF consumers (#638, ORFdata #61/#64)

Treat TranscriptORF as the coding placement, with optional catalogue provenance.
Use its coordinates, source identity, protein, code and evidence directly. Keep
annotated CDSs, source-defined starts and bounded fragments on their original RNA.

Complete models with ordinary, reference-validated decoding use the existing
transcript-model engine. Partial models retain the known fragment and explicit
unavailable operations; known internal coding changes can be calculated in a
validated frame without asserting full-protein reconstruction. Unsupported
recoding, circular paths, unknown frame or unmapped source models remain native
annotations with a reason, not dropped records or fabricated catalogue entries.

Persist the common placement with its archived RNA in the existing native effect
graph. Proteome export must preserve kind, identifiers, boundaries, unavailable
outcomes and evidence. Stamp producer method/version at the actual annotation
entry point (#639), not by guessing inside an exporter.

Validate ordinary public effects, source collisions, both strands, splice codons,
native non-AUG/multiple-start/partial examples, wrong references and phases, and
offline graph/proteome round trips. Existing catalogue behavior stays covered.

Source conventions: INSDC feature-table (codon_start and partial locations),
https://www.insdc.org/submitting-standards/feature-table/; NCBI genetic codes,
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
