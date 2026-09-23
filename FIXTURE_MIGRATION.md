# Osteosarc fixture compatibility

Pin the current optional snapshot adapter to Osteosarc 0.2.3. Preserve the
bundled historical snapshot, original VCF recipes, reference release and allele
expectations. The historical JSON export still records its original Osteosarc
version; regeneration of that artifact remains intentionally pinned separately.

Varcode currently has no BAM-subset builder. Ordinary allele/protein-only tests
continue without Osteosarc or BAM acquisition. Run the optional snapshot suite
offline against the shared release before publishing this patch.
