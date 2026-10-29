# Validate consumed SV INFO metadata (#657)

Normalize only consumed structural-variant INFO fields at the parser boundary.
Reject Boolean/invalid scalar SVTYPE and END values and malformed confidence
intervals with actionable ValueError messages. CIPOS/CIEND must be ordered pairs
of signed integer offsets (VCF 4.5, section 3); preserve valid negative offsets.
Accept existing scalar and singleton-list parser representations without
mutating the supplied mapping. Keep literal nonstructural ALT dispatch unchanged.

Use one helper for symbolic and breakend branches, with direct parser and VCF
loader regressions. Preserve loader strictness/row-error policy. Do not add an
Isovar-only workaround. Update relevant parser documentation and release notes,
bump a patch release, run lint/full tests, review the diff, require green CI,
merge and deploy from clean main; then consume the released dependency in the
Isovar source-validation PR.

Primary source: https://samtools.github.io/hts-specs/VCFv4.5.pdf (25 Feb 2026).
