# Real-world somatic-caller VCF fixtures

Tiny but format-faithful examples of what production cancer-pipeline callers
actually emit. Tests in `tests/test_real_world_callers.py` pin our parser
against hand-verified expected values for each fixture; the parametrized
oracle tests in `tests/test_vcf_parsing.py` additionally compare every row
against PyVCF3 element-by-element.

| File | Caller | Notable shapes exercised |
|---|---|---|
| `mutect2_example.vcf` | GATK4 MuTect2 | Per-alt arrays (`TLOD`/`NLOD`/`POPAF`/`MPOS` Number=A); `MBQ`/`MFRL`/`MMQ` Number=R; `RPA`/`RU`/`STR` repeat info; phasing (`PGT`/`PID`); `PL` Number=G; multi-allelic with triallelic GT `0/1/2`; full FILTER vocabulary |
| `strelka2_somatic_snvs.vcf` | Strelka2 SNV | Per-base tier counts (`AU:CU:GU:TU` Number=2); per-sample format with **no GT field**; `SGT=GG->AG` notation; `##cmdline=` metadata with paths and flags |
| `vep_annotated_csq.vcf` | VEP-annotated | `CSQ` Number=. with multi-transcript comma-split; per-record pipe-delimited subfields (23 in this fixture) including HGVS notation with `>`/`:`/`/`; Description containing pipes |

## Why these specific shapes

Each fixture was chosen for a behavior the older `mutect-example.vcf` /
`strelka-example.vcf` fixtures don't exercise:

- **MuTect2 vs. MuTect1**: The format families are unrelated — MuTect2 uses
  `TLOD`/`NLOD`/`POPAF`, `PGT`/`PID` phasing, `RPA`/`RU`/`STR` repeat info,
  and a richer FILTER vocabulary. None of that is in MuTect1 output.
- **Strelka2 vs. Strelka1**: Strelka2 dropped GT entirely from somatic SNV
  records — the per-sample format is `DP:FDP:SDP:SUBDP:AU:CU:GU:TU` and
  callers compute genotype downstream from the per-base counts.
- **VEP CSQ**: The single most common gnarly INFO field in clinical
  pipelines. Each record is comma-delimited transcripts × pipe-delimited
  fields × embedded `>`/`:`/`/`/`%xx`. Any parser that mangles commas or
  quotes ends up producing wrong output here; this fixture catches that.
