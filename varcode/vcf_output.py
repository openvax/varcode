# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

"""VCF export with explicit sample, field and source-allele ordering."""

from collections import defaultdict
from numbers import Integral, Real
from math import comb
import sys

from .genotype import parse_gt_string
from .vcf_parsing import FieldDef, RESERVED_FORMAT, RESERVED_INFO


_FORMAT_NUMBERS = {
    "GT": 1, "DP": 1, "GQ": 1, "PS": 1, "PQ": 1, "FT": 1, "MQ": 1,
    "AD": -3, "ADF": -3, "ADR": -3, "AF": -1,
    "PL": -2, "GL": -2, "GP": -2,
}
_FORMAT_TYPES = dict(RESERVED_FORMAT, AD="Integer", ADF="Integer", ADR="Integer", AF="Float")


def _value(value):
    if value is None:
        return "."
    if isinstance(value, (list, tuple)):
        return ",".join(_value(v) for v in value) if value else "."
    return str(value)


def _filter(value):
    if value is None or value == ".":
        return "."
    if isinstance(value, (list, tuple)):
        return ";".join(value) if value else "PASS"
    if isinstance(value, str):
        return value
    raise ValueError("Invalid VCF filter value: %r" % (value,))


def _field_def(name, values, kind):
    """Infer types for new metadata; source header declarations take priority."""
    flat = [v for value in values
            for v in (value if isinstance(value, (list, tuple)) else [value])
            if v is not None]
    reserved = _FORMAT_TYPES if kind == "FORMAT" else RESERVED_INFO
    if name in reserved:
        field_type = reserved[name]
    elif flat and all(isinstance(v, bool) for v in flat) and kind == "INFO":
        field_type = "Flag"
    elif flat and all(isinstance(v, Integral) and not isinstance(v, bool) for v in flat):
        field_type = "Integer"
    elif flat and all(isinstance(v, Real) and not isinstance(v, bool) for v in flat):
        field_type = "Float"
    else:
        field_type = "String"
    if field_type == "Flag":
        number = 0
    elif kind == "FORMAT" and name in _FORMAT_NUMBERS:
        number = _FORMAT_NUMBERS[name]
    else:
        number = None if any(isinstance(v, (list, tuple)) for v in values) else 1
    return FieldDef(name, number, field_type, "Exported by Varcode; inferred from metadata")


def _records(variants, metadata):
    """Restore complete source ALT lists without merging unrelated record IDs."""
    groups = defaultdict(list)
    for variant in variants:
        if not all(hasattr(variant, name) for name in ("original_ref", "original_alt")):
            raise ValueError("VCF export currently requires small Variant records")
        meta = metadata[variant]
        source_alts = tuple(meta.get("vcf_alt_alleles") or ())
        key = (variant.original_contig, variant.original_start, variant.original_ref,
               meta.get("id"), source_alts or (variant.original_alt,))
        groups[key].append((variant, meta))
    records = []
    for key, entries in groups.items():
        first, meta = entries[0]
        source_alts = tuple(meta.get("vcf_alt_alleles") or ())
        if source_alts:
            indexes = set()
            for variant, item in entries:
                index = item.get("alt_allele_index")
                if (not isinstance(index, int) or not 0 <= index < len(source_alts)
                        or source_alts[index] != variant.original_alt):
                    raise ValueError("VCF ALT provenance is inconsistent at %s" % (key,))
                indexes.add(index)
                for field in ("id", "qual", "filter", "info", "sample_info"):
                    if item.get(field) != meta.get(field):
                        raise ValueError("Conflicting source metadata at %s" % (key,))
            if indexes != set(range(len(source_alts))):
                raise ValueError(
                    "Cannot export an incomplete source ALT list at %s; "
                    "retain all alleles or explicitly remap GT and allele-indexed fields" % (key,))
            alts = source_alts
        else:
            # Old/manually constructed biallelic records remain supported.
            # An index alone cannot prove that the other source ALTs were
            # retained, even if GT only happens to refer to ALT 1.
            if "alt_allele_index" in meta:
                raise ValueError("Missing source ALT list; reload the VCF before exporting indexed alleles")
            alts = (first.original_alt,)
        for sample, fields in (meta.get("sample_info") or {}).items():
            gt = fields.get("GT")
            ploidy = None
            if gt is not None:
                alleles, _ = parse_gt_string(gt)
                if any(a is not None for a in alleles):
                    ploidy = len(alleles)
                if any(a is not None and (a < 0 or a > len(alts)) for a in alleles):
                    raise ValueError("GT for sample %s refers to an unavailable ALT allele" % sample)
            for field in ("AD", "ADF", "ADR", "AF"):
                value = fields.get(field)
                if value is not None and value != ".":
                    n = len(value) if isinstance(value, (tuple, list)) else 1
                    expected = len(alts) + (field != "AF")
                    if n != expected:
                        raise ValueError("%s for sample %s does not match the ALT list" % (field, sample))
            if ploidy is not None:
                for field in ("PL", "GL", "GP"):
                    value = fields.get(field)
                    if value is not None and value != ".":
                        n = len(value) if isinstance(value, (tuple, list)) else 1
                        if n != comb(len(alts) + ploidy, ploidy):
                            raise ValueError("%s for sample %s does not match ALT/ploidy" % (field, sample))
        records.append((first, alts, meta))
    return sorted(records, key=lambda row: (
        str(row[0].original_contig), int(row[0].original_start), row[0].original_ref, row[1]))


def variants_to_vcf(variants, variant_to_metadata, out=None, header=None):
    """Write small variants with stable sample identities and FORMAT fields.

    Parameters
    ----------
    variants : iterable of Variant
        Complete source ALT groups are required for multi-allelic records.
        Filtered groups are rejected rather than changing genotype meaning.
    variant_to_metadata : dict
        Variant-to-metadata mapping, normally from a loaded VariantCollection.
        ``load_vcf`` retains source ALT order in ``vcf_alt_alleles``.
    out : writable stream, optional
        Defaults to the current sys.stdout. Inputs are validated before writing.
    header : VCFHeader, optional
        Original INFO/FORMAT declarations to retain their types and cardinality.
        Otherwise declarations are inferred from values and standard FORMAT
        fields. Other source headers are not preserved; this is not a lossless
        archive of the source VCF.
    """
    variants = list(variants)
    records = _records(variants, variant_to_metadata)
    references = {v.original_reference_name for v in variants}
    if len(references) > 1:
        raise ValueError("Cannot create VCF for variants with multiple reference genomes: %s" % references)
    samples = sorted({sample for _, _, meta in records
                      for sample in (meta.get("sample_info") or {})})
    declarations = {"INFO": {}, "FORMAT": {}}
    for kind in declarations:
        values = defaultdict(list)
        for _, _, meta in records:
            dictionaries = ([meta.get("info") or {}] if kind == "INFO"
                            else (meta.get("sample_info") or {}).values())
            for fields in dictionaries:
                for name, value in fields.items():
                    values[name].append(value)
        original = (getattr(header, "info_fields" if kind == "INFO" else "format_fields")
                    if header is not None else {})
        for name, observed in sorted(values.items()):
            declarations[kind][name] = original.get(name) or _field_def(name, observed, kind)
    lines = ["##fileformat=VCFv4.2"]
    if references:
        lines.append("##reference=%s" % next(iter(references)))
    for kind, fields in declarations.items():
        for name, field in fields.items():
            number = {None: ".", -1: "A", -2: "G", -3: "R"}.get(field.number, field.number)
            description = field.description.replace("\\", "\\\\").replace('"', '\\"')
            lines.append('##%s=<ID=%s,Number=%s,Type=%s,Description="%s">' % (
                kind, name, number, field.type, description))
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        columns += ["FORMAT"] + samples
    lines.append("\t".join(columns))
    for variant, alts, meta in records:
        info = []
        for key, value in sorted((meta.get("info") or {}).items()):
            if declarations["INFO"][key].type == "Flag":
                if value:
                    info.append(key)
            else:
                info.append("%s=%s" % (key, _value(value)))
        row = [str(variant.original_contig), str(variant.original_start),
               _value(meta.get("id")), variant.original_ref, ",".join(alts),
               _value(meta.get("qual")), _filter(meta.get("filter")), ";".join(info) or "."]
        if samples:
            sample_info = meta.get("sample_info") or {}
            fields = sorted({key for info in sample_info.values() for key in info},
                            key=lambda key: (key != "GT", key))
            row.append(":".join(fields) or ".")
            for sample in samples:
                values = sample_info.get(sample) or {}
                row.append(":".join(_filter(values.get(key)) if key == "FT"
                                    else _value(values.get(key)) for key in fields) or ".")
        lines.append("\t".join(row))
    if out is None:
        out = sys.stdout
    out.write("\n".join(lines) + "\n")
