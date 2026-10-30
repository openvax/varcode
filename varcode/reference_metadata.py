"""Reference identity in flat exports, independent of local data paths.

An assembly label alone does not identify a transcript annotation. Readers
validate every recorded field before reannotating and never convert coordinates.
"""

from pyensembl import cached_release

from .reference import infer_genome


REFERENCE_COLUMNS = (
    "reference_name", "annotation_name", "annotation_version",
    "assembly_accession", "species",
)


def reference_metadata(genome):
    """Return portable assembly and annotation identity (unknown fields stay None)."""
    result = {name: getattr(genome, name, None) for name in REFERENCE_COLUMNS}
    species = result["species"]
    result["species"] = getattr(species, "latin_name", species)
    # CSV readers must not turn a release such as 2020_10 into a number.
    if result["annotation_version"] is not None:
        result["annotation_version"] = str(result["annotation_version"])
    return result


def common_reference_metadata(genomes):
    """Header fields shared by every row; heterogeneous identities stay on rows."""
    records = [reference_metadata(genome) for genome in genomes]
    if not records:
        return {}
    return {name: records[0][name] for name in REFERENCE_COLUMNS
            if all(record[name] == records[0][name] for record in records)}


def csv_genomes(frame, header, genome=None):
    """Resolve and validate one genome per CSV row, preserving annotation versions.

    Custom annotations require their installed Genome explicitly. Ordinary
    Ensembl releases can be reconstructed from their species and exact release.
    Old CSVs without annotation identity keep their previous assembly inference.
    """
    supplied = None if genome is None else infer_genome(genome)[0]
    resolved = {}
    result = []
    for _, row in frame.iterrows():
        identity = {}
        for name in REFERENCE_COLUMNS:
            value = row.get(name)
            if value is not None and str(value) not in ("", "nan", "<NA>"):
                identity[name] = str(value)
                if header.get(name) and str(header[name]) != identity[name]:
                    raise ValueError("CSV row %s=%r disagrees with header (%r)"
                                     % (name, identity[name], header[name]))
            elif header.get(name):
                identity[name] = header[name]
        key = tuple(sorted(identity.items()))
        if key not in resolved:
            found = supplied
            if found is None:
                if identity.get("annotation_name") == "ensembl" and identity.get("annotation_version"):
                    version = identity["annotation_version"]
                    found = cached_release(int(version) if version.isdecimal() else version,
                                           species=identity.get("species", "human"))
                elif identity.get("annotation_name") or identity.get("annotation_version"):
                    raise ValueError("CSV uses a custom annotation; pass its installed genome explicitly: %s"
                                     % identity)
                elif identity.get("reference_name"):
                    found = infer_genome(identity["reference_name"])[0]
                else:
                    raise ValueError("from_csv needs a genome or recorded reference_name")
            actual = reference_metadata(found)
            for name, expected in identity.items():
                if str(actual[name]) != expected:
                    raise ValueError("CSV %s=%r disagrees with supplied genome (%r); "
                                     "coordinate conversion and annotation changes must be explicit"
                                     % (name, expected, actual[name]))
            resolved[key] = found
        result.append(resolved[key])
    return result
