"""Build the offline genetic-code snapshot from NCBI gc.prt and Taxonomy XML.

Usage: python notes/build_genetic_code_data.py gc.prt taxonomy.xml output.json
Sources: https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt
and EUtils efetch.fcgi?db=taxonomy&id=<taxon IDs>&retmode=xml.
No network access at build or annotation time; review refreshed data in a PR.
"""

import hashlib
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


def build(codes_path, taxonomy_path):
    codes_data = Path(codes_path).read_bytes()
    taxonomy_data = Path(taxonomy_path).read_bytes()
    codes = {}
    for block in re.findall(r"\{([^{}]+)\}", codes_data.decode()):
        identifier = re.search(r"\bid\s+(\d+)", block)
        if identifier is None:
            continue
        codes[identifier[1]] = dict(
            name=re.search(r'name\s+"([^"]+)"', block)[1],
            amino_acids=re.search(r'ncbieaa\s+"([^"]+)"', block)[1],
            starts_stops=re.search(r'sncbieaa\s+"([^"]+)"', block)[1])
    taxa = {}
    for taxon in ET.fromstring(taxonomy_data).findall("Taxon"):
        aliases = {taxon.findtext("ScientificName")}
        for tag in ("Synonym", "GenbankCommonName", "CommonName", "EquivalentName"):
            aliases.update(n.text for n in taxon.findall("OtherNames/" + tag))
        codes_by_compartment = dict(
            nuclear=int(taxon.findtext("GeneticCode/GCId")),
            mitochondrial=int(taxon.findtext("MitoGeneticCode/MGCId")))
        for prop in taxon.findall("Properties/Property"):
            if prop.findtext("PropName") == "pgcode":
                codes_by_compartment["plastid"] = int(prop.findtext("PropValueInt"))
        taxa[taxon.findtext("TaxId")] = dict(
            name=taxon.findtext("ScientificName"),
            aliases=sorted(name.lower().replace(" ", "_") for name in aliases),
            **codes_by_compartment)
    return dict(
        source="NCBI gc.prt version 4.6 and Taxonomy EUtils; retrieved 2026-10-20",
        code_source_sha256=hashlib.sha256(codes_data).hexdigest(),
        taxonomy_source_sha256=hashlib.sha256(taxonomy_data).hexdigest(),
        codes=codes, taxa=taxa)


if __name__ == "__main__":
    Path(sys.argv[3]).write_text(json.dumps(build(*sys.argv[1:3]), indent=2, sort_keys=True) + "\n")
