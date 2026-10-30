"""Small native reference slices; genomic coordinates are never renumbered."""

import gzip
import json
from pathlib import Path

import pytest

from varcode import Genome

DATA = Path(__file__).parent / "data/native_references"
CASES = json.loads((DATA / "cases.json").read_text())


class ReferenceWindows:
    """Sparse, original-coordinate DNA; unavailable intervals remain unavailable."""

    def __init__(self, windows, contig=None):
        self.windows = windows
        self.contig = contig

    def __getitem__(self, key):
        if isinstance(key, str):
            return ReferenceWindows(self.windows, key)
        for window in self.windows:
            if (
                window["contig"] == self.contig
                and key.start >= window["start"] - 1
                and key.stop <= window["end"]
            ):
                offset = window["start"] - 1
                return window["sequence"][key.start - offset : key.stop - offset]
        return ""


@pytest.fixture(scope="module")
def native_genomes(tmp_path_factory):
    from orfdata._pinned_genome import PinnedEnsemblGenome

    windows = json.loads(gzip.decompress((DATA / "dna-windows.json.gz").read_bytes()))
    result = {}
    for name in ("mouse", "dog", "chm13"):
        metadata = json.loads((DATA / (name + ".json")).read_text())
        metadata.pop("sources")
        directory = tmp_path_factory.mktemp("native-" + name)
        genome = PinnedEnsemblGenome(
            **metadata,
            gtf_path_or_url=str(DATA / (name + ".gtf.gz")),
            transcript_fasta_paths_or_urls=[str(DATA / (name + ".cdna.fa.gz"))],
            protein_fasta_paths_or_urls=[str(DATA / (name + ".pep.fa.gz"))],
            cache_directory_path=str(directory),
        )
        genome.index()
        result[name] = genome
        available = [window for window in windows if window["genome"] == name]
        if available:
            result[name + "_dna"] = Genome(
                genome, fasta=ReferenceWindows(available), verify=False
            )
    yield result
    for name, genome in result.items():
        if not name.endswith("_dna"):
            genome.close()
