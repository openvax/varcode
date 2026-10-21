"""Native ROS_Cfam_1.0 calls, using independently published coding positions.

These selected observations have different levels of somatic evidence. See
data/public_somatic/README.md; this is not a complete somatic truth set.
"""

import hashlib
import json
from pathlib import Path
import re

from Bio.Seq import Seq
from pyensembl import cached_release
import pytest

from varcode import Variant
from varcode.effects import FrameShift, Substitution


DATA = Path(__file__).parent / "data/public_somatic"
BENCHMARK = json.loads((DATA / "canine-ros-native.json").read_text())


@pytest.fixture(scope="module")
def ros_genome():
    reference = BENCHMARK["reference"]
    genome = cached_release(reference["ensembl_annotation_date"], species="dog")
    assert genome.reference_name == reference["assembly"]
    assert genome.installed(), "Run pyensembl install --release 2020_10 --species dog"
    # Pin the actual reference files, not just a mutable download URL.
    urls = [genome.gtf_url, *genome.transcript_fasta_urls, *genome.protein_fasta_urls]
    assert set(urls) == {entry["url"] for entry in reference["files"].values()}
    for entry in reference["files"].values():
        path = Path(genome.download_cache.local_path(entry["url"]))
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    return genome


@pytest.mark.parametrize("record", BENCHMARK["records"], ids=lambda row: row["id"])
def test_native_ros_cfam_published_consequences(ros_genome, record):
    window = record["reference_window"]
    sequence = window["sequence"]
    assert hashlib.sha256(sequence.encode()).hexdigest() == window["sequence_sha256"]
    assert window["id"].split(":")[1:] == [
        "ROS_Cfam_1.0", record["contig"], str(window["start"]), str(window["end"]), "1"]
    offset = record["start"] - window["start"]
    assert sequence[offset:offset + len(record["ref"])] == record["ref"]

    variant = Variant(
        record["contig"], record["start"], record["ref"], record["alt"],
        genome=ros_genome)
    assert variant.effects(raise_on_error=True)
    compared = []
    for transcript in variant.transcripts:
        if transcript.gene_name != record["gene"] or not transcript.complete:
            continue
        # All ten selected loci are on the plus strand. The papers do not
        # specify transcript IDs; compare models supporting their coding position.
        assert transcript.strand == "+"
        cds_start = min(transcript.start_codon_spliced_offsets)
        rna_start = transcript.spliced_offset(record["start"])
        coding_offset = rna_start - cds_start
        if record["published_cdna_change"]:
            change = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", record["published_cdna_change"])
            published_offset = int(change[1]) - 1
            if coding_offset != published_offset:
                continue
            assert (record["ref"], record["alt"]) == change.group(2, 3)
            rna_start = cds_start + published_offset
        else:
            # The validated SPTBN1 deletion is published as p.Asn795fs.
            assert record["published_protein_change"] == "p.N795fs"
            if coding_offset // 3 != 794:
                continue
            assert transcript.protein_sequence[794] == "N"
            assert transcript.spliced_offset(record["end"]) == rna_start + 1

        rna = transcript.sequence
        assert rna[rna_start:rna_start + len(record["ref"])] == record["ref"]
        mutant_cds = (
            rna[cds_start:rna_start] + record["alt"]
            + rna[rna_start + len(record["ref"]):])
        # Independent translation of the published cDNA edit (SNVs) or
        # genomic deletion (SPTBN1), including available downstream RNA.
        expected = str(Seq(mutant_cds[:len(mutant_cds) // 3 * 3]).translate(to_stop=True))
        effect = variant.effect_on_transcript(transcript)
        conditional = getattr(effect, "effect_if_splicing_unchanged", None)
        if conditional is not None:
            effect = conditional
        assert effect.mutant_protein_sequence == expected
        if record["published_consequence"] == "Missense":
            assert isinstance(effect, Substitution)
            if record["published_protein_change"]:
                assert effect.short_description == record["published_protein_change"]
        else:
            assert isinstance(effect, FrameShift)
            assert effect.aa_mutation_start_offset == 794
            assert effect.short_description == "p.N795fs"
        compared.append(transcript.id)
    assert compared, "No current transcript supports the published consequence: %s" % record["id"]
