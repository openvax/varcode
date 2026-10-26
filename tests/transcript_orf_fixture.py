"""Pinned Ensembl 115 annotation/cDNA/protein slices shared with ORFdata's contract tests."""

import json
from pathlib import Path

import pytest
from pyensembl import Genome

DATA = json.loads((Path(__file__).parent / "data/transcript_orfs_ensembl115.json").read_text())


@pytest.fixture(scope="module")
def orf_genome(tmp_path_factory):
    path = tmp_path_factory.mktemp("transcript-orfs")
    rows, genes, cdna, pep = [], {}, [], []
    for t in DATA["transcripts"]:
        gene = (f'gene_id "{t["gene_id"]}"; gene_version "{t["gene_version"]}"; gene_name "{t["gene_name"]}"; '
                f'gene_biotype "{t["gene_biotype"]}";')
        attrs = gene + f' transcript_id "{t["id"]}"; transcript_version "{t["version"]}"; transcript_biotype "{t["biotype"]}";'
        start, end = t["exons"][0][0], t["exons"][-1][1]
        span = genes.get(t["gene_id"], (t["contig"], t["strand"], start, end, gene))
        genes[t["gene_id"]] = (t["contig"], t["strand"], min(start, span[2]), max(end, span[3]), gene)
        if t.get("native_gtf"):
            rows += [row["text"] + "\n" for row in t["native_gtf"]]
        else:
            line = f'{t["contig"]}\tensembl\t{{}}\t{{}}\t{{}}\t.\t{t["strand"]}\t.\t{attrs}{{}}\n'
            rows.append(line.format("transcript", start, end, ""))
            rows += [line.format("exon", a, b, f' exon_number "{i+1}"; exon_id "{t["id"]}-{i}";')
                     for i, (a, b) in enumerate(sorted(t["exons"], reverse=t["strand"] == "-"))]
            protein = f' protein_id "{t["protein_id"]}"; protein_version "{t["protein_version"]}";'
            rows += [line.format("CDS", a, b, protein) for a, b in t["cds"]]
            rows += [line.format(feature, a, b, "") for feature in ("start_codon", "stop_codon")
                     for a, b in t[feature]]
        cdna.append(f'>{t["id"]}.{t["version"]}\n{t["sequence"]}\n')
        pep.append(f'>{t["protein_id"]}.{t["protein_version"]}\n{t["protein"]}\n')
    rows += [f"{c}\tensembl\tgene\t{a}\t{b}\t.\t{s}\t.\t{attrs}\n"
             for c, s, a, b, attrs in genes.values()]
    (path / "slice.gtf").write_text("".join(rows))
    (path / "cdna.fa").write_text("".join(cdna))
    (path / "pep.fa").write_text("".join(pep))
    genome = Genome("GRCh38", "ensembl", "115", gtf_path_or_url=str(path / "slice.gtf"),
                    transcript_fasta_paths_or_urls=[str(path / "cdna.fa")],
                    protein_fasta_paths_or_urls=[str(path / "pep.fa")], cache_directory_path=str(path))
    genome.species = "homo_sapiens"
    genome.index()
    yield genome
    genome.close()
