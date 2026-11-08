"""Rebuild native reference slices and independent sequence expectations.

Requires installed mouse 115, dog 2020_10, ORFdata CHM13, Biopython and openpyxl.
See tests/data/native_references/README.md for public inputs and evidence limits.
No Varcode prediction is used to select transcripts or calculate expectations.
"""

import argparse
from pathlib import Path
import gzip, hashlib, json, re
from pyensembl import cached_release
from orfdata import load_genome
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sources", type=Path, required=True)
parser.add_argument("--chm13-cache", type=Path, required=True)
parser.add_argument("--chm13-dna", type=Path, required=True)
parser.add_argument("--dna-contig-prefix", default="")
parser.add_argument(
    "--output",
    type=Path,
    default=Path(__file__).resolve().parents[1] / "tests/data/native_references",
)
args = parser.parse_args()
root, src = args.output, args.sources
root.mkdir(parents=True, exist_ok=True)


def compressed(path):
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def sha256(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


from Bio import SeqIO


def build_extra_animals(root, src):
    # Separate mitochondrial reference: not a mitochondrial component of ROS_Cfam.
    r = SeqIO.read(src / "dog-mito-NC_002008.4.gb", "genbank")
    records = []
    cases = []
    for gene, p, ref, alt, published in [
        ("COX2", 7593, "T", "C", "p.M187T"),
        ("ND5", 11998, "T", "C", "p.M74T"),
        ("ND6", 13791, "T", "C", "p.I106V"),
        ("ATP6", 8369, "C", "T", "p.P136S"),
    ]:
        f = next(
            f
            for f in r.features
            if f.type == "CDS" and f.qualifiers.get("gene") == [gene]
        )
        sequence = str(f.extract(r.seq))
        protein = f.qualifiers["translation"][0]
        assert str(r.seq[p - 1]) == ref
        genomic = r.seq[: p - 1] + alt + r.seq[p:]
        mutant = str(f.extract(genomic))
        translated = str(Seq(mutant).translate(table=2, to_stop=True))
        translated = "M" + translated[1:]
        records.append(
            dict(
                id=r.id + ":" + gene,
                gene=gene,
                start=int(f.location.start) + 1,
                end=int(f.location.end),
                strand="+" if f.location.strand == 1 else "-",
                sequence=sequence,
                protein=protein,
                protein_id=f.qualifiers["protein_id"][0],
                translation_table=int(f.qualifiers["transl_table"][0]),
            )
        )
        cases.append(
            dict(
                id="dog-mito-" + gene,
                position=p,
                ref=ref,
                alt=alt,
                published_hgvs=published,
                expected_cdna=mutant,
                expected_protein=translated,
                source="https://doi.org/10.3390/ijms252413414",
                transcript=r.id + ":" + gene,
            )
        )
    (root / "dog-mito.json").write_text(
        json.dumps(
            dict(
                reference_name=r.id,
                annotation_name="refseq-genbank",
                annotation_version=r.id,
                species="canis_familiaris",
                source_url="https://www.ncbi.nlm.nih.gov/nuccore/NC_002008.4",
                source_sha256=hashlib.sha256(
                    (src / "dog-mito-NC_002008.4.gb").read_bytes()
                ).hexdigest(),
                sequence=str(r.seq),
                records=records,
                cases=cases,
            ),
            indent=2,
        )
        + "\n"
    )
    # Public curated structural placements; retain the source model and uncertainty.
    sv = []
    for name, release, species, contig, a, b, tid, source, kind in [
        (
            "Pax3-Sp2H",
            115,
            "mouse",
            "1",
            78108913,
            78108944,
            "ENSMUST00000004994",
            "https://archive.har.mrc.ac.uk/stock?id=239",
            "truncation",
        ),
        (
            "COL9A2-Samoyed",
            "2020_10",
            "dog",
            "15",
            2732685,
            2733955,
            "ENSCAFT00845017876",
            "https://omia.org/OMIA001523/9615/",
            "start_loss",
        ),
    ]:
        g = cached_release(release, species=species)
        t = g.transcript_by_id(tid)
        positions = [
            p
            for e in sorted(t.exons, key=lambda e: e.start, reverse=t.strand == "-")
            for p in (
                range(e.end, e.start - 1, -1)
                if t.strand == "-"
                else range(e.start, e.end + 1)
            )
        ]
        sequence = "".join(
            base for p, base in zip(positions, t.sequence) if not a <= p <= b
        )
        row = dict(
            id=name,
            genome=species,
            contig=contig,
            start=a,
            end=b,
            sv_type="DEL",
            transcript=tid,
            source=source,
            placement="public curated native-assembly placement; original experiment predates this assembly",
            reference_protein=t.protein_sequence,
            expected_cdna=sequence,
            expectation=kind,
        )
        if kind == "truncation":
            start = min(t.start_codon_spliced_offsets)
            coding = sequence[start:]
            row["expected_protein"] = str(
                Seq(coding[: len(coding) // 3 * 3]).translate(to_stop=True)
            )
            assert row["expected_protein"] == t.protein_sequence[:237]
        else:
            row["coordinate_caveat"] = (
                "Curated endpoints span 1271 bp; original report describes 1267 bp. Both delete the native annotated start. Exact patient breakpoint sequence is not asserted."
            )
        sv.append(row)
    (root / "structural-cases.json").write_text(json.dumps(sv, indent=2) + "\n")


build_extra_animals(root, src)

genomes = {
    "mouse": cached_release(115, species="mouse"),
    "dog": cached_release("2020_10", species="dog"),
    "chm13": load_genome("t2t-chm13", cache_root=str(args.chm13_cache)),
}
cases = [
    dict(
        id="Atp7a-A514T",
        genome="mouse",
        contig="X",
        position=105135710,
        ref="G",
        alt="A",
        transcript="ENSMUST00000055941",
        expected_hgvs="p.A514T",
        source="https://www.alliancegenome.org/allele/MGI%3A1856098",
        pmid="9147645",
        placement="curated GRCm39",
    ),
    dict(
        id="Atp7a-LA801del",
        genome="mouse",
        contig="X",
        position=105144150,
        ref="TCTTGC",
        alt="",
        transcript="ENSMUST00000055941",
        expected_hgvs="p.L801_A802del",
        source="https://www.alliancegenome.org/allele/MGI%3A1856098",
        pmid="9147645",
        placement="curated GRCm39",
    ),
    dict(
        id="mt-Nd6-m1Dwa",
        genome="mouse",
        contig="MT",
        position=13890,
        ref="C",
        alt="CC",
        transcript="ENSMUST00000082419",
        source="https://www.informatics.jax.org/allele/MGI%3A3777241",
        pmid="18276892",
        placement="curated GRCm39",
        published_stop_position=79,
        reference_stop_position=80,
        sequence_caveat="Published LA9/C3H termination at codon 79 differs from the 80th-codon stop independently calculated on pinned GRCm39; no strain-equivalence claim",
    ),
    dict(
        id="Ap2s1-del17",
        genome="mouse",
        contig="7",
        position=16477295,
        ref="T",
        alt="C",
        transcript="ENSMUST00000086112",
        source="https://pmc.ncbi.nlm.nih.gov/articles/PMC5824975/",
        pmid="29479578",
        placement="curated GRCm39",
        observed_deleted_amino_acids=[35, 51],
    ),
    dict(
        id="GRIK2-L197F",
        genome="dog",
        contig="12",
        position=60469829,
        ref="C",
        alt="T",
        transcript="ENSCAFT00845011819",
        expected_hgvs="p.L197F",
        source="https://doi.org/10.3390/ani13050810",
        placement="native ROS_Cfam Sanger comparison",
    ),
    dict(
        id="VPS54-D17E",
        genome="dog",
        contig="10",
        position=64662153,
        ref="A",
        alt="T",
        transcript="ENSCAFT00845035084",
        expected_hgvs="p.D17E",
        source="https://doi.org/10.3390/ani13050810",
        placement="published ROS_Cfam placement of CanFam3.1 WES call",
    ),
    dict(
        id="CCDC85A-P230L",
        genome="dog",
        contig="10",
        position=57941870,
        ref="C",
        alt="T",
        transcript="ENSCAFT00845021863",
        expected_hgvs="p.P230L",
        source="https://doi.org/10.3390/ani13050810",
        placement="published ROS_Cfam placement of CanFam3.1 WES call",
    ),
    dict(
        id="BRINP3-R193P",
        genome="dog",
        contig="38",
        position=8412424,
        ref="G",
        alt="C",
        transcript="ENSCAFT00845045693",
        expected_hgvs="p.R193P",
        source="https://doi.org/10.3390/ani13050810",
        placement="published ROS_Cfam placement of CanFam3.1 WES call",
    ),
    dict(
        id="FOXI3-dup7",
        genome="dog",
        contig="17",
        position=38764881,
        ref="C",
        alt="CGCCCGCC",
        transcript="ENSCAFT00845023151",
        source="https://omia.org/OMIA000323/9615/",
        pmid="18787161",
        placement="OMIA curated ROS_Cfam placement",
        refseq_protein="NP_001129118.1",
        refseq_hgvs="p.A23Rfs*219",
        ensembl_n_terminal_extension=152,
    ),
]
selected = {
    "SLC5A9": "ENST05220073379",
    "ERAP1": "ENST05220213652",
    "FBXW12": "ENST05220211096",
    "TCP11": "ENST05220061616",
    "PKD1L3": "ENST05220204757",
    "ZNF527": "ENST05220042385",
    "TRIM51G": "ENST05220063847",
    "DHRS4L2": "ENST05220112563",
    "ABCA10": "ENST05220229140",
}
s5 = json.loads((src / "swedish-s5.json").read_text())
seen = set()
for row in s5:
    if row["validated"] and row["gene"] in selected and row["gene"] not in seen:
        seen.add(row["gene"])
        cases.append(
            dict(
                row,
                id="SweGen-" + row["gene"],
                genome="chm13",
                contig=row["contig"][3:],
                transcript=selected[row["gene"]],
                source="https://doi.org/10.1101/gr.279320.124",
                placement="native CHM13 call; PacBio CLR validated",
                published_consequence="SnpEff RefSeq v5.1 predicted loss of function",
            )
        )
for line in (src / "chm13-new-calls.vcf").read_text().splitlines():
    f = line.split("\t")
    if f[1] not in ("21798486", "21801591"):
        continue
    for tid, name in [("ENST05220004619", "RIMBP3B"), ("ENST05220004799", "RIMBP3C")]:
        cases.append(
            dict(
                id="1KGP-" + name + "-" + f[1],
                genome="chm13",
                contig=f[0][3:],
                position=int(f[1]),
                ref=f[3],
                alt=f[4],
                transcript=tid,
                source="https://github.com/marbl/CHM13",
                placement="native 1000 Genomes CHM13 call; population frequency file; not individually validated",
                source_vcf_row=line,
                source_vcf_row_sha256=hashlib.sha256(line.encode()).hexdigest(),
            )
        )

# Read every native ROS_Cfam tumour-table row, including unmatchable ones.
# Match the paper's amino acid and position against independently edited RNA.
import openpyxl
from Bio.SeqUtils import seq1


def edited_cdna(t, position, ref, alt):
    positions = [
        p
        for e in sorted(t.exons, key=lambda e: e.start, reverse=t.strand == "-")
        for p in (
            range(e.end, e.start - 1, -1)
            if t.strand == "-"
            else range(e.start, e.end + 1)
        )
    ]
    if not all(p in positions for p in range(position, position + len(ref))):
        return None
    offsets = [positions.index(p) for p in range(position, position + len(ref))]
    lo, hi = min(offsets), max(offsets) + 1
    if t.strand == "-":
        ref, alt = (
            str(Seq(ref).reverse_complement()),
            str(Seq(alt).reverse_complement()),
        )
    if t.sequence[lo:hi] != ref:
        return None
    return t.sequence[:lo] + alt + t.sequence[hi:]


audit = []
for filename in ("dog-carcinoma.xlsx", "dog-sarcoma.xlsx"):
    workbook = openpyxl.load_workbook(src / filename, read_only=True, data_only=True)
    for number, values in enumerate(workbook.active.values, 1):
        if number < 3 or not values[0] or not isinstance(values[1], (float, int)):
            continue
        row = dict(
            source_file=filename,
            source_row=number,
            source_sha256=sha256(src / filename),
            contig=str(values[0]),
            position=int(values[1]),
            ref=values[2],
            alt=values[3],
            gene_id=values[8],
            gene=values[9] if filename == "dog-carcinoma.xlsx" else values[10],
            consequence=values[11],
            cdna_position=values[12],
            protein_position=values[13],
            amino_acids=values[14],
            filter=values[5],
            genotypes=list(values[6:8]),
            source="https://doi.org/10.1371/journal.pone.0307792",
        )
        transcripts = [
            t
            for t in genomes["dog"].transcripts_at_locus(row["contig"], row["position"])
            if t.gene_id == row["gene_id"]
        ]
        row["transcript_ids"] = [t.id for t in transcripts]
        row["independent_matches"] = []
        amino = row["amino_acids"]
        pos = row["protein_position"]
        if isinstance(amino, str) and "/" in amino and isinstance(pos, (int, float)):
            before, after = [seq1(a) if a != "*" else "*" for a in amino.split("/")]
            offset = int(pos) - 1
            for t in transcripts:
                if not t.complete or t.protein_sequence[offset : offset + 1] != before:
                    continue
                cdna = edited_cdna(t, row["position"], row["ref"], row["alt"])
                if cdna is None:
                    continue
                coding = cdna[min(t.start_codon_spliced_offsets) :]
                protein = str(
                    Seq(coding[: len(coding) // 3 * 3]).translate(to_stop=True)
                )
                wanted = (
                    t.protein_sequence[:offset]
                    + after
                    + t.protein_sequence[offset + 1 :]
                ).split("*")[0]
                if protein == wanted:
                    row["independent_matches"].append(t.id)
        if row["independent_matches"]:
            row["status"] = (
                "published amino acid and position reproduced from native RNA"
            )
            cases.append(
                dict(
                    row,
                    id="ROS-tumour-"
                    + filename.split("-")[1].split(".")[0]
                    + "-"
                    + str(number),
                    genome="dog",
                    transcript=sorted(row["independent_matches"])[0],
                    placement="native ROS_Cfam WGS; published prediction, not protein observation",
                )
            )
        elif "splice_" in row["consequence"]:
            row["status"] = (
                "splice annotation; no observed mutant RNA or protein supplied"
            )
            complete = sorted(t.id for t in transcripts if t.complete)
            if complete:
                cases.append(
                    dict(
                        row,
                        id="ROS-splice-"
                        + filename.split("-")[1].split(".")[0]
                        + "-"
                        + str(number),
                        genome="dog",
                        transcript=complete[0],
                        placement="native ROS_Cfam WGS; published splice prediction, not observed splicing",
                        protein_expectation="conditional normal splicing",
                    )
                )
        elif not transcripts:
            row["status"] = (
                "no overlapping transcript from the published gene in pinned annotation"
            )
        elif not any(t.complete for t in transcripts):
            row["status"] = "only incomplete native transcript models"
        elif not isinstance(amino, str) or "/" not in amino:
            row["status"] = (
                "published amino-acid field has no unambiguous substitution delimiter"
            )
        else:
            row["status"] = (
                "published amino-acid position or substitution differs from pinned native models"
            )
        audit.append(row)
    workbook.close()
(root / "dog-tumour-source-audit.json").write_text(json.dumps(audit, indent=2) + "\n")

# Expected edit products are computed from native RNA and an independent BioPython
# genetic-code implementation. No Varcode import or effect prediction is used.
for case in cases:
    g = genomes[case["genome"]]
    t = g.transcript_by_id(case["transcript"])
    case["strand"] = t.strand
    case["reference_protein"] = t.protein_sequence
    if "observed_deleted_amino_acids" in case:
        a, b = case["observed_deleted_amino_acids"]
        case["expected_protein"] = t.protein_sequence[: a - 1] + t.protein_sequence[b:]
        continue
    positions = [
        p
        for e in sorted(t.exons, key=lambda e: e.start, reverse=t.strand == "-")
        for p in (
            range(e.end, e.start - 1, -1)
            if t.strand == "-"
            else range(e.start, e.end + 1)
        )
    ]
    p, r, a = case["position"], case["ref"], case["alt"]
    if not all(x in positions for x in range(p, p + len(r))) or (
        a.startswith(r)
        and len(a) > len(r)
        and not any(e.start <= p + len(r) - 1 < p + len(r) <= e.end for e in t.exons)
    ):
        print("INTRONIC SOURCE", case["id"], flush=True)
        case["expected_cdna"] = t.sequence
        case["expected_protein"] = t.protein_sequence
        case["protein_expectation"] = "conditional normal splicing"
        continue
    offsets = [positions.index(x) for x in range(p, p + len(r))]
    lo, hi = min(offsets), max(offsets) + 1
    if t.strand == "-":
        r, a = str(Seq(r).reverse_complement()), str(Seq(a).reverse_complement())
    assert t.sequence[lo:hi] == r, (case["id"], r, t.sequence[lo:hi])
    sequence = t.sequence[:lo] + a + t.sequence[hi:]
    case["expected_cdna"] = sequence
    start = min(t.start_codon_spliced_offsets)
    code = 2 if t.contig == "MT" else 1
    coding = sequence[start:]
    coding = coding[: len(coding) // 3 * 3]
    case["expected_protein"] = str(Seq(coding).translate(table=code, to_stop=True))

for name, g in genomes.items():
    structural = json.loads((root / "structural-cases.json").read_text())
    tids = {c["transcript"] for c in cases + structural if c["genome"] == name}
    genes = {g.transcript_by_id(t).gene_id for t in tids}
    # EnsemblRelease constructor dictionaries name the release, so use installed public paths.
    g.index()
    source_g = g.gtf_path
    source_cdna = g.transcript_fasta_paths[0]
    source_pep = g.protein_fasta_paths[0]

    def lines(path):
        with gzip.open(path, "rt") if compressed(path) else open(path) as f:
            yield from f

    out = []
    for line in lines(source_g):
        if line.startswith("#"):
            out.append(line)
            continue
        tid = re.search(r'transcript_id "([^"]+)"', line)
        gene = re.search(r'gene_id "([^"]+)"', line)
        if (tid and tid[1] in tids) or (
            "\tgene\t" in line and gene and gene[1] in genes
        ):
            out.append(line)
    (root / (name + ".gtf.gz")).write_bytes(
        gzip.compress("".join(out).encode(), mtime=0)
    )
    from Bio import SeqIO

    for role, path, ids in [
        ("cdna", source_cdna, tids),
        ("pep", source_pep, {g.transcript_by_id(t).protein_id for t in tids}),
    ]:
        output = []
        found = set()
        with gzip.open(path, "rt") if compressed(path) else open(path) as f:
            for rec in SeqIO.parse(f, "fasta"):
                if rec.id.split(".")[0] in ids:
                    output.append(">" + rec.description + "\n" + str(rec.seq) + "\n")
                    found.add(rec.id.split(".")[0])
        assert found == ids, (name, role, ids - found)
        (root / (name + "." + role + ".fa.gz")).write_bytes(
            gzip.compress("".join(output).encode(), mtime=0)
        )
    meta = dict(
        reference_name=g.reference_name,
        annotation_name=g.annotation_name,
        annotation_version=str(g.annotation_version),
        species=g.species.latin_name,
        assembly_accession=g.assembly_accession,
        reference_metadata=g.reference_metadata,
    )
    meta["sources"] = [
        dict(role=role, path=Path(path).name, sha256=sha256(path))
        for role, path in [
            ("gtf", source_g),
            ("cdna", source_cdna),
            ("pep", source_pep),
        ]
    ]
    (root / (name + ".json")).write_text(json.dumps(meta, indent=2) + "\n")
windows = json.loads((src / "dog-splice-dna.json").read_text())
g = genomes["chm13"]
g = type(g)(**dict(g.to_dict(), genome_fasta_path_or_url=str(args.chm13_dna)))
for tid in sorted({c["transcript"] for c in cases if c["genome"] == "chm13"}):
    t = g.transcript_by_id(tid)
    a, b = max(1, t.start - 200), t.end + 200
    windows.append(
        dict(
            genome="chm13",
            contig=t.contig,
            start=a,
            end=b,
            source_contig=args.dna_contig_prefix + t.contig,
            sequence=g.sequence(args.dna_contig_prefix + t.contig, a, b),
            transcript=tid,
        )
    )
windows.append(
    dict(
        genome="mouse",
        contig="7",
        start=16472335,
        end=16483219,
        sequence="".join((src / "mouse-ap2s1-dna.fa").read_text().splitlines()[1:]),
    )
)
(root / "dna-windows.json.gz").write_bytes(
    gzip.compress((json.dumps(windows, indent=2) + "\n").encode(), mtime=0)
)
for c in cases:
    if c["genome"] != "chm13":
        continue
    for role, file in [
        ("segmental_duplication", "chm13v2.0_SD.bed"),
        ("newly_resolved", "chm13v2-unique_to_hg38.bed"),
    ]:
        c[role] = []
        for line in (src / file).read_text().splitlines():
            f = line.split()
            if f[0] == "chr" + c["contig"] and int(f[1]) <= c["position"] - 1 < int(
                f[2]
            ):
                c[role].append(line)
(root / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
(root / "FOXI3-NP_001129118.1.fa").write_bytes(
    (src / "dog-FOXI3-NP_001129118.1.fa").read_bytes()
)
print(
    "Built", len(cases), "cases", [(p.name, p.stat().st_size) for p in root.iterdir()]
)
