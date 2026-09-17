"""Keep the reader-facing examples executable, using their actual Markdown."""

from pathlib import Path
import re
import shutil

import pytest
from pyensembl import cached_release

import varcode
from varcode import EffectCandidate, StructuralVariant, Variant
from varcode.effect_hypotheses import RealizedEffectCandidate


ROOT = Path(__file__).resolve().parents[1]


def _blocks(path, heading=None):
    text = (ROOT / path).read_text()
    if heading is not None:
        text = text.split(heading + "\n", 1)[1]
    return re.findall(r"```python\n(.*?)```", text, re.S)


def _run(blocks, namespace=None):
    namespace = {} if namespace is None else namespace
    for block in blocks:
        exec(compile(block, "documentation example", "exec"), namespace)
    return namespace


@pytest.fixture
def input_directory(tmp_path, monkeypatch):
    for fixture in (ROOT / "tests/data/documentation").glob("*.vcf"):
        shutil.copy(fixture, tmp_path / fixture.name)
    shutil.copy(tmp_path / "variants.vcf", tmp_path / "tumor.vcf")
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.mark.parametrize("path,count", [
    ("README.md", 2),
    ("docs/getting_started.md", 4),
    ("docs/index.md", 1),
    ("docs/effect_annotation.md", 3),
    ("docs/structural_variants.md", 1),
    ("docs/transforms.md", 1),
    ("docs/germline.md", 1),
])
def test_documentation_first_steps(path, count, input_directory):
    namespace = _run(_blocks(path)[:count])
    if path == "docs/transforms.md":
        assert len(namespace["paired"]) == len(namespace["normalized"]) == 1
    else:
        assert len(namespace["effects"]) > 0
    if path in ("README.md", "docs/getting_started.md"):
        assert namespace["protein"][158] == "M"
    if path == "docs/getting_started.md":
        assert (input_directory / "effects.csv").is_file()
    if path == "docs/structural_variants.md":
        assert namespace["variants"][0].is_structural


def test_documentation_germline_scenarios():
    namespace = _run(_blocks(
        "docs/germline.md", "## Concrete example: same codon, three scenarios")[:4])
    assert namespace["eff_cis"].short_description == "p.S159T"
    assert namespace["eff_trans"].short_description == "p.L159M"
    assert {c.effect.short_description for c in namespace["eff"].candidates} == {
        "p.S159T", "p.L159M"}

    _run(_blocks("docs/germline.md", "## Loss of heterozygosity (LOH)")[:1],
         namespace)
    assert any(getattr(effect, "is_loh", False)
               for effect in namespace["overlap"].effects(germline=namespace["ctx"]))

    validation = _blocks("docs/germline.md", "## Cross-VCF build mismatch")[:1]
    _run(validation, namespace)
    namespace["ctx"] = varcode.GermlineContext.from_variants(
        [namespace["germline"]], reference_name="GRCh37")
    with pytest.raises(varcode.GenomeBuildMismatchError):
        _run(validation, namespace)


def test_documentation_splice_accessors():
    text = (ROOT / "docs/effect_annotation.md").read_text()
    section = text.split("### The `SpliceOutcomeSet` shape\n", 1)[1]
    section = section.split("### RNA evidence reconciliation\n", 1)[0]
    blocks = re.findall(r"```python\n(.*?)```", section, re.S)
    namespace = _run(blocks[:1])
    namespace["effect"] = namespace["splice_set"]
    _run(blocks[1:], namespace)
    assert isinstance(namespace["splice_set"], varcode.SpliceOutcomeSet)
    assert namespace["coding"] is not None


@pytest.mark.parametrize("kind", ["point", "BND", "CNV"])
def test_documentation_experimental_result_access(kind):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    if kind == "point":
        variant = Variant("7", 117531100, "T", "A", genome)
    elif kind == "BND":
        variant = StructuralVariant(
            "7", transcript.start + 50, "BND", alt="N]22:15500000]",
            mate_contig="22", mate_start=15500000, genome=genome)
    else:
        variant = StructuralVariant(
            "7", transcript.start + 50, "CNV", end=transcript.start + 100,
            genome=genome)
    namespace = dict(varcode=varcode, variant=variant, transcript=transcript)
    advanced_blocks = _blocks("docs/effect_annotation.md",
                              "## Advanced: annotators and implementation limits")
    assert "comparison =" in advanced_blocks[0]
    assert "candidate.outcomes" in advanced_blocks[1]
    _run(advanced_blocks[:2], namespace)
    _run(_blocks("docs/effect_annotation.md", "## Read an effect")[:1], namespace)
    candidates = getattr(namespace["experimental"], "candidates", ())
    if kind == "point":
        assert candidates and all(isinstance(c, RealizedEffectCandidate)
                                  for c in candidates)
        assert namespace["protein"][158] == "M"
    elif kind == "BND":
        assert candidates and all(isinstance(c, EffectCandidate) for c in candidates)
    else:
        assert not candidates
        assert isinstance(namespace["experimental"], varcode.Unresolved)
