"""Reference and ORF options reach annotation and preserve exported identity."""

import pandas as pd
import pytest

from varcode import EffectCollection
from varcode.cli.effects_script import arg_parser, main
from varcode.cli.variant_args import variant_collection_from_args


def test_cli_fasta_is_attached_without_downloading(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">99\nACGT\n")
    args = arg_parser.parse_args([
        "--genome", "GRCh38:95", "--genome-fasta", str(fasta),
        "--variant", "99", "2", "C", "A"])
    variants = variant_collection_from_args(args)
    try:
        assert variants[0].genome.reference_range("99", 1, 4) == "ACGT"
    finally:
        variants[0].genome.close()


@pytest.mark.parametrize("options,message", [
    (["--genome-fasta", "/missing/reference.fa"], "--genome-fasta"),
    (["--extra-orfs", "not-an-evidence-name"], "evidence"),
    (["--extra-orfs", "all", "--extra-orfs", "ms"], "cannot be combined"),
])
def test_bad_reference_or_orf_option_fails_before_annotating(options, message, monkeypatch, capsys):
    from varcode import VariantCollection

    def unexpected(*args, **kwargs):
        pytest.fail("Bad options reached annotation")

    monkeypatch.setattr(VariantCollection, "effects", unexpected)
    with pytest.raises(SystemExit) as error:
        main(["--genome", "GRCh38:95", "--variant", "7", "66205486", "C", "G", *options])
    assert error.value.code == 1
    assert message in capsys.readouterr().err


def test_cli_orf_selection_and_csv_round_trip(tmp_path):
    pytest.importorskip("orfdata")
    path = tmp_path / "effects.csv"
    main(["--genome", "GRCh38:95", "--variant", "7", "66205486", "C", "G",
          "--extra-orfs", "ms", "--output-csv", str(path)])
    frame = pd.read_csv(path)
    assert "c7riboseqorf56" in set(frame.orf_id.dropna())
    restored = EffectCollection.from_csv(path)
    annotation = next(a for a in restored.orf_annotations if a.orf_id == "c7riboseqorf56")
    assert annotation.effect.short_description == "p.P2A"
