"""Malformed source metadata must not become a plausible SV definition."""

from copy import deepcopy

import pytest
import numpy as np

from varcode import load_vcf, parse_symbolic_alt


@pytest.mark.parametrize('alt', ['<DEL>', 'N[2:300[', 'N.'])
@pytest.mark.parametrize('key,value', [('SVTYPE', True), ('SVTYPE', 1),
    ('CIPOS', [0]), ('CIPOS', [0, 1, 2]), ('CIPOS', [True, 3]),
    ('CIPOS', [1.5, 3]), ('CIPOS', [3, -2]), ('CIEND', [-2, None])])
def test_consumed_metadata_errors_are_actionable(alt, key, value):
    with pytest.raises(ValueError, match='INFO/' + key):
        parse_symbolic_alt('1', 100, 'N', alt, {'END': 200, key: value}, genome='GRCh38')


@pytest.mark.parametrize('value', [True, 200.5, '200.5', [200, 300]])
def test_symbolic_end_requires_one_integer(value):
    with pytest.raises(ValueError, match='INFO/END'):
        parse_symbolic_alt('1', 100, 'N', '<DEL>', {'END': value}, genome='GRCh38')


@pytest.mark.parametrize('alt', ['<DEL>', 'N[2:300[', '.N'])
def test_confidence_intervals_survive_without_mutating_input(alt):
    info = dict(END=['200'], SVTYPE=['DEL'], CIPOS=[-10, 5], CIEND='-4,8')
    original = deepcopy(info)
    variant = parse_symbolic_alt('1', 100, 'N', alt, info, genome='GRCh38')
    assert variant.ci_start == (-10, 5)
    assert variant.ci_end == (-4, 8)
    assert info == original


def test_missing_intervals_and_unconsumed_fields_do_not_invent_precision():
    variant = parse_symbolic_alt('1', 100, 'N', 'N[2:300[',
        dict(END='unused malformed symbolic endpoint', CIPOS=[None, None], CIEND='.'), genome='GRCh38')
    assert variant.ci_start is None and variant.ci_end is None
    assert parse_symbolic_alt('1', 100, 'A', 'T', dict(SVTYPE=True), genome='GRCh38') is None


@pytest.mark.parametrize('missing', [None, '.', [None], ['.'], [None, None]])
def test_parser_missing_interval_representations_remain_unknown(missing):
    variant = parse_symbolic_alt('1', 100, 'N', '<DEL>',
        dict(END=np.int64(200), CIPOS=missing, CIEND=[np.int64(-2), np.int64(2)]), genome='GRCh38')
    assert variant.end == 200 and variant.ci_start is None and variant.ci_end == (-2, 2)


def test_vcf_loader_selects_alt_specific_intervals_and_retains_original_metadata(tmp_path):
    path = tmp_path / 'multi.vcf'
    path.write_text('##fileformat=VCFv4.5\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Offsets">\n'
        '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Offsets">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        '1\t100\tevent\tN\t<DEL>,<INV>\t.\tPASS\tEND=200;CIPOS=-5,5,-8,8;CIEND=-2,2,-3,3\n')
    collection = load_vcf(str(path), genome='GRCh38', parse_structural_variants=True)
    assert sorted((v.end, v.ci_start, v.ci_end) for v in collection) == [
        (200, (-8, 8), (-3, 3)), (200, (-5, 5), (-2, 2))]
    for variant in collection:
        assert collection.metadata[variant]['info']['CIPOS'] == [-5, 5, -8, 8]


def test_vcf_loader_does_not_hide_malformed_interval(tmp_path):
    path = tmp_path / 'invalid.vcf'
    path.write_text('##fileformat=VCFv4.2\n'
        '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Offsets">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        '1\t100\tbad\tN\tN[2:300[\t.\tPASS\tCIPOS=0\n')
    with pytest.raises(ValueError, match='INFO/CIPOS'):
        load_vcf(str(path), genome='GRCh38', parse_structural_variants=True)
