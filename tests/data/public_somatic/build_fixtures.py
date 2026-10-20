"""Rebuild the public somatic fixtures from checksum-verified downloads.

Download manifest.json's source_url entries as their source_file names first.
This script is offline; pandas and openpyxl are required only for regeneration.
"""

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import shutil


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def seqc2(source, destination):
    filters = Counter()
    selected = 0
    total = 0
    with source.open('rb') as raw, gzip.open(raw, 'rt') as infile, destination.open('wb') as output:
        with gzip.GzipFile(filename='', fileobj=output, mode='wb', mtime=0) as compressed:
            def write(line):
                compressed.write((line + '\n').encode())
            write('##fileformat=VCFv4.2')
            write('##reference=GRCh38')
            write('##source=SEQC2_v1.2_complete_PASS_truth')
            write('##INFO=<ID=SOURCE_ROW,Number=1,Type=Integer,Description="One-based variant record in the pinned source VCF">')
            write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
            for line in infile:
                if line.startswith('#'):
                    continue
                total += 1
                fields = line.rstrip('\n').split('\t')
                filters[fields[6]] += 1
                labels = set(fields[6].split(';'))
                if 'PASS' not in labels or not labels.intersection({'HighConf', 'MedConf'}):
                    continue
                assert 'ArmLossInNormal' not in fields[7]
                write('\t'.join(fields[:7] + ['SOURCE_ROW=%d' % total]))
                selected += 1
    return dict(source_sha256=sha256(source), sha256=sha256(destination),
                records=selected, source_records=total, source_filters=dict(sorted(filters.items())))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_directory', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name('manifest.json').read_text())
    args.output_directory.mkdir(parents=True, exist_ok=True)
    result = {}
    for name, spec in manifest['fixtures'].items():
        source = args.source_directory / spec['source_file']
        destination = args.output_directory / name
        assert sha256(source) == spec['source_sha256'], source
        if name.startswith('hcc1395-seqc2-'):
            result[name] = seqc2(source, destination)
            assert result[name]['records'] == spec['records']
        elif name.endswith('.tsv'):
            import pandas as pd
            frame = pd.read_excel(
                source, sheet_name=spec['table'],
                header=2 if name.startswith('canine') else 0,
                keep_default_na=False)
            frame.columns = frame.columns.str.strip()
            if name.startswith('hcc1395'):
                # Remove exactly the source's final prose footnote, not any
                # data row that happens to have an empty or unexpected field.
                last = frame.iloc[-1]
                assert str(last.iloc[0]).startswith('Table S4:')
                assert all(value == '' for value in last.iloc[1:])
                frame = frame.iloc[:-1]
            assert len(frame) == spec['records']
            frame.to_csv(destination, sep='\t', index=False, lineterminator='\n')
        else:
            shutil.copyfile(source, destination)
        assert sha256(destination) == spec['sha256'], destination
    print(json.dumps(result, indent=2))
