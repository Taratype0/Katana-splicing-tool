#!/usr/bin/env python3
"""Gather per-comparison DESeq2 metadata without changing scientific tables."""
import argparse
import csv
from pathlib import Path


def merge(contrastsheet, model_files, metadata_files):
    with Path(contrastsheet).open(newline='', encoding='utf-8-sig') as handle:
        comparisons = list(csv.DictReader(handle))
    expected = [row['contrast'] for row in comparisons]
    if not expected or len(set(expected)) != len(expected):
        raise ValueError('Expected unique contrasts')
    records = {}
    fields = None
    for path in model_files:
        with Path(path).open(newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if fields is None:
                fields = reader.fieldnames
            if reader.fieldnames != fields:
                raise ValueError('Model schemas differ')
            rows = list(reader)
        if len(rows) != 1 or rows[0]['contrast'] in records:
            raise ValueError('Each task must report exactly one unique model')
        records[rows[0]['contrast']] = rows[0]
    if set(records) != set(expected):
        raise ValueError('Missing or unexpected contrast models')
    for row in comparisons:
        model = records[row['contrast']]
        if any(model[key] != row[key] for key in ['treatment', 'control']):
            raise ValueError('Model direction differs from declared contrast')
    metadata = [Path(path).read_bytes() for path in metadata_files]
    if len(metadata) != len(expected) or any(data != metadata[0] for data in metadata):
        raise ValueError('Missing or inconsistent differential metadata')
    with Path('models.tsv').open('x', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records[name] for name in expected)
    Path('differential_metadata.tsv').write_bytes(metadata[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contrasts', required=True)
    parser.add_argument('--models', nargs='+', required=True)
    parser.add_argument('--metadata', nargs='+', required=True)
    args = parser.parse_args()
    merge(args.contrasts, args.models, args.metadata)
