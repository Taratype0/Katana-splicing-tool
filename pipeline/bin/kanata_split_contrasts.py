#!/usr/bin/env python3
"""Write one CSV row per comparison, preserving field values and column order."""
import argparse
import csv
import re
from pathlib import Path


def split_contrasts(source, destination):
    with Path(source).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or not {'contrast', 'treatment', 'control'}.issubset(fields):
            raise ValueError('Contrasts require contrast,treatment,control columns')
        if len(fields) != len(set(fields)):
            raise ValueError('Duplicate contrast columns')
        rows = list(reader)
    if not rows:
        raise ValueError('No contrasts declared')
    seen = set()
    for row in rows:
        name = row['contrast']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name or ''):
            raise ValueError(f'Unsafe contrast identifier: {name!r}')
        if name in seen or None in row or any(value is None for value in row.values()):
            raise ValueError(f'Duplicate contrast or malformed row: {name!r}')
        if not row['treatment'] or not row['control'] or row['treatment'] == row['control']:
            raise ValueError(f'Invalid treatment/control: {name!r}')
        seen.add(name)
    output = Path(destination)
    output.mkdir(parents=True, exist_ok=False)
    for row in rows:
        with (output / f"contrast.{row['contrast']}.csv").open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
            writer.writeheader()
            writer.writerow(row)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('destination')
    args = parser.parse_args()
    split_contrasts(args.source, args.destination)
