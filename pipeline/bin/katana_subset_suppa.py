#!/usr/bin/env python3
"""Remove explicitly excluded sample columns without altering SUPPA values."""
import argparse
import csv


def subset(source, destination, excluded):
    with open(source, newline='') as src, open(destination, 'w', newline='') as dst:
        reader = csv.reader(src, delimiter='\t')
        writer = csv.writer(dst, delimiter='\t', lineterminator='\n')
        header = next(reader)
        first = next(reader, None)
        # R write.table can omit the row-ID header; SUPPA also accepts a named ID.
        offset = 1 if first is not None and len(first) == len(header) + 1 else 0
        if first is not None and len(first) != len(header) + offset:
            raise ValueError('Invalid SUPPA matrix width')
        remove = {i + offset for i, name in enumerate(header) if name in excluded}
        writer.writerow([name for name in header if name not in excluded])
        def emit(row):
            if len(row) != len(header) + offset:
                raise ValueError('Inconsistent SUPPA matrix width')
            writer.writerow([value for i, value in enumerate(row) if i not in remove])
        if first is not None:
            emit(first)
        for row in reader:
            emit(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exclude', default='')
    parser.add_argument('paths', nargs=8)
    args = parser.parse_args()
    excluded = {s.strip() for s in args.exclude.split(';') if s.strip()}
    for i in range(0, 8, 2):
        subset(args.paths[i], args.paths[i + 1], excluded)
