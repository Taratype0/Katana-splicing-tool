import base64
import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from katana_validate_inputs import InputError, validate_samples, validate_contrasts, write_sheet, SAMPLE_COLUMNS
from kanata_stage_fastq import StagingError, stage_reads


class LaneTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def rows(self, count=6):
        rows = []
        for index in range(count):
            for lane in ['L2', 'L1']:
                row = dict(sample=f'S{index}', condition='WT' if index % 2 == 0 else 'MUT',
                           subject=f'mouse{index // 2}', lane=lane, read_length='4', strandedness='auto')
                for mate in (1, 2):
                    path = self.root / f'S{index}_{lane}_{mate}.fq.gz'
                    with gzip.open(path, 'wb') as handle:
                        handle.write(f'@S{index}_{lane}/{mate}\nACGT\n+\nIIII\n'.encode())
                    row[f'fastq_{mate}'] = str(path)
                rows.append(row)
        return rows

    def validate(self, rows):
        return validate_samples(rows, self.root, check=True, read_length=None)

    def test_six_paired_animals_stay_six_samples(self):
        raw = self.rows()
        samples = self.validate(raw)
        self.assertEqual(len(samples), 6)
        _, audit = validate_contrasts([dict(contrast='MUT_vs_WT', treatment='MUT', control='WT', paired='true')], samples,
                                     design='~ subject + condition', analysis_purpose='analysis')
        self.assertEqual(audit[0]['pair_ids'], ['mouse0', 'mouse1', 'mouse2'])
        target = self.root / 'valid.csv'
        write_sheet(target, SAMPLE_COLUMNS, list(raw[0]), samples)
        with target.open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 6)

    def test_merges_both_mates_in_identical_deterministic_order(self):
        sample = self.validate(self.rows(1))[0]
        lanes = json.loads(base64.b64decode(sample['fastq_lanes']))
        out = self.root / 'merged'
        audit = stage_reads(dict(reads=[p for pair in lanes for p in pair], single_end=False, lane_count=2), sample='S0', output_dir=out)
        self.assertEqual(audit['lane_count'], 2)
        self.assertEqual(len(audit['reads']), 4)
        for mate in (1, 2):
            expected = b''.join(Path(pair[mate - 1]).read_bytes() for pair in lanes)
            self.assertEqual((out / f'S0_R{mate}.fastq.gz').read_bytes(), expected)
            with gzip.open(out / f'S0_R{mate}.fastq.gz', 'rt') as stream:
                self.assertEqual(stream.read().splitlines()[::4], [f'@S0_L1/{mate}', f'@S0_L2/{mate}'])

    def test_conflicting_metadata(self):
        for field, value in [('subject', 'different'), ('condition', 'different'), ('strandedness', 'reverse')]:
            rows = self.rows(1)
            rows[1][field] = value
            with self.assertRaisesRegex(InputError, 'identical biological'):
                self.validate(rows)

    def test_duplicate_or_missing_lane_ids(self):
        for value in ['', 'L2']:
            rows = self.rows(1)
            rows[1]['lane'] = value
            with self.assertRaises(InputError):
                self.validate(rows)

    def test_reused_source_file_rejected(self):
        rows = self.rows(1)
        rows[1]['fastq_1'] = rows[0]['fastq_1']
        with self.assertRaisesRegex(InputError, 'also assigned'):
            self.validate(rows)

    def test_bad_second_lane_pair_no_partial_outputs(self):
        rows = self.rows(1)
        with gzip.open(rows[1]['fastq_2'], 'wb') as stream:
            stream.write(b'@WRONG/2\nACGT\n+\nIIII\n')
        paths = [r[f'fastq_{m}'] for r in rows for m in (1, 2)]
        out = self.root / 'bad'
        with self.assertRaisesRegex(StagingError, 'Mate identifiers'):
            stage_reads(dict(reads=paths, single_end=False, lane_count=2), sample='S0', output_dir=out)
        self.assertFalse(out.exists())

    def test_single_end_multiple_lanes(self):
        rows = self.rows(1)
        for row in rows:
            row['fastq_2'] = ''
        sample = self.validate(rows)[0]
        paths = [lane[0] for lane in json.loads(base64.b64decode(sample['fastq_lanes']))]
        audit = stage_reads(dict(reads=paths, single_end=True, lane_count=2), sample='S0', output_dir=self.root / 'single')
        self.assertEqual(len(audit['reads']), 2)
        self.assertFalse((self.root / 'single/S0_R2.fastq.gz').exists())

if __name__ == '__main__':
    unittest.main()
