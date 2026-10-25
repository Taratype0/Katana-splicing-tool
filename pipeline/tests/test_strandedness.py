from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / 'bin'
sys.path.insert(0, str(BIN))
from katana_resolve_strandedness import infer, resolve, thresholds
from katana_validate_inputs import InputError, validate_samples, validate_contrasts


def native_text(forward=0.95, reverse=0.04, single=False):
    layout = 'SingleEnd' if single else 'PairEnd'
    rules = ('++,--', '+-,-+') if single else ('1++,1--,2+-,2-+', '1+-,1-+,2++,2--')
    return (f'This is {layout} Data\nFraction of reads failed to determine: {1-forward-reverse:.4f}\n'
            f'Fraction of reads explained by "{rules[0]}": {forward:.4f}\n'
            f'Fraction of reads explained by "{rules[1]}": {reverse:.4f}\n')


def samples(strands=None, source='fastq', single=False):
    strands = strands or ['auto'] * 4
    result=[]
    for i, strand in enumerate(strands):
        record = dict(sample=f'S{i}', condition='CONTROL' if i<2 else 'TREATED', strandedness=strand,
                      read_length='100', single_end=str(single).lower(), subject=f'animal{i}')
        if source=='fastq':
            record.update(fastq_1=f'R{i}_1.fastq.gz', fastq_2='' if single else f'R{i}_2.fastq.gz')
        else:
            record.update(genome_bam=f'S{i}.bam')
        result.append(record)
    return validate_samples(result,Path.cwd(),check=False,read_length=None)


CONTRAST=[dict(contrast='TREATED_vs_CONTROL',treatment='TREATED',control='CONTROL')]


def write_csv(path, rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


class StrandednessTests(unittest.TestCase):
    def test_native_signs_are_forward_reverse_and_unstranded_for_both_layouts(self):
        for single in (False,True):
            for forward,reverse,expected in [(0.95,0.04,'forward'),(0.04,0.95,'reverse'),(0.49,0.49,'unstranded')]:
                with self.subTest(single=single,expected=expected):
                    self.assertEqual(infer(native_text(forward,reverse,single),single_end=single)['strandedness'],expected)

    def test_normalizes_only_informative_reads_not_failed_fraction(self):
        result=infer(native_text(.09,.01),single_end=False)
        self.assertEqual(result['strandedness'],'forward')
        self.assertAlmostEqual(result['informative_forward_fraction'],.9)
        self.assertAlmostEqual(result['failed_to_determine_fraction'],.9)

    def test_inclusive_nfcore_threshold_boundaries(self):
        self.assertEqual(infer(native_text(.8,.2),single_end=False)['strandedness'],'forward')
        self.assertEqual(infer(native_text(.55,.45),single_end=False)['strandedness'],'unstranded')

    def test_ambiguous_zero_nan_malformed_and_wrong_layout_fail(self):
        for text in [native_text(.7,.3),native_text(0,0),native_text().replace('0.9500','nan'),
                     native_text()+native_text(),native_text(single=True),'Unknown Data type\n']:
            with self.subTest(text=text),self.assertRaises(InputError):
                infer(text,single_end=False)

    def test_invalid_or_overlapping_thresholds_fail(self):
        for a,b in [(.5,.1),(1.1,.1),(.8,-.1),(.8,.7),(float('nan'),.1),(.8,float('inf'))]:
            with self.subTest(thresholds=(a,b)),self.assertRaises(InputError):thresholds(a,b)

    def test_auto_accepted_for_fastq_and_bam_and_mixed_with_one_declared_type(self):
        for source in ('fastq','bam'):
            rows=samples(['auto','reverse','auto','reverse'],source)
            validate_contrasts(CONTRAST,rows,design='~ condition',analysis_purpose='analysis')
            self.assertEqual(rows[0]['strandedness'],'auto')

    def test_known_conflicts_still_fail_even_with_auto(self):
        rows=samples(['auto','reverse','auto','forward'])
        with self.assertRaisesRegex(InputError,'mixed strandedness'):
            validate_contrasts(CONTRAST,rows,design='~ condition',analysis_purpose='analysis')

    def prepare(self,root,rows,forward=.04,reverse=.95):
        sheet=root/'samples.csv';contrast=root/'contrasts.csv'
        write_csv(sheet,rows);write_csv(contrast,CONTRAST)
        evidence={}
        for row in rows:
            if row['strandedness']!='auto':continue
            report=root/f"{row['sample']}.infer_experiment.txt"
            log=root/f"{row['sample']}.infer_experiment.log"
            report.write_text(native_text(forward,reverse,row['single_end']=='true'))
            log.write_text('Total 1000 usable reads were sampled\n')
            evidence[row['sample']]={'report':str(report),'log':str(log)}
        return sheet,contrast,evidence

    def test_resolved_csv_and_audit_keep_declarations_native_evidence_and_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            sheet,contrast,evidence=self.prepare(root,samples(['auto','reverse','auto','reverse']))
            before=sheet.read_bytes()
            receipt=resolve(sheet,contrast,evidence,root/'resolved')
            self.assertEqual(sheet.read_bytes(),before)
            self.assertEqual(receipt['auto_samples'],2)
            with (root/'resolved/samplesheet.resolved.csv').open() as f:rows=list(csv.DictReader(f))
            self.assertEqual([r['sample'] for r in rows],['S0','S1','S2','S3'])
            self.assertEqual({r['strandedness'] for r in rows},{'reverse'})
            self.assertEqual(rows[0]['strandedness_requested'],'auto')
            self.assertEqual(rows[1]['strandedness_method'],'user_declared_not_inferred')
            self.assertEqual(receipt['samples'][0]['usable_reads_sampled'],1000)
            self.assertEqual(len(receipt['samples'][0]['native_report']['sha256']),64)
            self.assertEqual(json.loads((root/'resolved/strandedness.audit.json').read_text()),receipt)

    def test_actual_mixed_libraries_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            sheet,contrast,evidence=self.prepare(root,samples(['auto','reverse','auto','reverse']),.95,.04)
            with self.assertRaisesRegex(InputError,'mixed strandedness'):
                resolve(sheet,contrast,evidence,root/'out')
            self.assertFalse((root/'out').exists())

    def test_explicit_only_requires_no_inference_and_is_not_relabelled_as_measured(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            sheet,contrast,evidence=self.prepare(root,samples(['forward']*4))
            receipt=resolve(sheet,contrast,evidence,root/'out')
            self.assertEqual(receipt['auto_samples'],0)
            self.assertTrue(all(r['method']=='user_declared_not_inferred' for r in receipt['samples']))

    def test_missing_extra_or_empty_native_evidence_cannot_pass(self):
        for kind in ('missing','extra','zero','empty_log','ambiguous'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as temporary:
                root=Path(temporary)
                sheet,contrast,evidence=self.prepare(root,samples())
                if kind=='missing':evidence.pop('S0')
                if kind=='extra':evidence['bogus']=evidence['S0']
                if kind=='zero':Path(evidence['S0']['log']).write_text('Total 0 usable reads were sampled\n')
                if kind=='empty_log':Path(evidence['S0']['log']).write_text('')
                if kind=='ambiguous':Path(evidence['S0']['report']).write_text(native_text(.6,.4))
                with self.assertRaises(InputError):resolve(sheet,contrast,evidence,root/'out')
                self.assertFalse((root/'out').exists())


if __name__=='__main__':unittest.main()
