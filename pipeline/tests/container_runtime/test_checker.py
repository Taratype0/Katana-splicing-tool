"""Pure receipt/trace validation tests; these never launch Nextflow or Docker."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


spec = importlib.util.spec_from_file_location('container_checker', Path(__file__).resolve().parents[2] / 'tools/check_container_runtime.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class EvidenceChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)
        self.image = 'sha256:' + 'a' * 64

    def trace(self, status='COMPLETED', kinds=('network', 'mapp', 'motif'), exit_code='0'):
        path = self.root / 'trace.tsv'
        with path.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['name', 'status', 'exit', 'container', 'workdir'], delimiter='\t')
            writer.writeheader()
            for kind in kinds:
                writer.writerow({'name': 'IMAGE_RUNTIME (' + kind + ')', 'status': status,
                    'exit': exit_code, 'container': self.image, 'workdir': str(self.root)})
        return path

    def test_actual_and_cached_trace_contract(self):
        for status in ('COMPLETED', 'CACHED'):
            self.assertEqual(len(checker.read_trace(self.trace(status), status, self.image)), 3)

    def test_failed_or_duplicate_task_rejected(self):
        for path in (self.trace(exit_code='1'),):
            with self.assertRaises(AssertionError):
                checker.read_trace(path, 'COMPLETED', self.image)
        with self.assertRaises(AssertionError):
            checker.read_trace(self.trace(kinds=('network', 'network', 'mapp')), 'COMPLETED', self.image)

    def test_wrong_container_identity_rejected(self):
        with self.assertRaises(AssertionError):
            checker.read_trace(self.trace(), 'COMPLETED', 'sha256:' + 'b' * 64)

    def products(self):
        for kind in checker.KINDS:
            out = self.root / kind; out.mkdir()
            value = {'status': 'ready', 'probe_kind': kind, 'installed_prefix': '/opt/conda',
                'task_uid': 1000, 'task_gid': 1000, 'host_runtime_rejected': True,
                'existing_mounted_runtime_tested': True, 'host_runtime_executed': False,
                'native_cli_runtime_check': True, 'scientific_computation': False,
                'host_runtime_error': ['external host path']}
            (out / 'acceptance.json').write_text(json.dumps(value), encoding='utf-8')
            (out / 'katana-flow-help.stdout.txt').write_text('usage: katana-flow', encoding='utf-8')

    def test_host_rejection_proof_required(self):
        self.products()
        self.assertEqual(set(checker.validate_products(self.root, 1000, 1000)), checker.KINDS)
        path = self.root / 'motif/acceptance.json'
        value = json.loads(path.read_text()); value['existing_mounted_runtime_tested'] = False
        path.write_text(json.dumps(value))
        with self.assertRaises(AssertionError):
            checker.validate_products(self.root, 1000, 1000)

    def test_output_hash_detects_changed_data(self):
        file = self.root / 'runtime.json'; file.write_text('first')
        before = checker.inventory(self.root)
        file.write_text('other')
        self.assertNotEqual(before, checker.inventory(self.root))


if __name__ == '__main__':
    unittest.main()
