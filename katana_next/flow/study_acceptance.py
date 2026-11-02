"""Positive acceptance of the installed whole-study CoSplice engine."""
from pathlib import Path
import json
import re
import pandas as pd

from .model import read_json, file_sha256


def whole_study_report_gaps(root, plan, report, report_dir):
    root, report_dir = Path(root).resolve(), Path(report_dir).resolve()
    try:
        tasks = [t for t in plan['tasks'] if t['workflow_id'] == 'programs.cosplice']
        if len(tasks) != 1:
            raise ValueError('Exactly one whole-study network task must be planned')
        task = tasks[0]
        collection = read_json(root/'collection.receipt.json')
        records = [r for r in collection['task_receipts'] if r['task_id'] == task['id']]
        if len(records) != 1:
            raise ValueError('The whole-study task completion record is absent')
        receipt = read_json(root/records[0]['collection_path'])
        if receipt.get('status') not in {'completed','reused'} or receipt.get('job',{}).get('status') not in {'complete','cache_hit'}:
            raise ValueError('Whole-study network computation has not completed')
        transport = {r['source']['path']: r for r in read_json(root/'transport-provenance.json')['products']}
        def member(source):
            record = transport[str(source)]
            path = (root/record['collection_path']).resolve()
            if not path.is_relative_to(root) or file_sha256(path) != record['source']['sha256']:
                raise ValueError('A native whole-study output is absent or changed')
            return path
        artifacts = [a for a in receipt['generated_artifacts'] if a.get('metadata',{}).get('report_role') == 'whole_study_chapter']
        if len(artifacts) != 1:
            raise ValueError('The completed whole-study module browser is missing')
        browser = artifacts[0]
        from .presentation import _producer_path
        base = _producer_path(browser['path']).parent
        # All these products must come from this exact completed task, not a preview.
        files = ['status.json','summary.json','study/manifest.json','samples.tsv','model/status.json',
                 'model/native_model.rds','model/assignments.tsv','study/scores.tsv','study/module_basis.tsv',
                 'interpretations/default/associations.tsv','results.zip','cosplice_chapter6.html']
        paths = {name: member(base/name) for name in files}
        if read_json(paths['status.json']).get('status') != 'completed' or read_json(paths['model/status.json']).get('status') != 'completed':
            raise ValueError('Native WGCNA did not complete')
        meta = pd.read_csv(paths['samples.tsv'],sep='\t',dtype=str)
        project = read_json(root/'project.json')
        expected = {s['sample_id']:s['biological_sample_id'] for s in project['samples']}
        if not meta.sample_id.is_unique or meta.set_index('sample_id').biological_sample_id.to_dict() != expected:
            raise ValueError('Whole-study network sample or biological-unit identities differ')
        study = read_json(paths['study/manifest.json'])
        if study['samples'] != len(expected) or study['biological_units'] != len(set(expected.values())):
            raise ValueError('Repeated biological units were miscounted')
        for name, digest in study['frozen_files'].items():
            if file_sha256(member(base/name)) != digest:
                raise ValueError('Frozen whole-study scientific output changed: '+name)
        copies = [r for r in report.get('complete_results',[]) if r.get('metadata',{}).get('report_role') == 'whole_study_chapter']
        if len(copies) != 1:
            raise ValueError('The final report omits the whole-study network browser')
        copied = (report_dir/copies[0]['source']).resolve()
        if not copied.is_relative_to(report_dir) or not copied.is_file() or not copied.stat().st_size:
            raise ValueError('The report module browser is missing or empty')
        payload = re.search(r'<script type="application/json" id="reportData">(.*?)</script>',
                            copied.read_text(encoding='utf8'), re.S)
        if payload is None:
            raise ValueError('The network chapter is missing its interactive figure data')
        data = json.loads(payload[1])
        if data['manifest']['network_id'] != study['network_id']:
            raise ValueError('The network chapter belongs to a different native fit')
        if {s['sample_id'] for s in data['samples']} != set(expected):
            raise ValueError('The network chapter omits study samples')
        for module in data['modules']:
            key = module['module']
            graph = data['graphs'].get(key)
            if graph is None or graph['total_nodes'] != module['events']:
                raise ValueError('The network chapter omits a module graph: '+key)
            if key not in data['module_files'] or not any(g['module'] == key for g in data['genes']):
                raise ValueError('The network chapter omits a module gene list: '+key)
        browsers = [b for b in report.get('network_browsers', []) if b.get('whole_study')]
        if len(browsers) != 1 or browsers[0]['source'] != copies[0]['source']:
            raise ValueError('The completed module browser is not displayed in the final report')
        # Copying HTML may rewrite relative links; native numerical files must remain exact.
        if not paths['model/native_model.rds'].stat().st_size or not paths['results.zip'].stat().st_size:
            raise ValueError('Whole-study native model or download is empty')
        return []
    except (ValueError,OSError,KeyError,TypeError) as exc:
        return [{'component':'network','reason':str(exc)}]
