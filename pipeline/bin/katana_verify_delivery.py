#!/usr/bin/env python3
"""Final Nextflow task: verify retained client data and the unchanged report."""
import argparse
import hashlib
import json
from pathlib import Path


def verify(project, report, out):
    from katana_next.flow.storage_project import validate_retained_project
    project, report = Path(project).resolve(), Path(report).resolve()
    validate_retained_project(project)
    retention = json.loads((project / 'retention.receipt.json').read_text())
    acceptance = json.loads((report / 'flow.acceptance.json').read_text())
    plan = json.loads((project / 'collection/plan.json').read_text())
    if (acceptance.get('schema') != 'katana.flow.report_acceptance/v1'
            or acceptance.get('status') != 'completed'
            or acceptance.get('plan_id') != plan['plan_id']):
        raise ValueError('Report and retained project must have the same completed plan')
    files = []
    for name in ('REPORT.html', 'REPORT.zip', 'report-manifest.json', 'flow.acceptance.json'):
        path = report / name
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('Missing report artifact: ' + name)
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        files.append(dict(path='report/' + name, bytes=path.stat().st_size, sha256=digest))
        if name == 'REPORT.zip' and digest != acceptance['archive_sha256']:
            raise ValueError('Report archive changed after acceptance')
    result = dict(schema='katana.nextflow.delivery/v1', status='completed',
                  plan_id=plan['plan_id'], project='katana/project/project.json',
                  project_receipt_sha256=retention['receipt_sha256'], report_files=files,
                  orchestration='nextflow', scientific_recalculation=False)
    with Path(out).open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    verify(args.project, args.report, args.out)
