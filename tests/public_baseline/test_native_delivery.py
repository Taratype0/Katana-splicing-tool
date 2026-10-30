"""Final native gate must bind the exact retained plan and accepted archive."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('native_delivery', ROOT / 'pipeline/bin/katana_verify_delivery.py')
delivery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(delivery)


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    import katana_next.flow.storage_project as storage
    observed = []
    # Actual retention, relocation and original-work deletion are exercised by
    # test_flow_storage_project and installed native delivery integration.
    monkeypatch.setattr(storage, 'validate_retained_project', lambda path: observed.append(path))
    project, report = tmp_path / 'project', tmp_path / 'report'
    (project / 'collection').mkdir(parents=True)
    report.mkdir()
    (project / 'collection/plan.json').write_text(json.dumps({'plan_id': 'generated-plan'}))
    (project / 'retention.receipt.json').write_text(json.dumps({'receipt_sha256': 'generated-receipt'}))
    for name in ('REPORT.zip', 'REPORT.html', 'report-manifest.json'):
        (report / name).write_text('generated gate fixture')
    accepted = {'schema': 'katana.flow.report_acceptance/v1', 'status': 'completed',
                'plan_id': 'generated-plan', 'archive_sha256': hashlib.sha256((report / 'REPORT.zip').read_bytes()).hexdigest()}
    (report / 'flow.acceptance.json').write_text(json.dumps(accepted))
    return project, report, tmp_path / 'receipt.json', observed


def test_delivery_records_relative_paths_and_validates_client_project(inputs):
    project, report, out, observed = inputs
    result = delivery.verify(project, report, out)
    assert result['status'] == 'completed'
    assert observed == [project.resolve()]
    assert all(not Path(row['path']).is_absolute() for row in result['report_files'])


@pytest.mark.parametrize('damage', ['plan', 'archive', 'missing_html', 'failed_report'])
def test_incomplete_or_mixed_delivery_has_no_success_receipt(inputs, damage):
    project, report, out, _ = inputs
    if damage == 'plan':
        (project / 'collection/plan.json').write_text(json.dumps({'plan_id': 'another-plan'}))
    elif damage == 'archive':
        (report / 'REPORT.zip').write_text('changed bytes')
    elif damage == 'missing_html':
        (report / 'REPORT.html').unlink()
    else:
        acceptance = report / 'flow.acceptance.json'
        data = json.loads(acceptance.read_text()); data['status'] = 'failed'
        acceptance.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        delivery.verify(project, report, out)
    assert not out.exists()
