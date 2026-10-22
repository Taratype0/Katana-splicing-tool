"""Saved measurements, portable all-row views and explicit provenance boundaries."""
from pathlib import Path
import json
import pandas as pd
import pytest
from katana_next.reports.landscape import build_landscape_views
from katana_next.reports.sources import load_source, digest


def test_landscape_preserves_missing_measurements_and_all_matrix_rows(tmp_path):
    samples = [{'sample_id': x, 'biological_sample_id': x, 'condition': 'A' if x.startswith('a') else 'B'} for x in ['a1', 'a2', 'b1', 'b2']]
    matrix = tmp_path/'saved.parquet'
    pd.DataFrame({'event_id': ['one', 'two', 'three'], 'a1': [.1, None, .8], 'a2': [.2, .3, .7], 'b1': [.5, .4, .9], 'b2': [.6, .5, .8]}).to_parquet(matrix)
    before = (digest(matrix), matrix.stat().st_mtime_ns)
    result = build_landscape_views({'samples': samples}, [{'name': 'Sample-level PSI matrix', 'path': str(matrix), 'kind': 'table'}], tmp_path/'new')
    browser = next(Path(t['path']) for t in result['tables'] if t.get('kind') == 'html')
    raw = browser.read_text('utf-8').split('<script id="matrix-data" type="application/json">')[1].split('</script>')[0]
    data = json.loads(raw)
    assert data['events'] == ['one', 'two', 'three']
    assert data['values'][1][0] is None and data['values'][0][0] == .1
    assert before == (digest(matrix), matrix.stat().st_mtime_ns)
    assert not any('pca' == p['id'] for p in result['panels'])


def test_missing_matrix_does_not_generate_fake_results(tmp_path):
    result = build_landscape_views({'samples': []}, [], tmp_path/'new')
    assert result['panels'] == []
    assert result['unavailable'][0]['chapter'] == 2


def test_supplemental_results_require_exact_parent_cohort_and_hash(tmp_path):
    original = {'project_id': 'parent', 'output_root': str(tmp_path/'old'), 'samples': [{'sample_id': 's1', 'biological_sample_id': 's1', 'condition': 'A'}], 'contrasts': [], 'metadata': {}}
    manifest = tmp_path/'katana.project.json'; manifest.write_text(json.dumps(original), encoding='utf-8')
    child = {**original, 'project_id': 'child', 'metadata': {'source_project_manifest_SHA256': digest(manifest), 'source_project_id': 'parent'}}
    descriptor = {'project_manifest': str(manifest), 'artifacts': [], 'supplemental_sources': [{'project': child, 'artifacts': []}]}
    assert load_source(descriptor)['supplemental'][0]['project']['project_id'] == 'child'
    child['samples'] = [{'sample_id': 'another', 'condition': 'A'}]
    with pytest.raises(ValueError, match='exact parent'):
        load_source(descriptor)
