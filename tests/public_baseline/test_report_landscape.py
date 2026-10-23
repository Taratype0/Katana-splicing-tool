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
    heatmap = next(p for p in result['panels'] if p['id'] == 'psi-heatmap')
    plot_path = next(link['path'] for link in heatmap['data_links'] if link['role'] == 'plot_data')
    plotted = pd.read_csv(plot_path, sep='\t')
    original = pd.read_parquet(matrix).set_index('event_id')
    pd.testing.assert_frame_equal(plotted.set_index('event_id')[list(original.columns)], original.loc[plotted.event_id], check_dtype=False)
    assert [link['role'] for link in heatmap['data_links']] == ['plot_data', 'full_results', 'full_results', 'methods']
    assert all(Path(link['path']).is_file() for panel in result['panels'] for link in panel['data_links'])


def test_pca_and_consistency_data_are_saved_coordinates_and_actual_plot_rows(tmp_path):
    samples = [{'sample_id': s, 'condition': 'A' if s.startswith('a') else 'B'} for s in ['a1', 'a2', 'b1', 'b2']]
    matrix = tmp_path / 'matrix.parquet'
    pd.DataFrame({'event_id': ['e1', 'e2'], 'a1': [.1, .8], 'a2': [.2, .7], 'b1': [.8, .2], 'b2': [.9, .1]}).to_parquet(matrix)
    pca = tmp_path / 'pca.parquet'
    coordinates = pd.DataFrame({'sample_id': ['a1', 'a2', 'b1', 'b2'], 'condition': ['A', 'A', 'B', 'B'],
                                'PC1': [-2.125, -1.875, 1.875, 2.125], 'PC2': [.1, -.1, .1, -.1],
                                'explained_PC1': [.82] * 4, 'explained_PC2': [.13] * 4})
    coordinates.to_parquet(pca)
    similarity = tmp_path / 'similarity.parquet'
    pairs = pd.DataFrame([{'sample_a': a, 'sample_b': b, 'spearman_r': 1. if a == b else (.9 if a[0] == b[0] else .4)}
                          for a in coordinates.sample_id for b in coordinates.sample_id])
    pairs.to_parquet(similarity)
    artifacts = [{'name': name, 'path': str(path)} for name, path in [('Sample-level PSI matrix', matrix), ('Sample PSI PCA', pca), ('Sample PSI similarity', similarity)]]
    original = {p: (digest(p), p.stat().st_mtime_ns) for p in (matrix, pca, similarity)}
    result = build_landscape_views({'samples': samples}, artifacts, tmp_path / 'new')
    def plotted(identity):
        panel = next(p for p in result['panels'] if p['id'] == identity)
        return pd.read_csv(next(link['path'] for link in panel['data_links'] if link['role'] == 'plot_data'), sep='\t')
    pd.testing.assert_frame_equal(plotted('psi-pca'), coordinates)
    rows = plotted('replicate-consistency')
    assert len(rows) == 12 and not (rows.sample_a == rows.sample_b).any()
    assert set(rows.pair_type) == {'Within condition', 'Between conditions'}
    assert set(plotted('sample-consistency-diagnostic').median_within_condition_similarity) == {.9}
    assert all((digest(p), p.stat().st_mtime_ns) == identity for p, identity in original.items())


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
