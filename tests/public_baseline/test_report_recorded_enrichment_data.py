import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_next.reports.rbp_enrichment import build_rbp_summaries
from katana_next.reports.saved_plot_data import recorded_plot_data_links
from katana_next.reports.bootstrap_plot_data import bootstrap_plot_links
from katana_next.reports.layout import _rbp_summary


def test_boundary_motif_contributors_use_only_native_matched_foreground(tmp_path):
    from katana_next.reports.rbp_enrichment import build_motif_hit_summaries
    frames = {
        'motif_enrichment.parquet':pd.DataFrame([dict(comparison_id='A', RBP='Rbp1', evidence_class='direct', enrichment_pvalue=.001, enrichment_fdr=.01, foreground_hit_genes=1, hit_fraction_difference=.5)]),
        'motif_catalog.parquet':pd.DataFrame([dict(RBP='Rbp1', evidence_class='direct', motif_token='m')]),
        'motif_matching.parquet':pd.DataFrame([dict(comparison_id='A',status='matched',foreground_event='e1',foreground_gene='GeneA'),dict(comparison_id='A',status='unmatched',foreground_event='e2',foreground_gene='Excluded')]),
        'fimo_sites.parquet':pd.DataFrame([dict(motif_token='m',event_id='e1'),dict(motif_token='m',event_id='e2'),dict(motif_token='m',event_id='background')]),
    }
    artifacts=[]
    for name, data in frames.items():
        path=tmp_path/name; data.to_parquet(path)
        artifacts.append(dict(path=str(path),analysis_id='prioritize.event_motifs',run_fingerprint='native'))
    records=build_motif_hit_summaries(artifacts)
    assert records[0]['rows'][0]['hit_genes'] == ['GeneA']
    assert records[0]['rows'][0]['hit_events'] == ['e1']
    assert records[0]['rows'][0]['fdr'] == .01
    assert 'Boundary motif enrichment' in _rbp_summary(records[0])
    with pytest.raises(ValueError,match='registered'):
        build_motif_hit_summaries(artifacts[:-1])


def test_rbp_hit_genes_and_significance_stay_comparison_specific(tmp_path):
    path = tmp_path / 'rbp_enrichment.parquet'
    pd.DataFrame([
        dict(comparison_id='A', term='RBP membership', genes='GeneA, GeneB', overlap_genes=2, pvalue=.001, fdr=.01, fold_enrichment=4),
        dict(comparison_id='A', term='nonsignificant', genes='Other', overlap_genes=1, pvalue=.2, fdr=.5, fold_enrichment=1),
        dict(comparison_id='B', term='RBP membership', genes='GeneC', overlap_genes=1, pvalue=.004, fdr=.05, fold_enrichment=3),
    ]).to_parquet(path)
    artifact = dict(path=str(path), analysis_id='prioritize.rbp_support', run_fingerprint='native', artifact_id='r')
    before = path.read_bytes()
    a, b = build_rbp_summaries([artifact, artifact])
    assert a['tested_sets'] == 2 and len(a['rows']) == 1
    assert a['rows'][0]['hit_genes'] == ['GeneA', 'GeneB'] and b['rows'][0]['hit_genes'] == ['GeneC']
    html = _rbp_summary(a)
    assert 'GeneA, GeneB' in html and 'GeneC' not in html and 'nonsignificant' not in html
    assert '0.001' in html and '0.01' in html and '<details>' in html
    assert path.read_bytes() == before
    invalid = pd.read_parquet(path); invalid.loc[0, 'overlap_genes'] = 3; invalid.to_parquet(path)
    with pytest.raises(ValueError, match='hit count'):
        build_rbp_summaries([artifact])


def test_saved_vega_exact_rows_and_exports_are_bound_without_guessing(tmp_path):
    spec = tmp_path / 'saved.vl.json'
    rows = [{'term': 'set', 'pvalue': .001, 'fdr': .02, 'comparison_id': 'A'}]
    spec.write_text(json.dumps({'$schema': 'https://vega.github.io/schema/vega-lite/v5.json', 'mark': 'point', 'data': {'values': rows}, 'encoding': {'x': {'field': 'fdr'}}}))
    panel = {'contrast_id': 'A', 'policy': None}
    links = recorded_plot_data_links({'path': str(spec)}, panel, tmp_path / 'report')
    assert json.loads(Path(links[0]['path']).read_text()) == rows
    image = tmp_path / 'export.png'; image.write_bytes(b'owned export')
    binding = dict(path=str(spec), source_sha256=hashlib.sha256(spec.read_bytes()).hexdigest(), export_sha256=hashlib.sha256(image.read_bytes()).hexdigest())
    assert recorded_plot_data_links({'path': str(image), 'metadata': {'saved_plot_specification': binding}}, panel, tmp_path / 'export-report')
    image.write_bytes(b'changed')
    with pytest.raises(ValueError, match='export SHA'):
        recorded_plot_data_links({'path': str(image), 'metadata': {'saved_plot_specification': binding}}, panel, tmp_path / 'bad')
    with pytest.raises(ValueError, match='different comparison'):
        recorded_plot_data_links({'path': str(spec)}, {'contrast_id': 'B'}, tmp_path / 'bad-scope')


def test_saved_transform_is_not_falsely_promoted_to_exact_values(tmp_path):
    spec = tmp_path / 'transformed.vl.json'
    spec.write_text(json.dumps({'$schema': 'vega-lite', 'mark': 'bar', 'data': {'values': [{'x': 1}]}, 'transform': [{'filter': 'datum.x > 0'}]}))
    assert recorded_plot_data_links({'path': str(spec)}, {}, tmp_path / 'report') == []


def test_native_missing_statistic_is_null_never_zero(tmp_path):
    spec = tmp_path / 'missing.vl.json'
    spec.write_text('{"$schema":"vega-lite", "mark":"point", "data":{"values":[{"pvalue":NaN}]}}')
    before = spec.read_bytes()
    links = recorded_plot_data_links({'path':str(spec)}, {}, tmp_path/'report')
    assert json.loads(Path(links[0]['path']).read_text()) == [{'pvalue':None}]
    assert spec.read_bytes() == before


@pytest.mark.parametrize('key_dtype', [str, object])
def test_bootstrap_plot_bindings_match_native_subset_order_and_sizes(tmp_path, key_dtype):
    root = tmp_path / '04c_consensus_modules'; (root / 'plots').mkdir(parents=True)
    pd.DataFrame({'event_key': ['e0','e1','e2','e3'], 'consensus_module': ['M2','M1','grey','M1'], 'consensus_membership_confidence': [.5,.3,.9,.8]}).to_csv(root / 'consensus_event_module_assignments.tsv', sep='\t', index=False)
    np.savez(root / 'bootstrap_coassignment_matrix.npz', event_keys=np.array(['e0','e1','e2','e3'], dtype=key_dtype), coassignment=np.arange(16).reshape(4,4)/16)
    pd.DataFrame({'consensus_module':['M2','M1'], 'n_events':[1,2], 'mean_event_confidence':[.5,.55]}).to_csv(root / 'consensus_module_summary.tsv', sep='\t', index=False)
    for key, filename in [('heatmap','bootstrap_coassignment_heatmap.png'), ('sizes','consensus_module_sizes.png')]:
        image = root / 'plots' / filename; image.write_bytes(b'native plot')
        links = bootstrap_plot_links({'path':str(image), 'metadata':{'flow_stage':'network.bootstrap','flow_artifact_key':key}}, {'contrast_id':'A','policy':None}, tmp_path/'report', [])
        if key == 'heatmap':
            data = json.loads(Path(links[0]['path']).read_text())
            assert data['event_keys'] == ['e3','e1','e0']
            assert data['coassignment'] == (np.arange(16).reshape(4,4)[np.ix_([3,1,0],[3,1,0])]/16).tolist()
        else:
            assert pd.read_csv(links[0]['path'],sep='\t').n_events.tolist() == [2,1]
