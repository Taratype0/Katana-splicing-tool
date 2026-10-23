"""An old caller result can feed new work without relabelling or refitting it."""
import pandas as pd
import pytest

from katana_next.flow.eligibility_projection import project_existing_as
from katana_next.flow.model import file_sha256, write_json


def fixture(tmp_path):
    numerator = ['b1', 'b2', 'b3']
    denominator = ['a1', 'a2', 'a3']
    project = {'project_id': 'same-study', 'name': 'Projection fixture',
        'input_root': str(tmp_path / 'input'), 'output_root': str(tmp_path / 'old-results'),
        'annotation': {'genome_build': 'test'},
        'samples': [{'sample_id': s, 'biological_sample_id': s, 'condition': s[0]}
                    for s in numerator + denominator],
        'contrasts': [{'contrast_id': 'B_vs_A', 'label': 'B versus A',
            'numerator': {'label': 'B', 'sample_ids': numerator},
            'denominator': {'label': 'A', 'sample_ids': denominator}}],
        'assets': [{'asset_id': 'rmats-SE', 'kind': 'rmats_result', 'path': str(tmp_path/'original.txt'),
            'metadata': {'contrast_id': 'B_vs_A', 'event_type': 'SE', 'reverse_direction': False,
                         'numerator_sample_ids': list(numerator), 'denominator_sample_ids': list(denominator)}}]}
    rows = []
    for index, counts in enumerate(['10,10,9', '10,10,10', '100,100,100']):
        rows.append({'row_id': 'B_vs_A::e'+str(index), 'event_id': 'e'+str(index),
            'comparison_id': 'B_vs_A', 'event_type': 'SE', 'gene_id': 'g'+str(index),
            'chr': 'chr1', 'strand': '+', 'numerator': 'B', 'denominator': 'A',
            'IncLevel1': '0.2,0.9,0.8', 'IncLevel2': '0.4,0.3,0.1',
            'IJC_SAMPLE_1': counts, 'SJC_SAMPLE_1': '0,0,0',
            'IJC_SAMPLE_2': '20,20,20', 'SJC_SAMPLE_2': '0,0,0',
            'significant': index != 2, 'stable_significant': index != 2,
            'measurement_eligible': True, 'passes_direction_consistency': False,
            'criterion_min_supported_replicates': 2, 'reverse_direction': False,
            'PValue': 0.001 + index/10, 'FDR': 0.02 + index/10,
            'dpsi': -0.15, 'IncLevelDifference': -0.15})
    table = tmp_path / 'old.parquet'
    pd.DataFrame(rows).to_parquet(table, index=False)
    manifest = write_json(tmp_path/'project.json', project)
    return table, manifest, project


def test_all_sample_projection_preserves_statistics_and_legacy_flags(tmp_path):
    table, manifest, _ = fixture(tmp_path)
    before = (table.read_bytes(), manifest.read_bytes())
    receipt = project_existing_as(table, manifest, tmp_path/'projected', contrast_ids=['B_vs_A'])
    result = pd.read_parquet(receipt['output']['path'])
    assert result.measurement_eligible.tolist() == [False, True, True]
    assert result.stable_significant.tolist() == [False, True, False]
    assert result.source_stable_significant.tolist() == [True, True, False]
    assert result.source_passes_direction_consistency.eq(False).all()
    assert 'passes_direction_consistency' not in result
    assert receipt['native_statistics_reestimated'] is False
    assert receipt['count_vector_order_verified_against_declared_samples'] is True
    native = pd.read_parquet(table)
    for column in ('PValue', 'FDR', 'dpsi', 'IncLevelDifference', 'significant', 'IJC_SAMPLE_1'):
        pd.testing.assert_series_equal(native[column], result[column], check_exact=True)
    assert before == (table.read_bytes(), manifest.read_bytes())
    with pytest.raises(FileExistsError):
        project_existing_as(table, manifest, tmp_path/'projected', contrast_ids=['B_vs_A'])


@pytest.mark.parametrize('change', ['order', 'direction', 'unknown', 'duplicate'])
def test_projection_rejects_ambiguous_or_different_identity(tmp_path, change):
    table, manifest, project = fixture(tmp_path)
    contrasts = ['B_vs_A']
    if change == 'order':
        project['assets'][0]['metadata']['numerator_sample_ids'].reverse()
    elif change == 'direction':
        project['assets'][0]['metadata']['reverse_direction'] = True
    elif change == 'unknown':
        contrasts = ['A_vs_B']
    else:
        frame = pd.read_parquet(table)
        pd.concat([frame, frame.iloc[:1]]).to_parquet(table, index=False)
    write_json(manifest, project)
    with pytest.raises(ValueError):
        project_existing_as(table, manifest, tmp_path/'rejected', contrast_ids=contrasts)
    assert not (tmp_path/'rejected/projection.receipt.json').exists()


def test_stage_publishes_input_preparation_identity_not_new_as_workflow(tmp_path):
    from katana_next.flow.stages import execute_stage
    table, manifest, _ = fixture(tmp_path)
    result = execute_stage({'adapter': 'as.input_projection',
        'inputs': {'source': str(table), 'project': str(manifest)},
        'scope': {'contrast_ids': ['B_vs_A']}}, tmp_path/'task')
    output = next(p for p in result['outputs'] if p.get('kind') == 'table')
    assert output['analysis_id'] == 'flow.as_input_projection'
    assert output['metadata']['scientific_model_fitted'] is False
    assert output['metadata']['native_statistics_reestimated'] is False
    assert output['sha256'] == file_sha256(output['path'])
