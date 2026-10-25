"""Fresh and refreshed reports share the approved compact RBP presentation."""
from copy import deepcopy
from types import SimpleNamespace

from katana_next.reports.presentation import prepare_presentation


def _builder(tmp_path):
    return SimpleNamespace(
        source={'project': {'contrasts': [{'contrast_id': 'B_vs_A'}]}}, output=tmp_path,
        tables=[], network_browsers=[], panels=[], exclusions=[], derived_audits=[],
        enrichment_summaries=[{'id': 'rbp-enrichment-saved', 'contrast_id': 'B_vs_A',
                               'cell_type': 'all', 'rows': [], 'tested_sets': 3}],
    )


def test_saved_rbp_charts_stay_in_outputs_with_same_scope_summary(tmp_path):
    builder = _builder(tmp_path)
    builder.panels = [{'id': layer, 'chapter': 7, 'title': f'RBP identity | {layer} | B_vs_A',
                       'contrast_id': 'B_vs_A', 'cell_type': 'all', 'policy': None,
                       'image': f'outputs/{layer}.png',
                       'data_links': [{'role': 'full_results', 'path': f'outputs/{layer}.parquet'}]}
                      for layer in ('DEG', 'DAS', 'RS-DAS', 'DTU', 'DEXSeq')]
    before = deepcopy(builder.panels)
    prepare_presentation(builder)
    assert len(builder.panels) == 5
    assert all(row['output_only'] for row in builder.panels)
    assert [(row['image'], row['data_links'], row['policy']) for row in builder.panels] == [
        (row['image'], row['data_links'], row['policy']) for row in before]
    assert builder.enrichment_summaries[0]['rows'] == []  # Zero significant sets is still a valid summary.
    refreshed = deepcopy(builder)
    prepare_presentation(refreshed)
    assert refreshed.panels == builder.panels


def test_missing_or_other_scope_summary_never_hides_rbp_chart(tmp_path):
    builder = _builder(tmp_path)
    builder.panels = [
        {'chapter': 7, 'title': 'RBP identity | DEG | D_vs_C', 'contrast_id': 'D_vs_C'},
        {'chapter': 7, 'title': 'RBP identity | DEG | B_vs_A', 'contrast_id': 'B_vs_A', 'cell_type': 'T'},
        {'chapter': 7, 'title': 'Boundary motif enrichment | B_vs_A', 'contrast_id': 'B_vs_A'},
    ]
    prepare_presentation(builder)
    assert not any(row.get('output_only') for row in builder.panels)
    builder.enrichment_summaries = []
    builder.panels = [{'chapter': 7, 'title': 'RBP identity | DEG | B_vs_A', 'contrast_id': 'B_vs_A'}]
    prepare_presentation(builder)
    assert not builder.panels[0].get('output_only')


def test_unselected_summaries_are_not_inserted_into_report_body(tmp_path):
    builder = _builder(tmp_path)
    builder.enrichment_summaries += [{'id': 'rbp-enrichment-other', 'contrast_id': 'D_vs_C',
                                     'full_path': 'outputs/native-all-comparisons.parquet'}]
    prepare_presentation(builder)
    assert [row['contrast_id'] for row in builder.enrichment_summaries] == ['B_vs_A']
