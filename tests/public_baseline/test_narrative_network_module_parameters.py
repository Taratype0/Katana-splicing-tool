"""Report only the settings of the network that actually supplied the cards."""
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import zipfile

import pandas as pd
import pytest

from katana_next.full_fixture import create_full_fixture
from katana_next.services.project_service import ProjectService
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows import narrative_report


@pytest.fixture
def context(tmp_path):
    return ProjectService(application_root=tmp_path / 'application').open_project(
        create_full_fixture(tmp_path / 'project'))


def network(context, run, size, height):
    writer = ArtifactWriter(context, 'programs.cosplice', run * 64, {}, 'test',
        common_metadata={'network_module_parameters': {
            'min_module_size': size, 'merge_cut_height': height}})
    scores = writer.table('Scores', 'scores.parquet',
        pd.DataFrame({'sample_id': ['s1'], 'C1': [.1]}), entity_type='module')
    confidence = writer.table('Confidence', 'confidence.parquet',
        pd.DataFrame({'consensus_module': ['C1'], 'n_events': [size]}), entity_type='module')
    writer.finish()
    return scores, confidence


def delivery(context, scores, confidence):
    writer = ArtifactWriter(context, 'programs.explain_expression', 'c' * 64, {}, 'test',
        common_metadata={'input_files': {
            'scores': scores.metadata['exports']['tsv'],
            'confidence': confidence.metadata['exports']['tsv']}})
    cards = writer.table('Module interpretation cards', 'module_cards.parquet',
        pd.DataFrame({'module': ['C1'], 'n_events': [37]}), entity_type='module')
    writer.finish()
    root = context.manifest.output_root / 'delivery'
    target = root / '23_Module_explanation' / 'selected' / 'module_cards.tsv'
    target.parent.mkdir(parents=True)
    shutil.copy2(cards.metadata['exports']['tsv'], target)
    return {'study': 'Toy', 'result_root': str(root),
            'manifest_path': str(context.manifest_path)}, target, cards


def test_actual_bound_parameters_do_not_follow_a_newer_network(context):
    config, cards, card_artifact = delivery(context, *network(context, 'a', 37, .12))
    unrelated = network(context, 'b', 20, .25)
    record = narrative_report._network_module_parameter_record(context, config)
    assert record['status'] == 'recorded'
    assert record['network_module_parameters'] == {'min_module_size': 37, 'merge_cut_height': .12}
    assert record['network_run_fingerprint'] == 'a' * 64
    assert record['module_cards_artifact_id'] == card_artifact.artifact_id
    assert set(record['network_source_artifact_ids']).isdisjoint(a.artifact_id for a in unrelated)


@pytest.mark.parametrize('change', ['unrecorded_cards', 'modified_export', 'modified_network', 'mixed_runs'])
def test_missing_or_broken_source_chain_does_not_substitute_defaults(context, change):
    sources = network(context, 'a', 37, .12)
    if change == 'mixed_runs':
        sources = sources[0], network(context, 'b', 37, .12)[1]
    config, cards, _ = delivery(context, *sources)
    if change == 'unrecorded_cards':
        cards.write_text('module\tn_events\nC2\t90\n', encoding='utf-8')
    elif change == 'modified_export':
        Path(sources[0].metadata['exports']['tsv']).write_text('sample_id\tC1\ns1\t9\n', encoding='utf-8')
    elif change == 'modified_network':
        pd.DataFrame({'sample_id': ['s1'], 'C1': [9.0]}).to_parquet(sources[0].path, index=False)
    record = narrative_report._network_module_parameter_record(context, config)
    assert record['status'] == 'not_recorded'
    assert 'network_module_parameters' not in record


def test_new_report_renders_bound_custom_values_and_packages_provenance(context, tmp_path, monkeypatch):
    config, _, _ = delivery(context, *network(context, 'a', 37, .12))
    root = Path(config['result_root'])
    summary = root / '20_Presentation_report' / 'comparison_totals.tsv'
    summary.parent.mkdir()
    summary.write_text('comparison_id\nA_vs_B\n', encoding='utf-8')
    config_path = tmp_path / 'config.json'
    config_path.write_text(json.dumps(config), encoding='utf-8')
    historical = root / '24_Narrative_report' / 'REPORT_EN.html'
    historical.parent.mkdir()
    historical.write_text('<html>Previously approved report</html>', encoding='utf-8')
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file()}

    def build(study, base, manifest, out):
        out.mkdir(parents=True)
        (out / 'REPORT_EN.html').write_text('<html><body><main><section>Original study chapters</section></main></body></html>', encoding='utf-8')
        (out / 'SOURCE_MANIFEST_SHA256.tsv').write_text('path\tsha256\n', encoding='utf-8')
        (out / 'report_status.json').write_text('{}', encoding='utf-8')

    monkeypatch.setattr(narrative_report, 'build_report', build)
    workflow = narrative_report.NarrativeReportWorkflow()
    parents = workflow.parent_artifacts(context, {'config_path': str(config_path)})
    assert len(parents) == 3
    workflow.run(context, SimpleNamespace(), 'd' * 64, {'config_path': str(config_path)})
    out = root / '24_Narrative_report' / ('d' * 16)
    text = (out / 'REPORT_EN.html').read_text(encoding='utf-8')
    assert 'Original study chapters' in text
    assert 'Minimum module size: 37 splice events' in text
    assert 'Module merge height: 0.12' in text
    assert 'NETWORK_MODULE_PARAMETERS.json' in text
    assert before == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before}
    archive = next(context.results.artifact_directory(workflow.definition.workflow_id, 'd' * 64).glob('*.zip'))
    with zipfile.ZipFile(archive) as packaged:
        attachment = next(n for n in packaged.namelist() if n.endswith('/NETWORK_MODULE_PARAMETERS.json'))
        assert json.loads(packaged.read(attachment))['network_module_parameters']['min_module_size'] == 37


@pytest.mark.parametrize('existing_provenance', [0, 1, 2])
def test_approved_reader_copy_does_not_claim_current_network_defaults(context, tmp_path, existing_provenance):
    network(context, 'a', 20, .25)
    source = tmp_path / 'approved'
    (source / 'tables').mkdir(parents=True)
    (source / 'figures').mkdir()
    (source / 'figures' / 'image.png').write_bytes(b'figure')
    (source / 'tables' / 'all.tsv').write_text('x\n1\n', encoding='utf-8')
    original = source / 'REPORT_EN.html'
    original.write_text('<html><body><section id="results"><img src="figures/image.png"><a href="tables/all.tsv">All</a></section></body></html>', encoding='utf-8')
    if existing_provenance:
        (source / 'NETWORK_MODULE_PARAMETERS.json').write_text('{"historical": true}', encoding='utf-8')
    if existing_provenance == 2:
        (source / 'NETWORK_MODULE_PARAMETERS_COPY.json').write_text('{"historical_copy": true}', encoding='utf-8')
        original.write_text(original.read_text(encoding='utf-8').replace('</body>',
            '<section id="network-module-parameters">Historical methods</section>'
            '<section id="network-module-parameters-copy">Historical copy methods</section></body>'), encoding='utf-8')
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.rglob('*') if p.is_file()}
    config = {'study': 'Approved', 'approved_reader_root': str(source),
              'manifest_path': str(context.manifest_path)}
    config_path = tmp_path / 'approved-config.json'
    config_path.write_text(json.dumps(config), encoding='utf-8')
    workflow = narrative_report.NarrativeReportWorkflow()
    workflow.run(context, SimpleNamespace(), 'e' * 64, {'config_path': str(config_path)})
    out = context.results.artifact_directory(workflow.definition.workflow_id, 'e' * 64) / 'reader_report'
    text = (out / 'REPORT_EN.html').read_text(encoding='utf-8')
    assert 'Not recorded.' in text
    assert 'Minimum module size: 20' not in text
    assert 'Module merge height: 0.25' not in text
    attachment = ('NETWORK_MODULE_PARAMETERS_COPY_2.json' if existing_provenance == 2 else
                  'NETWORK_MODULE_PARAMETERS_COPY.json' if existing_provenance else 'NETWORK_MODULE_PARAMETERS.json')
    assert json.loads((out / attachment).read_text())['status'] == 'not_recorded'
    if existing_provenance:
        assert (out / 'NETWORK_MODULE_PARAMETERS.json').read_bytes() == (source / 'NETWORK_MODULE_PARAMETERS.json').read_bytes()
    if existing_provenance == 2:
        assert (out / 'NETWORK_MODULE_PARAMETERS_COPY.json').read_bytes() == (source / 'NETWORK_MODULE_PARAMETERS_COPY.json').read_bytes()
        assert 'Historical methods' in text and 'Historical copy methods' in text
    assert before == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before}


def _real_report_render_fixture(base):
    """Supply genuine inputs for both polished heatmap rendering branches."""
    landscape = base / '03_AS_landscape/Native_PSI_clustering/A_vs_B/significant'
    landscape.mkdir(parents=True)
    pd.DataFrame({'A_rep1': [.1, .7], 'A_rep2': [.2, .8], 'B_rep1': [.5, .3]},
        index=pd.Index(['e1', 'e2'], name='event_id')).to_csv(
            landscape / 'display_PSI_clustered.tsv', sep='\t')
    pd.DataFrame({'gene_symbol': ['Gene1', 'Gene2'], 'event_type': ['SE', 'RI']},
        index=pd.Index(['e1', 'e2'], name='event_id')).to_csv(
            landscape / 'event_annotations.tsv', sep='\t')
    (landscape / 'parameters.json').write_text(json.dumps({'eligible_events': 2}))
    selected = base / '12_Gene_cards_heatmaps/gene-run'
    selected.mkdir(parents=True)
    pd.DataFrame({'event_id': ['e1', 'e2'], 'display_event_id': ['SE1', 'RI1'],
                  'A_rep1': [.1, .7], 'A_rep2': [.2, .8], 'B_rep1': [.5, .3]}).to_csv(
                      selected / 'selected_gene_psi.tsv', sep='\t', index=False)
    pd.DataFrame({'event_id': ['e1', 'e2'], 'gene_symbol': ['Gene1', 'Gene2'],
                  'event_type': ['SE', 'RI']}).to_csv(
                      selected / 'selected_gene_events.tsv', sep='\t', index=False)
    (selected / 'parameters.json').write_text(json.dumps({'parameters': {'rows_per_page': 40}}))

    def table(relative, records):
        target = base / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(records).to_csv(target, sep='\t', index=False)

    table('20_Presentation_report/comparison_totals.tsv',
          [{'comparison_id': 'A_vs_B', 'significant_events': 1, 'DEG': 1}])
    table('20_Presentation_report/DEG_DAS_gene_intersection.tsv',
          [{'comparison_id': 'A_vs_B', 'gene_id': 'gene1', 'P': .03, 'FDR': .06}])
    table('22_Expression_interaction_pathways/design/run/interaction.tsv',
          [{'gene_id': 'gene1', 'gene_symbol': 'Gene1', 'padj': .06}])
    table('22_Expression_interaction_pathways/design/run/ranked_pathways.tsv',
          [{'contrast': 'interaction', 'pathway': 'KEGG_EXAMPLE', 'study_family_FDR': .2,
            'NES': -1.2, 'leadingEdge': 'Gene1'}])
    table('16_Interaction/design/run/dexseq_interaction.tsv',
          [{'groupID': 'gene1', 'bin_id': 'E001', 'interaction_FDR': .08, 'pvalue': .04}])
    table('23_Module_explanation/run/module_cards.tsv', [{'module': 'C1', 'n_events': 2}])

    names = [Path('A_vs_B/significant') / f'{mode}.{ext}'
             for mode in ('PSI', 'rowZ') for ext in ('png', 'svg', 'pdf')]
    names += [Path('A_vs_B/significant/display_policy.json')]
    names += [Path('selected_genes/gene-run') / f'selected_gene_heatmap_001.{ext}'
              for ext in ('png', 'svg', 'pdf')]
    for relative in names:
        frozen = base / '24_Narrative_report/figures' / relative
        frozen.parent.mkdir(parents=True, exist_ok=True)
        frozen.write_bytes(b'Previously approved report attachment: ' + str(relative).encode())
    return names


@pytest.mark.parametrize('entry', ['polished_heatmaps', 'build_report'])
def test_real_rendering_keeps_frozen_figures_and_uses_scoped_output(tmp_path, entry):
    from katana_next.science.report_complete import polished_heatmaps
    from katana_next.science.narrative_report import build_report

    base = tmp_path / 'delivery'
    expected = _real_report_render_fixture(base)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in base.rglob('*') if p.is_file()}
    out = base / '24_Narrative_report' / 'new-run'
    if entry == 'polished_heatmaps':
        polished_heatmaps(base, output_dir=out)
    else:
        manifest = tmp_path / 'manifest.json'
        manifest.write_text(json.dumps({'samples': [
            {'sample_id': 'A_rep1', 'condition': 'A', 'biological_sample_id': 'mouse1'}]}))
        build_report('Real-render fixture', base, manifest, out)
        page = (out / 'REPORT_EN.html').read_text(encoding='utf-8')
        assert 'figures/A_vs_B/significant/rowZ.png' in page
    assert before == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before}
    for relative in expected:
        target = out / 'figures' / relative
        assert target.is_file(), target
        if target.suffix == '.png':
            assert target.read_bytes().startswith(b'\x89PNG')
        elif target.suffix == '.pdf':
            assert target.read_bytes().startswith(b'%PDF')
        elif target.suffix == '.svg':
            assert '<svg' in target.read_text(encoding='utf-8')
    policy = json.loads((out / 'figures/A_vs_B/significant/display_policy.json').read_text())
    assert policy['render_only'] is True
    assert policy['eligible_events'] == 2
