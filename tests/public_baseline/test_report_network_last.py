"""The approved shared network closes the report without relabelling old data."""
import json
import shutil
from hashlib import sha256
from pathlib import Path

from bs4 import BeautifulSoup

from katana_next.reports.engine import ASSET_ROOT, Builder
from katana_next.reports.refresh import refresh_report


def test_whole_study_network_is_last_and_shared_across_comparisons(tmp_path):
    output = tmp_path / 'report'
    shutil.copytree(ASSET_ROOT, output / 'assets')
    source = {'project': {'project_id': 'shared-network', 'samples': [], 'contrasts': [
        {'contrast_id': c, 'numerator': {'label': c}, 'denominator': {'label': 'A'}}
        for c in ('B_vs_A', 'C_vs_A')]}, 'artifacts': [], 'mode': 'bulk'}
    builder = Builder(source, output, 'stable_significant')
    network = output / 'outputs' / 'shared-network.html'
    network.parent.mkdir()
    network.write_text('<html><h1>Saved whole-study modules</h1><p>Complete genes</p></html>')
    original_hash = sha256(network.read_bytes()).hexdigest()
    builder.tables.append({'name': 'Shared network', 'chapter': 6,
        'analysis_id': 'rich_browser', 'contrast_id': 'all',
        'source': 'outputs/shared-network.html', 'metadata': {'report_role': 'whole_study_chapter'}})
    manifest = builder.write()
    page = BeautifulSoup((output / 'REPORT.html').read_text('utf8'), 'html.parser')
    chapters = page.select('section.chapter:not(#appendix)')
    assert chapters[-2].h2.text == 'RBP and MotifTrace'
    assert chapters[-1].h2.text == 'Katana-Splice-Network'
    assert [s.select_one('.chapter-number').text for s in chapters] == [f'{i:02}' for i in range(1, 8)]
    assert page.select('.toc a')[-2].get_text() == '07Katana-Splice-Network'
    browser = chapters[-1].select_one('.network-browser')
    assert browser['data-contrast'] == 'all' and browser['data-policy'] == 'common'
    assert 'does not refit or subset this network' in browser.get_text()
    assert manifest['presentation']['chapter_order'] == [1, 2, 3, 4, 5, 7, 6]
    assert manifest['network_browsers'][0]['whole_study'] is True
    index = BeautifulSoup((output / 'outputs/index.html').read_text('utf8'), 'html.parser')
    assert index.select('section h2')[-1].text == '7. Katana-Splice-Network'
    refreshed = tmp_path / 'refreshed'
    result = refresh_report(output, refreshed)
    assert result['network_browsers'][0]['whole_study'] is True
    assert sha256(network.read_bytes()).hexdigest() == original_hash
    assert sha256((refreshed / 'outputs/shared-network.html').read_bytes()).hexdigest() == original_hash


def test_network_opt_out_does_not_reappear_when_chapters_reordered(tmp_path):
    output = tmp_path / 'report'
    shutil.copytree(ASSET_ROOT, output / 'assets')
    source = {'project': {'project_id': 'optional-network', 'samples': [], 'contrasts': []},
              'artifacts': [], 'mode': 'bulk', 'report_analyses': {'network': False}}
    builder = Builder(source, output, 'stable_significant')
    builder.write()
    page = BeautifulSoup((output / 'REPORT.html').read_text('utf8'), 'html.parser')
    assert page.select_one('#chapter-6') is None
    assert page.select('section.chapter:not(#appendix)')[-1].h2.text == 'RBP and MotifTrace'


def test_motif_annotations_are_packaged_and_statistical_endpoints_stay_separate():
    from importlib.resources import files
    from katana_motiftrace.annotations import candidate_genes
    from katana_next.reports.motif_sequences import evidence
    reference = json.loads(files('katana_motiftrace').joinpath('motif_gene_reference.json').read_text('utf8'))
    assert reference['records']
    models = {'demo': {'kind': 'PWM', 'matrix': [[.7, .1, .1, .1]] * 5,
                       'assignments': [{'RBP': 'PTBP1', 'species': 'Homo sapiens'}]}}
    assert candidate_genes(models['demo'], 'homo_sapiens', reference)
    report = evidence(models, [{'model_id': 'demo', 'cluster_P_FWER_new_family': .003,
        'layer': 'affected_sequence', 'region': 'affected_inclusion',
        'event_type': 'SE', 'direction': 'inclusion_up'}], [], regional_available=True)
    assert report['displayed_models'] == 1
    assert report['window_significant_models'] == 0
    assert report['regional_significant_regions'] == 1
