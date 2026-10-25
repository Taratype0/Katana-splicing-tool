"""Empty display groups never alter measured values or claim missing analyses."""
import hashlib
import copy
from pathlib import Path
import re
import shutil

import pandas as pd
import pytest

from katana_next.reports.engine import ASSET_ROOT, Builder, CHAPTERS
from katana_next.reports.sources import load_source


def _builder(tmp_path):
    output = tmp_path / 'report'
    output.mkdir()
    shutil.copytree(ASSET_ROOT, output / 'assets')
    project = {'project_id': 'completed-display-fixture', 'name': 'Recorded results', 'samples': [],
               'contrasts': [{'contrast_id': 'B_vs_A', 'numerator': {'label': 'B'}, 'denominator': {'label': 'A'}}]}
    return Builder({'project': project, 'mode': 'bulk', 'artifacts': []}, output, 'stable_significant')


def _chapter_five(builder):
    builder.write()
    html = (builder.output / 'REPORT.html').read_text('utf-8')
    return re.search(r'<section class="chapter" id="chapter-5">(.*?)</section>', html, re.S).group(1)


@pytest.mark.parametrize('state', ['significant', 'nonsignificant', 'empty'])
def test_completed_custom_enrichment_without_go_kegg_panel_is_a_view_status(tmp_path, state):
    native = tmp_path / 'custom-gmt-results.parquet'
    frame = pd.DataFrame({'term': ['custom set'], 'ontology_category': ['Unclassified'],
                          'pvalue': [1e-13 if state == 'significant' else .7],
                          'fdr': [2e-13 if state == 'significant' else .9]})
    if state == 'empty':
        frame = frame.iloc[:0]
    frame.to_parquet(native)
    before = (native.read_bytes(), native.stat().st_mtime_ns)
    b = _builder(tmp_path)
    b.table({'path': str(native), 'analysis_id': 'prioritize.enrichment',
             'contrast_id': 'B_vs_A', 'name': 'Complete recorded custom GMT results'})
    # An independently absent resource must retain its explicit source status.
    b.exclusions.append({'chapter': 5, 'title': 'Separate database not computed',
                         'reason': 'No matching optional database was supplied.'})
    html = _chapter_five(b)
    assert 'No GO/KEGG panel is available for this selection.' in html
    assert 'See complete output tables for recorded enrichment results, input genes and backgrounds.' in html
    assert 'This component is not computed' not in html
    assert 'No matching optional database was supplied.' in html
    assert 'No significant' not in html
    assert not b.panels  # Custom GMT does not become an unapproved GO/KEGG panel.
    saved = pd.read_parquet(b.output / b.tables[0]['source'])
    pd.testing.assert_frame_equal(saved, frame)
    assert before == (native.read_bytes(), native.stat().st_mtime_ns)
    assert all((b.output / 'assets' / name).read_bytes() == (ASSET_ROOT / name).read_bytes()
               for name in ('report.css', 'report.js'))
    assert len(CHAPTERS) == 7


@pytest.mark.parametrize('attachment', [None, 'gene', 'term_mapping'])
def test_absent_enrichment_statistics_keep_not_computed_status(tmp_path, attachment):
    b = _builder(tmp_path)
    if attachment:
        genes = tmp_path / 'input-genes.tsv'
        genes.write_text('gene_id\ng1\n' if attachment == 'gene' else
                         'term\tinput_identifiers\tmapped_unique_genes\nset1\t2\t2\n', encoding='utf-8')
        b.table({'path': str(genes), 'analysis_id': 'prioritize.enrichment', 'chapter': 5, 'name': 'Input genes only'})
    html = _chapter_five(b)
    assert 'No completed result available' in html and 'This component is not computed.' in html
    assert 'No GO/KEGG panel is available' not in html


def test_same_named_enrichment_views_preserve_four_foregrounds_and_native_comparisons(tmp_path):
    records = []
    for foreground in ('das', 'deg', 'dtu', 'dexseq'):
        for comparison in ('B_vs_A', 'C_vs_A'):
            for kind, entity in (('plot_spec', 'gene_set_summary'), ('table', 'gene_set')):
                path = tmp_path / f'{foreground}-{comparison}-{kind}.tsv'
                path.write_text('term\tpvalue\tfdr\nset1\t1e-13\t2e-13\n', encoding='utf-8')
                records.append({'analysis_id': 'prioritize.enrichment', 'name': 'Gene-set enrichment',
                    'kind': kind, 'entity_type': entity, 'path': str(path), 'status': 'complete',
                    'parameters': {'foreground_source': foreground, 'event_set': 'stable_significant'},
                    'contrast_ids': [], 'metadata': {'native_background_provenance': {
                        'native_inputs': [{'contrast_id': comparison}, {'contrast_id': comparison}]}}})
    original = copy.deepcopy(records)
    paths = {row['path'] for row in records}
    before = {path: (Path(path).read_bytes(), Path(path).stat().st_mtime_ns) for path in paths}
    loaded = load_source({'project': {'project_id': 'display-study', 'samples': [], 'contrasts': []},
                          'artifacts': records + [copy.deepcopy(records[0])]})['artifacts']
    assert len(loaded) == 16  # Two comparisons, four evidence families, plot + complete statistics.
    assert {r['parameters']['foreground_source'] for r in loaded if r['kind'] == 'table'} == {'das', 'deg', 'dtu', 'dexseq'}
    assert len({r['path'] for r in loaded}) == 16
    assert records == original
    assert before == {path: (Path(path).read_bytes(), Path(path).stat().st_mtime_ns) for path in paths}


def test_zero_pie_categories_hide_only_their_labels_and_keep_complete_counts(tmp_path, monkeypatch):
    from matplotlib.colors import to_rgba
    import matplotlib.pyplot as plt
    from katana_next.reports import comparisons
    from test_report_comparisons import _contrast, _events
    native = tmp_path / 'measured-events.parquet'
    original = _events()
    original.to_parquet(native)
    before = hashlib.sha256(native.read_bytes()).hexdigest()
    observed = {}
    def capture(self, fig, identity, chapter, title, comparison, caption, policy=None):
        if '-AS-composition-' in identity:
            axis = fig.axes[0]
            observed[title] = {'texts': [text.get_text() for text in axis.texts],
                               'colors': [patch.get_facecolor() for patch in axis.patches]}
        self.panels.append({'id': identity, 'chapter': chapter, 'title': title,
                            'contrast_id': comparison, 'policy': policy})
        plt.close(fig)
    monkeypatch.setattr(comparisons._Views, 'figure', capture)
    result = comparisons.build_comparison_views({'contrasts': [_contrast()], 'assets': []},
        [{'path': str(native), 'analysis_id': 'differential.as', 'name': 'Differential AS events', 'parameters': {}}],
        tmp_path / 'views')
    assert len(observed) == 2
    for record in observed.values():
        assert [text for text in record['texts'] if text] == ['SE', '100.0%']
        assert record['colors'] == [to_rgba(color) for color in comparisons.EVENT_COLORS]
    for panel in result['panels']:
        if '-AS-composition-' not in panel['id']:
            continue
        plotted = pd.read_csv(tmp_path / 'views' / panel['data_links'][0]['path'], sep='\t')
        assert plotted.event_type.tolist() == list(comparisons.EVENTS)
        assert plotted.events.tolist()[1:] == [0, 0, 0, 0]
        assert plotted.percent.tolist() == [100., 0., 0., 0., 0.]
        assert plotted.display_color.tolist() == comparisons.EVENT_COLORS
    assert result['new_scientific_fits'] == 0
    assert hashlib.sha256(native.read_bytes()).hexdigest() == before
