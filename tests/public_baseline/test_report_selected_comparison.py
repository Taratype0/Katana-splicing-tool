"""Report refresh changes presentation, while native files remain immutable."""
from hashlib import sha256
import json
from pathlib import Path

from PIL import Image
import pytest

from katana_next.reports.engine import ASSET_ROOT, Builder
from katana_next.reports.refresh import refresh_report


def builder(tmp_path, name='report'):
    output = tmp_path / name
    output.mkdir()
    import shutil
    shutil.copytree(ASSET_ROOT, output / 'assets')
    project = {'project_id': 'shared-project', 'name': 'Display check', 'samples': [],
               'contrasts': [{'contrast_id': cid, 'numerator': {'label': 'B'},
                              'denominator': {'label': 'A'}} for cid in ['B_vs_A', 'C_vs_A']]}
    return Builder({'project': project, 'artifacts': [], 'mode': 'bulk'}, output, 'stable_significant')


def test_report_shows_selected_comparison_network_and_keeps_dendrogram_in_outputs(tmp_path):
    b = builder(tmp_path)
    browser = b.output / 'outputs/rich/browser.html'
    browser.parent.mkdir(parents=True)
    browser.write_text('<html><h1>Network</h1><p>Complete genes</p></html>')
    b.tables.append({'name': 'All network modules: example', 'source': browser.relative_to(b.output).as_posix(),
                     'chapter': 6, 'analysis_id': 'rich_browser', 'contrast_id': 'B_vs_A',
                     'metadata': {'network_view_role': 'module_browser'}})
    image = b.output / 'tree.png'
    Image.new('RGB', (32, 32), 'white').save(image)
    b.panels.append({'id': 'tree', 'title': 'Observed network dendrogram', 'chapter': 6,
                     'image': 'tree.png', 'contrast_id': 'B_vs_A', 'metadata': {}})
    m = b.write()
    page = (b.output / 'REPORT.html').read_text('utf8')
    assert 'All comparisons</option>' not in page
    assert '<iframe loading="lazy" data-src="outputs/rich/browser.html"' in page
    assert 'id="tree"' not in page
    assert m['panels'][0]['output_only'] is True
    assert 'Observed network dendrogram' in (b.output / 'outputs/index.html').read_text('utf8')
    assert m['presentation']['default_comparison'] == 'B_vs_A'
    assert len(m['network_browsers']) == 1
    assert m['display_project']['contrasts'] == b.source['project']['contrasts']


def test_refresh_is_public_and_keeps_previous_report_and_native_values(tmp_path):
    b = builder(tmp_path)
    native = b.output / 'outputs/native.tsv'
    native.parent.mkdir()
    native.write_text('event_id\tfdr\ndisplay-event\t0.5925\n')
    b.tables.append({'name': 'Native statistical values', 'source': 'outputs/native.tsv', 'chapter': 7})
    b.source['artifacts'] = [{'analysis_id': 'prioritize.motiftrace_consequences', 'artifact_id': 'native', 'parameters': {'native_setting': 17}}]
    b.write()
    before = {p.relative_to(b.output).as_posix(): sha256(p.read_bytes()).hexdigest()
              for p in b.output.rglob('*') if p.is_file()}
    out = tmp_path / 'refreshed'
    m = refresh_report(b.output, out)
    assert (out / 'REPORT.zip').is_file()
    assert (out / 'outputs/native.tsv').read_bytes() == native.read_bytes()
    assert m['contrasts'] == ['B_vs_A', 'C_vs_A']
    receipt = json.loads((out / 'presentation-refresh.json').read_text('utf8'))
    assert json.loads((out/'provenance.json').read_text())['source_projects'] == json.loads((b.output/'provenance.json').read_text())['source_projects']
    assert receipt['original_report_unchanged'] is True and receipt['completed_analyses_rerun'] == 0
    assert all(sha256((b.output / p).read_bytes()).hexdigest() == v for p, v in before.items())
    with pytest.raises(FileExistsError): refresh_report(b.output, out)
    with pytest.raises(ValueError): refresh_report(b.output, b.output / 'overwrite-source')
    with pytest.raises(ValueError): refresh_report(b.output, tmp_path / 'foreign', project={'project_id': 'foreign'})
    changed_design = json.loads(json.dumps(b.source['project']))
    changed_design['contrasts'][0]['denominator']['sample_ids'] = ['substitute-biological-unit']
    with pytest.raises(ValueError, match='contrast arms'):
        refresh_report(b.output, tmp_path / 'changed-design', project=changed_design)


def test_group_scope_is_explicit_and_reset_returns_to_first_comparison(tmp_path):
    from katana_next.reports.layout import _scope
    assert 'data-comparisons="[&quot;B_vs_A&quot;]"' in _scope({'comparison_ids': ['B_vs_A']})
    js = (ASSET_ROOT / 'report.js').read_text()
    assert 'group.includes(contrast)' in js
    assert "$('contrastSelect').value=initialComparison" in js
    assert "if(!show&&frame.getAttribute('src'))frame.removeAttribute('src')" in js
