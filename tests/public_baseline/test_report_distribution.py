"""Installed report templates and console command work outside the checkout."""
from importlib import metadata, resources
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
import zipfile


def test_report_resources_are_in_the_distribution():
    for name, token in [('report.css', '--blue'), ('report.js', 'reportEventSet'),
                        ('table-template.html', '__DATA__'), ('matrix-browser.html', '__DATA__')]:
        assert token in resources.files('katana_next.reports').joinpath('assets', name).read_text('utf-8')
    entry = next(e for e in metadata.distribution('katana-splicing-tool').entry_points if e.name == 'katana-report')
    assert entry.value == 'katana_next.reports.cli:main'


def test_installed_report_cli_from_an_unrelated_directory(tmp_path):
    import katana_next.reports
    if os.environ.get('KATANA_WHEEL_ONLY') == '1':
        assert Path(katana_next.reports.__file__).resolve().is_relative_to(Path(sysconfig.get_path('purelib')).resolve())
    project = tmp_path/'source.json'
    project.write_text(json.dumps({'project': {'project_id': 'installed-report-fixture', 'name': 'Completed fixture',
                                              'output_root': str(tmp_path/'original'), 'samples': [], 'contrasts': []},
                                   'artifacts': []}), encoding='utf-8')
    args = [sys.executable, *(['-I'] if os.environ.get('KATANA_WHEEL_ONLY') == '1' else []),
            '-m', 'katana_next.reports.cli', '--project', str(project), '--outdir', str(tmp_path/'report')]
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'MPLBACKEND': 'Agg'}
    if os.environ.get('KATANA_WHEEL_ONLY') != '1':
        env['PYTHONPATH'] = str(Path(katana_next.reports.__file__).resolve().parents[2])
    process = subprocess.run(args, cwd=tmp_path, capture_output=True, text=True, env=env, timeout=120)
    assert process.returncode == 0, process.stderr
    result = json.loads(process.stdout)
    assert result['status'] == 'completed' and result['no_upstream_rerun']
    page = Path(result['report']).read_text('utf-8')
    assert all(f'id="chapter-{n}"' in page for n in range(1, 8))
    assert 'No completed result' in page
    with zipfile.ZipFile(result['archive']) as archive:
        assert archive.testzip() is None and 'assets/matrix-browser.html' in archive.namelist()
