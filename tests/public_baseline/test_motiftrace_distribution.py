"""MotifTrace public distribution and real installed-resource smoke checks.

The release acceptance also executes this file from an unrelated directory in a
fresh wheel-only environment, with KATANA_WHEEL_ONLY=1 and isolated Python mode.
These checks do not fit external tools or validate a complete scientific study.
"""
from importlib import metadata, resources
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig

import pandas as pd
import pytest


@pytest.mark.parametrize('name, marker', [
    ('rbp_sequence_browser.html', '__RBP_DATA__'),
    ('motiftrace_windows.html', '__WINDOW_DATA__'),
    ('motiftrace_consequences.html', '__MOTIFTRACE_CONSEQUENCES__'),
    ('motiftrace_parallel.html', '__PARALLEL_PAYLOAD__'),
])
def test_report_templates_are_real_package_resources(name, marker):
    text = resources.files('katana_next.science').joinpath(name).read_text(encoding='utf-8')
    assert marker in text and '<html' in text.lower()


def test_installed_namespace_and_cli_inventory_do_not_initialize_gui_or_server():
    import katana_motiftrace
    from katana_motiftrace import mode_catalog

    if os.environ.get('KATANA_WHEEL_ONLY') == '1':
        installed = Path(katana_motiftrace.__file__).resolve()
        assert installed.is_relative_to(Path(sysconfig.get_path('purelib')).resolve())
        distribution = metadata.distribution('katana-splicing-tool')
        entry = next(item for item in distribution.entry_points if item.name == 'katana-motiftrace')
        assert entry.value == 'katana_motiftrace.cli:main'
    catalog = mode_catalog()
    assert catalog['package_version'] == katana_motiftrace.__version__ == '1.0.0.post1'
    assert catalog['scientific_release'] == katana_motiftrace.SCIENTIFIC_RELEASE == '1.0.0-frozen.20261017'
    assert {item['mode_id'] for item in catalog['modes']} == {
        'boundary', 'se_impact', 'completed_maps', 'consequences', 'usage', 'parallel'}
    command = [sys.executable, *(['-I'] if os.environ.get('KATANA_WHEEL_ONLY') == '1' else [])]
    module = subprocess.run([*command, '-m', 'katana_motiftrace', 'modes'],
                            text=True, capture_output=True, check=True)
    assert json.loads(module.stdout) == catalog
    child = subprocess.run([*command, '-c',
        "import json,sys; import katana_motiftrace; from katana_next.client import KatanaProject; "
        "print(json.dumps([m for m in sys.modules if m.startswith(('PySide6','PyQt','fastapi','uvicorn'))]))"],
        text=True, capture_output=True, check=True)
    assert json.loads(child.stdout) == []


@pytest.mark.parametrize('renderer', ['browser', 'windows', 'parallel', 'consequences'])
def test_installed_templates_render_without_a_checkout_or_external_tools(tmp_path, renderer):
    output = tmp_path / f'{renderer}.html'
    payload = {'p': 1e-12, 'not_computed': None, 'text': '</script><unsafe>'}
    if renderer == 'browser':
        from katana_next.science.rbp_sequence_results import render_browser
        render_browser(payload, output)
        marker = '__RBP_DATA__'
    elif renderer == 'windows':
        from katana_next.science.motiftrace_windows import render_windows
        render_windows(payload, output)
        marker = '__WINDOW_DATA__'
    elif renderer == 'parallel':
        from katana_next.science.motiftrace_parallel import render_parallel
        render_parallel(payload, output)
        marker = '__PARALLEL_PAYLOAD__'
    else:
        from katana_next.science.motiftrace_consequences import render_consequences
        empty = pd.DataFrame()
        render_consequences(empty, empty, empty, empty, payload, output)
        marker = '__MOTIFTRACE_CONSEQUENCES__'
    page = output.read_text(encoding='utf-8')
    assert marker not in page and '<html' in page.lower()
    assert '1e-12' in page and 'null' in page
    assert '</script><unsafe>' not in page


def test_public_project_inventory_and_plan_are_available_from_the_distribution(tmp_path):
    from katana_motiftrace import MotifTrace
    from katana_next.full_fixture import create_full_fixture

    manifest = create_full_fixture(tmp_path / 'synthetic-project')
    with MotifTrace(manifest) as motiftrace:
        catalog = motiftrace.inventory()
        assert len(catalog['modes']) == 6
        assert all(item['workflow']['workflow_id'] == item['workflow_id'] for item in catalog['modes'])
        plan = motiftrace.plan('parallel')
        assert plan['motiftrace']['mode_id'] == 'parallel'
        assert motiftrace.project.jobs.list(motiftrace.project.project_id) == []
    from katana_next.client import KatanaProject
    with KatanaProject(manifest) as project:
        borrowed = project.motiftrace
        assert borrowed.project is project
        assert len(borrowed.inventory()['modes']) == 6
        borrowed.close()
        assert project.inventory()  # Closing the borrowed facade leaves its owner usable.


def test_shared_domain_and_adapter_modules_are_distributed():
    from src.domain.models.gene_filter import GeneFilterConfig
    from importlib.util import find_spec

    assert GeneFilterConfig is not None
    assert find_spec('src.adapters.native_sashimi_adapter') is not None


def test_boundary_static_plot_exports_are_available_without_a_checkout(tmp_path):
    from katana_next.workflows.figure_exports import export_vega

    specification = {
        '$schema': 'https://vega.github.io/schema/vega-lite/v5.json',
        'data': {'values': [
            {'region': 'upstream', 'count': 2},
            {'region': 'exonic', 'count': 1},
        ]},
        'mark': 'bar',
        'encoding': {
            'x': {'field': 'region', 'type': 'nominal'},
            'y': {'field': 'count', 'type': 'quantitative'},
        },
    }
    exported = export_vega(specification, tmp_path / 'plots', 'motif-boundary-smoke')
    png = Path(exported['png']).read_bytes()
    svg = Path(exported['svg']).read_text(encoding='utf-8')
    assert png.startswith(b'\x89PNG\r\n\x1a\n') and len(png) > 1000
    assert '<svg' in svg and 'upstream' in svg and 'exonic' in svg
