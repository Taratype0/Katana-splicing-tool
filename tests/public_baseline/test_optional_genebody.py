"""Normal QC remains on without implicitly scheduling an extra BAM scan."""
import copy

from katana_next.flow.model import validate_plan
from katana_next.flow.planner import create_plan, preflight_plan
from katana_next.flow.selection import active_chapters, analysis_selection, component_enabled
from test_flow_planner import study


def test_default_plan_does_not_require_extra_rseqc_or_bams(study):
    plan = create_plan(study, {})
    enabled = plan['configuration']['analyses']
    assert enabled['sequencing_qc'] is True and enabled['coverage'] is False
    assert active_chapters(enabled) == set(range(1, 8))
    assert not component_enabled('gene_body_coverage', enabled)
    assert any(t.get('adapter') == 'qc.sequencing' for t in plan['tasks'])
    assert not any(t.get('adapter') in {'qc.genebody', 'reference.bed12'} for t in plan['tasks'])
    assert not any(t['workflow_id'] == 'quality.read_coverage_bias' for t in plan['tasks'])
    assert not any('genebody' in e.get('task_id', '') or 'Gene-body' in e.get('error', '')
                   for e in preflight_plan(plan)['errors'])


def test_explicit_coverage_still_schedules_one_scan_per_sample(study):
    from katana_next.infrastructure.manifest_store import ManifestStore
    _, manifest = ManifestStore().load(study)
    config = {'analyses': {'coverage': True}}
    original = copy.deepcopy(config)
    plan = create_plan(study, config)
    assert config == original
    assert len([t for t in plan['tasks'] if t.get('adapter') == 'qc.genebody']) == len(manifest.samples)
    assert any(t['id'] == 'reference-bed' for t in plan['tasks'])
    assert any(t['id'] == 'coverage-summary' for t in plan['tasks'])
    # Explicit selections in a sealed plan survive the changed new-project default.
    validate_plan(plan)
    assert analysis_selection(plan['configuration'])['coverage'] is True
