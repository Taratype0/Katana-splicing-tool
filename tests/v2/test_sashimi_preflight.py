from pathlib import Path

import pandas as pd
import pytest

from katana_next.client import KatanaProject
from katana_next.domain.models import SelectionSet
from katana_next.sashimi_fixture import create_sashimi_fixture


@pytest.mark.parametrize('missing', ['bam', 'index'])
def test_count_cache_cannot_make_bam_preflight_ready(tmp_path, missing):
    manifest = create_sashimi_fixture(tmp_path/'fixture')
    with KatanaProject(manifest) as project:
        project.run('single_cell.quantify', {'source':'bam','min_cell_event_reads':1})
        source = next(a for a in project.artifacts('single_cell.quantify') if a.name == 'Single-cell junction choices')
        row = pd.read_parquet(source.path).row_id.iloc[0]
        selection = SelectionSet(project_id=project.project_id, name='preflight test',
            source_artifact_id=source.artifact_id, entity_type=source.entity_type, entity_ids=(str(row),))
        project.projects.save_selection(selection)
        parameters = {'selection_id':selection.selection_id}
        assert project.plan('single_cell.sashimi', parameters)['ready']
        sample = project.context.manifest.samples[0]
        bam = sample.assets.get('single_cell_bam') or sample.assets.get('bam')
        removed = bam if missing == 'bam' else Path(str(bam)+'.bai')
        backup = removed.with_suffix(removed.suffix+'.test-disabled')
        removed.rename(backup)
        try:
            before = len(project.jobs.list())
            plan = project.plan('single_cell.sashimi', parameters)
            assert not plan['ready']
            assert any(('BAM' in e or 'source' in e.lower() or 'snapshot' in e.lower()) for e in plan['errors'])
            assert len(project.jobs.list()) == before
        finally:
            backup.rename(removed)


def test_sashimi_preflight_requires_saved_selection(tmp_path):
    manifest = create_sashimi_fixture(tmp_path/'fixture')
    with KatanaProject(manifest) as project:
        plan = project.plan('single_cell.sashimi')
        assert not plan['ready'] and any('selection' in error for error in plan['errors'])
        assert not project.jobs.list()
