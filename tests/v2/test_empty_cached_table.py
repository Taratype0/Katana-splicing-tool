import pandas as pd
import pytest

from katana_next.client import KatanaProject
from katana_next.full_fixture import create_full_fixture
from katana_next.infrastructure.result_store import ResultStore
from katana_next.workflows.base import ArtifactWriter


def test_legacy_schema_less_not_run_table_can_be_previewed(tmp_path):
    path = tmp_path/'legacy.parquet'
    pd.DataFrame().to_parquet(path)
    store = ResultStore(tmp_path)
    assert store.table_page(path, search='anything', offset=-2, limit=1000) == {
        'columns':[], 'rows':[], 'total':0, 'offset':0, 'limit':500}
    with pytest.raises(ValueError, match='Unknown sort column'):
        store.table_page(path, sort_by='not_real')


def test_new_empty_artifacts_have_readable_table_and_export_headers(tmp_path):
    manifest = create_full_fixture(tmp_path/'fixture')
    with KatanaProject(manifest) as project:
        writer = ArtifactWriter(project.context, 'test.empty', 'empty-acceptance', {}, '1')
        artifact = writer.table('Not-run empty members', 'empty.parquet', pd.DataFrame(), entity_type='diagnostic')
        writer.finish()
        result = project.context.results.table_page(artifact.path)
        assert result['columns'] == ['row_id'] and result['total'] == 0
        assert artifact.path.with_suffix('.tsv').read_text().strip() == 'row_id'
        assert not project.jobs.list()
