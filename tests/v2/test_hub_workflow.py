from pathlib import Path
import json
import numpy as np
import pandas as pd
from fastapi.testclient import TestClient

from katana_next.client import KatanaProject
from katana_next.full_fixture import create_full_fixture
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows.registry import WorkflowRegistry
from katana_next.api.app import PlatformRuntime, create_app
from test_hub_events import planted, RSCRIPT


def real_R_fixture(path):
    manifest=create_full_fixture(path)
    value=json.loads(manifest.read_text(encoding='utf-8'))
    value.setdefault('metadata',{})['rscript_path']=str(RSCRIPT)
    manifest.write_text(json.dumps(value,indent=2),encoding='utf-8')
    return manifest


def make_discovery(project, matrix, assignments, samples, events, *, edge_cap=False):
    writer=ArtifactWriter(project.context,'programs.cosplice','fixed-fixture-discovery',{},WorkflowRegistry().get('programs.cosplice').version)
    parent=writer.table('Consensus assignments','consensus_modules__assignments.parquet',assignments.rename(columns={'module':'consensus_module'}),entity_type='module')
    writer.table('Measured PSI','matrix__psi.parquet',matrix.rename_axis('event_key').reset_index(),entity_type='matrix')
    writer.table('Event annotation','matrix__events.parquet',events,entity_type='event')
    writer.table('Discovery samples','matrix__samples.parquet',samples,entity_type='sample')
    if edge_cap:
        writer.table('Truncated edges','network_observed__edges.parquet',pd.DataFrame({'event_a':['e1'],'event_b':['e2'],'adjacency':[.8],'correlation':[.9]}),entity_type='module')
    config=writer.directory/'resolved_cosplice.yaml'
    import yaml
    config.write_text(yaml.safe_dump({'network':{'method':'pearson','power':6,'signed':True,'edge_threshold':.1,'wgcna_rscript':str(RSCRIPT)}}),encoding='utf-8')
    writer.existing_file('Resolved configuration',config,entity_type='configuration')
    writer.finish()
    return parent


def test_end_to_end_immutable_hub_job_cache_dashboard_and_API(tmp_path, monkeypatch):
    manifest=real_R_fixture(tmp_path/'project')
    matrix,assignments,samples,events=planted()
    with KatanaProject(manifest) as project:
        parent=make_discovery(project,matrix,assignments,samples,events,edge_cap=True)
        digest=parent.path.read_bytes()
        before=project.module_dashboard(parent.artifact_id,'M1')
        assert before['hub_events']==[] and before['hub_status']=='not_computed'
        parameters={'source_artifact_id':parent.artifact_id,'bootstrap_iterations':20,'top_fraction':.2,'minimum_selection_frequency':.5,'apply_threshold_labels':True}
        assert project.plan('programs.hub_events',parameters)['ready']
        result=project.run('programs.hub_events',parameters)
        assert result.status=='complete'
        assert project.run('programs.hub_events',parameters).status=='cache_hit'
        assert parent.path.read_bytes()==digest
        outputs=project.artifacts('programs.hub_events')
        metrics=next(a for a in outputs if a.path.stem=='hub_events__metrics')
        rows=pd.read_parquet(metrics.path)
        assert len(rows)==len(matrix)*3 and rows.loc[rows.event_key.eq('e0'),'hub_candidate'].all()
        assert parent.artifact_id in metrics.parent_artifact_ids
        assert metrics.metadata['source_discovery_run']==parent.run_fingerprint
        assert metrics.path.with_suffix('.tsv').is_file()
        view=project.module_dashboard(metrics.artifact_id,'M1')
        assert view['hub_status']=='available' and view['hub_events'][0]['event_key']=='e0'
        assert view['hub_events'][0]['connectivity']>1 # The one exported edge did not define hub connectivity.
        assert len(view['hub_comparison'])==len(matrix)
        hidden=project.module_dashboard(metrics.artifact_id,'M1',hide_same_gene=True)
        assert hidden['hub_events']==view['hub_events'] # Explicit side-by-side ranks exist independently of drawing toggle.
        assert 'rank_no_same_gene' in hidden['hub_comparison'][0]
        # The same science with labels off must reuse its completed native output.
        # In particular, toggling a descriptive label must not rerun R/bootstrap.
        old_metrics_bytes = metrics.path.read_bytes()
        monkeypatch.setattr('katana_next.workflows.hub_events.run_hub_analysis',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('Label switch reran the scientific engine')))
        unlabelled = project.run('programs.hub_events', {**parameters, 'apply_threshold_labels': False})
        assert unlabelled.status == 'complete', unlabelled.error
        off_metrics = next(a for a in project.artifacts('programs.hub_events')
                           if a.run_fingerprint == unlabelled.run_fingerprint and a.path.stem == 'hub_events__metrics')
        off_rows = pd.read_parquet(off_metrics.path)
        assert len(off_rows) == len(rows)
        assert off_rows.hub_candidate.isna().all()
        assert off_rows.label_status.eq('not_applied').all()
        keys = ['event_key', 'module', 'condition']
        continuous = ['kWithin', 'kWithin_no_same_gene', 'kME', 'rank_all_edges', 'rank_no_same_gene']
        pd.testing.assert_frame_equal(rows.set_index(keys)[continuous].sort_index(),
                                      off_rows.set_index(keys)[continuous].sort_index())
        assert metrics.path.read_bytes() == old_metrics_bytes
        off_view = project.module_dashboard(off_metrics.artifact_id, 'M1')
        assert off_view['hub_status'] == 'labels_not_applied'
        assert len(off_view['hub_comparison']) == len(matrix)
        assert off_view['hub_events'] == []  # No candidate labels were requested.
        assert off_view['hub_label_status'] == 'not_applied'
        assert {j.analysis_id for j in project.jobs.list()}=={'programs.hub_events'}
        project_id=project.project_id
        source_id=metrics.artifact_id
        off_source_id=off_metrics.artifact_id
    runtime=PlatformRuntime()
    with TestClient(create_app(runtime)) as client:
        opened=client.post('/api/projects/open',json={'manifest_path':str(manifest)})
        assert opened.status_code==200
        response=client.get(f'/api/projects/{project_id}/artifacts/{source_id}/module-dashboard',params={'module':'M1','condition':'A'})
        assert response.status_code==200,response.text
        assert response.json()['hub_events'][0]['event_key']=='e0'
        assert response.json()['hub_analysis_run']
        unlabelled_response=client.get(f'/api/projects/{project_id}/artifacts/{off_source_id}/module-dashboard',params={'module':'M1','condition':'all'})
        assert unlabelled_response.status_code==200,unlabelled_response.text
        assert unlabelled_response.json()['hub_status']=='labels_not_applied'
        assert len(unlabelled_response.json()['hub_comparison'])==len(matrix)


def test_no_edge_analysis_exports_reasons_and_returns_no_hubs(tmp_path):
    manifest=real_R_fixture(tmp_path/'project')
    matrix,a,s,e=planted();matrix.loc[:,:]=.5
    with KatanaProject(manifest) as project:
        parent=make_discovery(project,matrix,a,s,e)
        project.run('programs.hub_events',{'source_artifact_id':parent.artifact_id,'bootstrap_iterations':20,'condition_scope':'all','apply_threshold_labels':True})
        view=project.module_dashboard(parent.artifact_id,'M1')
        assert view['hub_events']==[] and view['hub_status']=='no_supported_edges'
        assert len(view['hub_comparison'])==len(matrix)
        assert all('no_supported_edges' in row['selection_reason'] for row in view['hub_comparison'])
