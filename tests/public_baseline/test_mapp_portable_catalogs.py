"""Fresh native-result collections need usable, exact full-fit report links."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import zipfile

import pytest

from katana_next.reports.mapp_catalogs import prepare_portable_mapp_panel, validate_portable_mapp_panel
from katana_next.reports.sources import digest


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    return path


def _family(tmp_path, *, count=3, cid='B_vs_A', collection=None):
    collection = collection or tmp_path/'collection'
    rows = []
    for number in range(count):
        task = f'{cid}-fit-{number}'
        path = collection/'products'/task/'AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f'model\tcombined.Zscore\tcorr_pval_A1\nAAA\t{number}.25\t1e-200\n'.encode())
        rows.append({'task_id':task, 'path':path.relative_to(collection).as_posix(), 'sha256':digest(path),
                     'scope':{'contrast_ids':[cid], 'site':'3ss', 'window_start':number, 'window_end':number+20}})
    index = _write(collection/'presentation'/cid/'complete_native_fits.json', {
        'comparison_id':cid, 'fits':rows, 'correction':'Native MAEI per-window/sample Bonferroni; no pooled cross-window FDR',
        'calibrated_motiftrace_parallel':False})
    panel = {'id':'native-maei-'+cid, 'chapter':7, 'contrast_id':cid, 'metadata':{
        'component':'mapp_native', 'analysis_id':'prioritize.mapp_native', 'fit_task_ids':[r['task_id'] for r in rows],
        'native_fit_index_sha256':digest(index), 'calibrated_motiftrace_parallel':False},
        'data_links':[{'role':'full_results', 'path':str(index), 'source_collection_root':str(collection.resolve()),
                       'label':'Every completed native fit and original statistics', 'contrast_id':cid,
                       'policy':None, 'unit':'native motif/window fit'}]}
    return panel,collection,index


def _copy_into(report, *, prefix_override=None):
    def copy(source, *, prefix='results'):
        target = report/'outputs'/(prefix_override or prefix)/(digest(source)[:12]+'-'+source.name)
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copy2(source,target)
        return target.relative_to(report).as_posix()
    return copy


def _bind(panel, collection, report, *, copy=None):
    report.mkdir(exist_ok=True,parents=True)
    return prepare_portable_mapp_panel(panel,report_root=report,copy_file=copy or _copy_into(report))


def test_complete_portable_catalog_retains_original_family_values_and_source_bytes(tmp_path):
    panel,collection,index = _family(tmp_path)
    before = {str(p):(p.read_bytes(),p.stat().st_mtime_ns) for p in collection.rglob('*') if p.is_file()}
    original = deepcopy(panel)
    report = tmp_path/'report'
    result,audit = _bind(panel,collection,report,copy=_copy_into(report,prefix_override='already-copied'))
    assert panel == original
    assert audit['passed'] and audit['fit_count'] == 3
    assert audit['native_statistics_unchanged'] and audit['new_fits'] == 0
    portable = report/audit['portable_catalog']
    payload = json.loads(portable.read_text('utf-8'))
    assert payload['correction'] == json.loads(index.read_text('utf-8'))['correction']
    assert payload['native_source_catalog_sha256'] == digest(index)
    assert payload['path_base'] == 'index_file_parent' and payload['statistics_or_scope_changed'] is False
    assert result['metadata']['native_fit_index_sha256'] == digest(index)
    assert result['metadata']['original_native_fit_index_sha256'] == digest(index)
    for row, native in zip(payload['fits'],json.loads(index.read_text('utf-8'))['fits'],strict=True):
        assert not Path(row['path']).is_absolute() and '\\' not in row['path']
        copied=(portable.parent/row['path']).resolve()
        assert copied.is_relative_to(report) and digest(copied) == native['sha256']
        assert copied.read_bytes() == (collection/native['path']).read_bytes()
        assert row['scope'] == native['scope'] and row['task_id'] == native['task_id']
        assert row['native_source_collection_path'] == native['path']
    original_link = next(link for link in result['data_links'] if link.get('path')==str(index))
    assert 'provenance' in original_link['label']
    assert before == {str(p):(p.read_bytes(),p.stat().st_mtime_ns) for p in collection.rglob('*') if p.is_file()}


def test_twelve_comparisons_and_168_fits_have_complete_separate_indexes(tmp_path):
    report=tmp_path/'report';report.mkdir()
    collection=tmp_path/'collection'
    catalogs=[]
    for number in range(12):
        cid=f'comparison_{number:02}_B_vs_A'
        panel,root,_=_family(tmp_path,count=14,cid=cid,collection=collection)
        result,audit=prepare_portable_mapp_panel(panel,report_root=report,copy_file=_copy_into(report))
        catalogs.append(audit)
        portable=report/audit['portable_catalog']
        data=json.loads(portable.read_text('utf-8'))
        assert data['comparison_id']==cid and len(data['fits'])==14
        assert len({r['task_id'] for r in data['fits']})==14
        assert all(r['scope']['contrast_ids']==[cid] for r in data['fits'])
        assert result['metadata']['portable_fit_index_sha256']==digest(portable)
    assert len({a['portable_catalog'] for a in catalogs})==12
    assert sum(a['fit_count'] for a in catalogs)==168
    for audit in catalogs:
        path=report/audit['portable_catalog']
        for row in json.loads(path.read_text('utf-8'))['fits']:
            assert digest((path.parent/row['path']).resolve())==row['sha256']


@pytest.mark.parametrize('mutation', ['missing_fit','changed_fit','changed_catalog','catalog_comparison',
    'fit_comparison','duplicate_task','missing_task','foreign_task','parent_path','absolute_path','windows_drive',
    'backslash_path','invalid_sha','missing_root','nested_output'])
def test_invalid_catalog_or_fit_never_gets_accepted(tmp_path,mutation):
    panel,collection,index=_family(tmp_path)
    data=json.loads(index.read_text('utf-8'))
    if mutation=='missing_fit':
        (collection/data['fits'][0]['path']).unlink()
    elif mutation=='changed_fit':
        (collection/data['fits'][0]['path']).write_bytes(b'changed native statistics')
    elif mutation=='changed_catalog':
        index.write_text('{}',encoding='utf-8')
    elif mutation=='missing_root':
        panel['data_links'][0].pop('source_collection_root')
    elif mutation=='catalog_comparison':
        data['comparison_id']='OTHER'
    elif mutation=='fit_comparison':
        data['fits'][0]['scope']['contrast_ids']=['OTHER']
    elif mutation=='duplicate_task':
        data['fits'][1]['task_id']=data['fits'][0]['task_id']
    elif mutation=='missing_task':
        data['fits'].pop()
    elif mutation=='foreign_task':
        data['fits'][0]['task_id']='not-declared'
    elif mutation=='parent_path':
        data['fits'][0]['path']='../outside.tsv'
    elif mutation=='absolute_path':
        data['fits'][0]['path']='/outside.tsv'
    elif mutation=='windows_drive':
        data['fits'][0]['path']='C:/outside.tsv'
    elif mutation=='backslash_path':
        data['fits'][0]['path']='..\\outside.tsv'
    elif mutation=='invalid_sha':
        data['fits'][0]['sha256']='not-a-digest'
    if mutation not in {'missing_fit','changed_fit','changed_catalog','missing_root','nested_output'}:
        _write(index,data)
        panel['metadata']['native_fit_index_sha256']=digest(index)
    report=collection/'nested-report' if mutation=='nested_output' else tmp_path/'report'
    calls=[]
    def forbidden_copy(*args,**kwargs):
        calls.append(args)
        raise AssertionError('Invalid source must be rejected before copying fits')
    with pytest.raises((ValueError,FileNotFoundError)):
        _bind(panel,collection,report,copy=forbidden_copy)
    assert calls==[]


def test_wrong_copy_bytes_cannot_be_certified(tmp_path):
    panel,collection,_=_family(tmp_path)
    report=tmp_path/'report';report.mkdir()
    def wrong_copy(source,*,prefix):
        target=report/'wrong.tsv';target.write_bytes(b'wrong numerical bytes')
        return target.name
    with pytest.raises(ValueError,match='does not match'):
        prepare_portable_mapp_panel(panel,report_root=report,copy_file=wrong_copy)


def _bound_family(tmp_path):
    panel,collection,_=_family(tmp_path)
    report=tmp_path/'report';report.mkdir()
    copy=_copy_into(report)
    panel,audit=prepare_portable_mapp_panel(panel,report_root=report,copy_file=copy)
    for link in panel['data_links']:
        source=Path(link['path'])
        link['path']=source.relative_to(report).as_posix() if source.is_relative_to(report) else copy(source)
        link['source_sha256']=link['portable_sha256']=digest(report/link['path'])
    return panel,report,audit


def test_final_bound_catalog_can_be_checked_without_source_collection(tmp_path):
    panel,report,audit=_bound_family(tmp_path)
    shutil.rmtree(tmp_path/'collection')  # controlled fixture; report is now independent
    result=validate_portable_mapp_panel(panel,report_root=report)
    assert result['passed'] and result['fit_count']==3
    assert result['source_machine_paths_required'] is False


@pytest.mark.parametrize('mutation',['missing_index','changed_fit','missing_fit','wrong_comparison',
    'changed_scope','missing_row','duplicated_row','traversal','changed_correction','changed_original_catalog',
    'wrong_parallel_family'])
def test_final_acceptance_rejects_missing_or_changed_portable_evidence(tmp_path,mutation):
    panel,report,audit=_bound_family(tmp_path)
    portable=report/audit['portable_catalog']
    payload=json.loads(portable.read_text('utf-8'))
    if mutation=='missing_index':
        portable.unlink()
    elif mutation=='changed_fit':
        (portable.parent/payload['fits'][0]['path']).resolve().write_bytes(b'changed native numbers')
    elif mutation=='missing_fit':
        (portable.parent/payload['fits'][0]['path']).resolve().unlink()
    elif mutation=='changed_original_catalog':
        link=next(l for l in panel['data_links'] if l['source_sha256']==panel['metadata']['native_fit_index_sha256'])
        (report/link['path']).write_bytes(b'changed native source catalogue')
    elif mutation=='wrong_parallel_family':
        panel['metadata']['calibrated_motiftrace_parallel']=True
    else:
        if mutation=='wrong_comparison':payload['comparison_id']='OTHER'
        elif mutation=='changed_scope':payload['fits'][0]['scope']['contrast_ids']=['OTHER']
        elif mutation=='missing_row':payload['fits'].pop()
        elif mutation=='duplicated_row':payload['fits'][1]=deepcopy(payload['fits'][0])
        elif mutation=='traversal':payload['fits'][0]['path']='../../../../outside.tsv'
        elif mutation=='changed_correction':payload['correction']='Invented pooled FDR'
        _write(portable,payload)
        old=panel['metadata']['portable_fit_index_sha256']
        panel['metadata']['portable_fit_index_sha256']=digest(portable)
        for link in panel['data_links']:
            if link['source_sha256']==old:
                link['source_sha256']=link['portable_sha256']=digest(portable)
    with pytest.raises((ValueError,FileNotFoundError)):
        validate_portable_mapp_panel(panel,report_root=report)


def test_fresh_completed_view_builds_a_portable_report_without_previous_snapshot(tmp_path):
    from test_flow_mapp_views import inputs
    from katana_next.flow.mapp_views import build_mapp_views
    from katana_next.reports import build_report
    root,plan,receipts,transport,table=inputs(tmp_path)
    views=build_mapp_views(root,plan,receipts,transport)
    before={str(p):digest(p) for p in root.rglob('*') if p.is_file()}
    out=tmp_path/'fresh-report'
    manifest=build_report({'project_manifest':str(root/'project.json'),'artifacts':[],
                           'extra_views':views['extra_views'],'dummy':True,'report_require_plot_data':True},out)
    panel=next(p for p in manifest['panels'] if p.get('metadata',{}).get('component')=='mapp_native')
    assert validate_portable_mapp_panel(panel,report_root=out)['passed']
    link=next(link for link in panel['data_links'] if link.get('source_sha256')==panel['metadata']['portable_fit_index_sha256'])
    index=out/link['path']
    native=json.loads(index.read_text('utf-8'))
    for row in native['fits']:
        target=(index.parent/row['path']).resolve()
        assert target.is_relative_to(out) and digest(target)==row['sha256']
        assert target.read_bytes()==table.read_bytes()
    assert all(digest(Path(p))==sha for p,sha in before.items())
    with zipfile.ZipFile(out/'REPORT.zip') as archive:
        assert archive.testzip() is None
        for row in native['fits']:
            target=(index.parent/row['path']).resolve()
            assert archive.read(target.relative_to(out).as_posix())==table.read_bytes()


def test_moved_collection_relative_descriptor_needs_no_old_machine_path(tmp_path):
    from test_flow_mapp_views import inputs
    from katana_next.flow.mapp_views import build_mapp_views
    from katana_next.reports import build_report
    root,plan,receipts,transport,table=inputs(tmp_path)
    view=build_mapp_views(root,plan,receipts,transport)['extra_views'][0]
    for section in ('panels','tables'):
        for record in view[section]:
            record['path']=Path(record['path']).relative_to(root).as_posix()
            for link in record.get('data_links',[]):
                link['path']=Path(link['path']).relative_to(root).as_posix()
                if 'source_collection_root' in link:
                    link['source_collection_root']='.'
    descriptor=_write(root/'report-project.json',{'project_manifest':'project.json','artifacts':[],
                       'extra_views':[view],'dummy':True,'report_require_plot_data':True})
    moved=tmp_path/'moved-to-other-machine';root.rename(moved)
    out=tmp_path/'report'
    manifest=build_report(moved/'report-project.json',out)
    panel=next(p for p in manifest['panels'] if p.get('metadata',{}).get('component')=='mapp_native')
    index=out/next(l['path'] for l in panel['data_links'] if l.get('source_sha256')==panel['metadata']['portable_fit_index_sha256'])
    assert all((index.parent/r['path']).resolve().is_relative_to(out) and digest((index.parent/r['path']).resolve())==r['sha256']
               for r in json.loads(index.read_text('utf-8'))['fits'])
    assert validate_portable_mapp_panel(panel,report_root=out)['passed']


def test_compacted_fresh_mapp_report_preserves_native_catalog_and_usable_fit_index(tmp_path):
    from test_flow_mapp_views import inputs
    from katana_next.flow.mapp_views import build_mapp_views
    from katana_next.reports import build_report
    from katana_next.reports.compact import compact_report
    root,plan,receipts,transport,table=inputs(tmp_path)
    views=build_mapp_views(root,plan,receipts,transport)
    out=tmp_path/'report'
    manifest=build_report({'project_manifest':str(root/'project.json'),'artifacts':[],
                          'extra_views':views['extra_views'],'dummy':True,'report_require_plot_data':True},out)
    panel=next(p for p in manifest['panels'] if p.get('metadata',{}).get('component')=='mapp_native')
    before_native=panel['metadata']['native_fit_index_sha256']
    (out/'REPORT.zip').unlink()  # owned generated fixture, before its separate storage pass
    alias=out/'outputs/results/a.tsv'
    shutil.copyfile(table,alias)
    compact_report(out)
    final=json.loads((out/'report-manifest.json').read_text('utf-8'))
    panel=next(p for p in final['panels'] if p.get('metadata',{}).get('component')=='mapp_native')
    assert panel['metadata']['native_fit_index_sha256']==before_native
    assert validate_portable_mapp_panel(panel,report_root=out)['passed']
