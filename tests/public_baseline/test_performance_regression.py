"""Generated-input regressions against the pre-optimization implementation.

References are byte-for-byte sources from commit 04941d6, not study outputs.
No reference source or accepted scientific parameter is rewritten by this test.
"""
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / 'tests/fixtures/performance_reference_20261026'


def reference(path):
    record = json.loads((REFERENCE / 'sources.json').read_text())
    source = REFERENCE / Path(path).name
    assert hashlib.sha256(source.read_bytes()).hexdigest() == record['sources'][path]
    name = path[:-3].replace('/', '.')
    name = name.rsplit('.', 1)[0] + '._performance_reference_' + source.stem
    spec = importlib.util.spec_from_file_location(name, source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


MODELS = {'literal': {'kind': 'consensus', 'pattern': 'ACG', 'RBPs': 'fixture'},
          'pwm': {'kind': 'PWM', 'matrix': [[.7,.1,.1,.1],[.1,.6,.2,.1],[.1,.1,.1,.7]]}}


@pytest.mark.parametrize('sequence', ['ACGT'*17, 'ACGNNTACG'*8, '', 'AC', 'acguACGTN'])
@pytest.mark.parametrize('width', [0, 10, 20, 50])
def test_cached_profiles_equal_reference_at_every_position(sequence, width):
    from katana_next.science import motiftrace_clusters as current
    from katana_next.science.motiftrace_events import ScoringCache
    old = reference('katana_next/science/motiftrace_clusters.py')
    cache = ScoringCache()
    for model in MODELS.values():
        methods = ['PWM_mean_odds'] if model['kind'] == 'PWM' else ['finite_nonoverlap_counts', 'overlapping_union_coverage']
        for method in methods:
            np.testing.assert_array_equal(current._profile(sequence, model, method, width, scorer=cache),
                                          old._profile(sequence, model, method, width))


def test_scans_are_reused_and_eviction_changes_no_scores(monkeypatch):
    from katana_next.science import motiftrace_events as events, motiftrace_clusters as clusters
    original = events.score_positions
    calls = []
    def counted(*args):
        calls.append(args)
        return original(*args)
    monkeypatch.setattr(events, 'score_positions', counted)
    cache = events.ScoringCache()
    for width in [0,10,20,50]:
        for method in ['finite_nonoverlap_counts','overlapping_union_coverage']:
            clusters._profile('ACGT'*20, MODELS['literal'], method, width, scorer=cache)
    assert len(calls) == 1  # previously eight complete sequence scans
    uncached = events.ScoringCache(max_bytes=1)
    for model in MODELS.values():
        observed = uncached('ACGT'*20, model)
        expected = original('ACGT'*20, model)
        for a,b in zip(observed[:2], expected[:2]):
            np.testing.assert_array_equal(a,b)
    assert not uncached.values


def paired_fixture():
    rows, matches = [], []
    rng = np.random.default_rng(33)
    for index in range(22):
        fg, bg = f'f{index}', f'b{index}'
        matches.append(dict(status='matched',foreground_gene=fg,background_gene=bg,
                            foreground_event=fg,background_event=bg,event_type='SE',direction='inclusion_up'))
        for event, seq in [(fg, 'ACG'*22), (bg, ''.join(rng.choice(list('ACGT'),66)))]:
            for region in ['body', 'upstream_acceptor_intron']:
                rows.append(dict(event_id=event,event_type='SE',layer='affected_sequence',region=region,sequence=seq))
    return pd.DataFrame(rows), pd.DataFrame(matches)


def test_streamed_family_preserves_all_positions_null_P_BH_and_FWER():
    from katana_next.science import motiftrace_clusters as current
    old = reference('katana_next/science/motiftrace_clusters.py')
    rna, matching = paired_fixture()
    expected, ids = old.build_paired_datasets(rna,matching,MODELS,widths=(10,20))
    stream, actual_ids = current.build_paired_datasets(rna,matching,MODELS,widths=(10,20),stream=True)
    assert iter(stream) is stream and actual_ids == ids
    for observed, baseline in zip(stream, expected, strict=True):
        assert {k:v for k,v in observed.items() if k!='differences'} == {k:v for k,v in baseline.items() if k!='differences'}
        np.testing.assert_array_equal(observed['differences'], baseline['differences'])
    stream, _ = current.build_paired_datasets(rna,matching,MODELS,widths=(10,20),stream=True)
    actual = current.paired_cluster_family(stream,ids,permutations=99,seed=71)
    baseline = old.paired_cluster_family(expected,ids,permutations=99,seed=71)
    for observed, wanted in zip(actual[:2],baseline[:2]):
        pd.testing.assert_frame_equal(observed,wanted,check_exact=True)
    np.testing.assert_array_equal(actual[2],baseline[2])
    assert len(actual[1]) > 0  # exercise actual family-wise corrected clusters


def test_native_network_precedes_python_dense_allocations_without_result_changes(monkeypatch):
    from katana_cosplice import network, r_wgcna
    from katana_cosplice.config import NetworkConfig
    old = reference('katana_cosplice/network.py')
    rng = np.random.default_rng(17)
    matrix = pd.DataFrame(rng.uniform(size=(42,16)),index=[f'e{i}' for i in range(42)])
    qc = pd.DataFrame({'stability_score':1.},index=matrix.index)
    settings = replace(NetworkConfig(),module_detection='wgcna_blockwise',min_module_size=20)
    order = []
    def labels(selected, **kwargs):
        order.append('native')
        return np.array(['M01']*21+['M02']*21)
    monkeypatch.setattr(r_wgcna,'exact_wgcna_labels',labels)
    original = network._tom_hierarchy
    def dense(*args,**kwargs):
        order.append('python_dense')
        return original(*args,**kwargs)
    monkeypatch.setattr(network,'_tom_hierarchy',dense)
    observed = network.fit_network(matrix,qc,settings)
    assert order[:2] == ['native','python_dense']
    expected = old.fit_network(matrix,qc,settings)
    for field, value in vars(expected).items():
        actual = getattr(observed,field)
        if isinstance(value,pd.DataFrame): pd.testing.assert_frame_equal(actual,value,check_exact=True)
        elif isinstance(value,pd.Series): pd.testing.assert_series_equal(actual,value,check_exact=True)
        elif isinstance(value,np.ndarray): np.testing.assert_array_equal(actual,value)
        else: assert actual == value


def test_plan_uses_independent_core_branches_and_ancestry_only_cache(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.flow.planner import create_plan
    from katana_next.flow.model import receipt_plan_id, seal_plan, validate_plan
    source = create_full_fixture(tmp_path/'study')
    plan = create_plan(source,{})
    molecular = [t for t in plan['tasks'] if t['execution_group'] in {'as','deg','dtu'}]
    assert len(molecular) == 3*len(plan['contrast_ids'])
    assert all(len(t['scope']['contrast_ids'])==1 for t in molecular)
    alternative = create_plan(source,{'analyses':{'motiftrace':False}})
    for task in molecular:
        assert receipt_plan_id(plan,task['id']) == receipt_plan_id(alternative,task['id'])
    changed = json.loads(json.dumps(plan))
    target = molecular[0]
    next(t for t in changed['tasks'] if t['id']==target['id'])['parameters']['read_support'] = 77
    changed = seal_plan(changed)
    validate_plan(changed)
    assert receipt_plan_id(changed,target['id']) != receipt_plan_id(plan,target['id'])
    assert receipt_plan_id(changed,'audit') == receipt_plan_id(plan,'audit')
    other = next(t for t in molecular if t['scope']!=target['scope'])
    assert receipt_plan_id(changed,other['id']) == receipt_plan_id(plan,other['id'])


@pytest.mark.parametrize('eligible',[False,True])
def test_mapp_shared_preparation_equals_original_for_all_fourteen_windows(tmp_path,monkeypatch,eligible):
    from katana_next.flow import mapp_inputs as current, stages
    from test_flow_mapp_applicability import fixture
    old = reference('katana_next/flow/mapp_inputs.py')
    spec = fixture(tmp_path,mapped=True)
    # Test both measurable and TPM-ineligible events; no native scientific fit.
    if not eligible:
        for path in spec['inputs']['quantification'].values():
            table=pd.read_csv(path,sep='\t');table['TPM']=0.001;table.to_csv(path,sep='\t',index=False)
    old_shared=old.prepare_shared(spec['inputs'],spec['parameters'])
    shared=current.prepare_shared(spec['inputs'],spec['parameters'])
    prepared = []
    original=current.prepare_event_quantities
    def tracked(*args):
        value=original(*args);prepared.append(id(value));return value
    monkeypatch.setattr(current,'prepare_event_quantities',tracked)
    def counter_fixture(command,out,label):
        assert eligible
        # Only downstream preparation is under test; native kmer math is unchanged.
        pd.DataFrame(np.zeros((1,5440)),index=['SE:1']).to_csv(out/'NATIVE_SITECOUNTS.tsv',sep='\t')
        (out/(label+'.command.json')).write_text('{}')
    for index,window in enumerate(current.native_windows()):
        task={**spec,'parameters':{**spec['parameters'],**window}}
        before=tmp_path/f'old-{index}';after=tmp_path/f'new-{index}';before.mkdir();after.mkdir()
        old.prepare_inputs(task,before,run=counter_fixture,sha=stages._sha,executable=lambda *args:args[0],shared=old_shared)
        current.prepare_inputs(task,after,run=counter_fixture,sha=stages._sha,executable=lambda *args:args[0],shared=shared)
        for path in before.glob('*.tsv'):
            pd.testing.assert_frame_equal(pd.read_csv(after/path.name,sep='\t'),pd.read_csv(path,sep='\t'),check_exact=True)
        if eligible:
            assert (after/'sequences.fasta').read_bytes()==(before/'sequences.fasta').read_bytes()
    assert len(prepared)==14 and len(set(prepared))==1


def test_owned_payload_is_retained_without_copy_and_borrowed_input_is_isolated(tmp_path):
    import os
    from katana_next.reports.storage import copy_payload
    owned=tmp_path/'collection';owned.mkdir()
    source=owned/'immutable.tsv';source.write_bytes(b'gene\tvalue\nG1\t2\n')
    destination=tmp_path/'report/payload.tsv'
    copy_payload(source,destination,owned_root=owned)
    assert os.path.samefile(source,destination)
    source.unlink()  # removing an owned work pathname cannot break delivery
    assert destination.read_bytes()==b'gene\tvalue\nG1\t2\n'
    borrowed=tmp_path/'user-input.tsv';borrowed.write_bytes(destination.read_bytes())
    other=tmp_path/'report/borrowed.tsv'
    copy_payload(borrowed,other,owned_root=owned)
    assert not os.path.samefile(borrowed,other)
    borrowed.write_bytes(b'changed')
    assert other.read_bytes()==destination.read_bytes()


def test_archive_keeps_file_bytes_without_recompressing_compressed_payloads(tmp_path):
    import gzip
    import zipfile
    from katana_next.reports.refresh import write_archive
    root=tmp_path/'report';root.mkdir()
    content={'REPORT.html':b'<html>fixed report</html>',
             'table.tsv.gz':gzip.compress(b'gene\tvalue\nG1\t2\n'),
             'image.png':b'generated compression fixture'}
    for name,data in content.items(): (root/name).write_bytes(data)
    archive=write_archive(root)
    with zipfile.ZipFile(archive) as reader:
        assert set(reader.namelist())==set(content)
        for name,data in content.items():
            assert reader.read(name)==data
            assert reader.getinfo(name).compress_type==(zipfile.ZIP_DEFLATED if name.endswith('.html') else zipfile.ZIP_STORED)


def test_comparison_branches_preserve_aggregated_reference_tables(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.flow.planner import create_plan
    from katana_next.flow.selection import DEFAULTS
    from katana_next.flow.model import write_json, read_json
    from katana_next.flow.runtime import execute_task
    old=reference('katana_next/flow/planner.py')
    source=create_full_fixture(tmp_path/'study')
    config={'design_factors':['condition'],'analyses':dict.fromkeys(DEFAULTS,False)}
    tables=[]
    for label,planner in [('old',old.create_plan),('new',create_plan)]:
        plan=planner(source,config)
        location=write_json(tmp_path/(label+'-plan.json'),plan)
        receipts={}; result={}
        for task in sorted(plan['tasks'],key=lambda t:t['level']):
            receipt=execute_task(location,task['id'],tmp_path/label/task['id'],
                                 dependency_receipts=[receipts[p] for p in task['dependencies']])
            receipts[task['id']]=receipt
            for artifact in read_json(receipt)['generated_artifacts']:
                if artifact['kind']=='table' and artifact['analysis_id'].startswith(('quality.','differential.','prioritize.evidence')):
                    result.setdefault((artifact['analysis_id'],artifact['name']),[]).append(pd.read_parquet(artifact['path']))
        tables.append(result)
    assert tables[0].keys()==tables[1].keys()
    for key in tables[0]:
        old_table=pd.concat(tables[0][key],ignore_index=True)
        new_table=pd.concat(tables[1][key],ignore_index=True)
        assert len(old_table)==len(new_table),key
        if 'row_id' in old_table:
            old_table=old_table.sort_values('row_id').reset_index(drop=True)
            new_table=new_table.sort_values('row_id').reset_index(drop=True)
        pd.testing.assert_frame_equal(new_table,old_table,check_exact=True,obj=str(key))


def test_public_execution_adapter_streams_and_matches_preoptimization_complete_tables(tmp_path,monkeypatch):
    from test_motiftrace_execution_budget import _generated
    from katana_next.science import motiftrace_execution as adapter
    old=reference('katana_next/science/motiftrace_consequences.py')
    old_clusters=reference('katana_next/science/motiftrace_clusters.py')
    old_events=reference('katana_next/science/motiftrace_events.py')
    old.build_paired_datasets=old_clusters.build_paired_datasets
    old.paired_cluster_family=old_clusters.paired_cluster_family
    old.trace_instances=old_events.trace_instances
    old.paired_form_scores=old_events.paired_form_scores
    old.render_consequences=lambda *args,**kwargs: None
    monkeypatch.setattr(adapter,'render_consequences',lambda *args,**kwargs:None)
    build=adapter.build_paired_datasets
    observed=[]
    def streaming(*args,**kwargs):
        observed.append(kwargs.get('stream'))
        datasets,ids=build(*args,**kwargs)
        assert iter(datasets) is datasets
        return datasets,ids
    monkeypatch.setattr(adapter,'build_paired_datasets',streaming)
    frame,ids,models,fasta=_generated(tmp_path)
    old.analyze_consequences(frame,ids,models,fasta,tmp_path/'old',permutations=99)
    adapter.analyze_consequences(frame,ids,models,fasta,tmp_path/'new',permutations=99)
    assert observed==[True]
    tables=list((tmp_path/'old').glob('*.tsv*'))
    assert len(tables)>=8
    for path in tables:
        pd.testing.assert_frame_equal(pd.read_csv(tmp_path/'new'/path.name,sep='\t'),pd.read_csv(path,sep='\t'),check_exact=True,obj=path.name)
    np.testing.assert_array_equal(np.load(tmp_path/'new/MAX_CLUSTER_NULL.npy'),np.load(tmp_path/'old/MAX_CLUSTER_NULL.npy'))
    assert (tmp_path/'new/ACTUAL_LOCAL_RNA.fa').read_bytes()==(tmp_path/'old/ACTUAL_LOCAL_RNA.fa').read_bytes()
