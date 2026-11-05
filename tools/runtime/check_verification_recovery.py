#!/usr/bin/env python3
"""Installed image regression: generated DTU bytes and native-table report gate.

The native MAPP lifecycle is a synthetic receipt fixture, not a new model fit.
The renderer, portable data transport and final acceptance gate are real.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
from unittest.mock import patch


def verify(root):
    import pandas as pd
    from katana_next.flow import file_io
    from katana_next.flow.runtime import file_record, canonical_hash
    from katana_next.flow.model import read_json, write_json, file_sha256, seal_plan
    from katana_next.flow.mapp_views import build_mapp_views
    from katana_next.flow.acceptance import required_report_gaps
    from katana_next.full_fixture import create_full_fixture
    from katana_next.reports import build_report

    data = root / 'dtu_features.tsv'; data.write_text('transcript_id\tfdr\nT1\t0.02\n')
    read = file_io.hashlib.file_digest
    calls = []
    def transition(handle, algorithm):
        result = read(handle, algorithm); calls.append(1)
        if len(calls) == 1:
            stamp = data.stat()
            os.utime(data, ns=(stamp.st_atime_ns, stamp.st_mtime_ns + 1000000000))
        return result
    with patch.dict(os.environ, KATANA_VERIFY_CONTENT='1'), patch.object(file_io.hashlib, 'file_digest', transition):
        record = file_record(data)
    assert record['sha256'] == file_sha256(data)
    assert record['file_identity'] == file_io.identity(data)

    source = create_full_fixture(root / 'study')
    collection = root / 'collection'; collection.mkdir()
    shutil.copyfile(source, collection / 'project.json')
    project = read_json(source); cid = project['contrasts'][0]['contrast_id']
    task = dict(id='native', workflow_id='flow.mapp.maei', adapter='mapp.maei', dependencies=[],
                parameters={}, scope=dict(contrast_ids=[cid], site='3ss', window_start=-50, window_end=0))
    plan = seal_plan(dict(project_manifest=str(source), project_manifest_sha256=file_sha256(source),
        contrast_ids=[cid], configuration={'design_factors':['condition']}, tasks=[task]))
    task = plan['tasks'][0]
    write_json(collection / 'plan.json', plan)
    original = root / 'native/AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv'; original.parent.mkdir()
    pd.DataFrame({'activities_A1':[-.2,.4], 'activities_B1':[.3,-.1], 'combined.Zscore':[2.,3.],
        'corr_pval_A1':[.1,.02], 'corr_pval_B1':[.04,.2]}, index=['window|AAA','window|CCC']).to_csv(original, sep='\t')
    native_hash = file_sha256(original)
    receipt = dict(schema='katana.flow.task_receipt/v1', task_id='native', plan_id=plan['plan_id'],
        task_contract_sha256=canonical_hash(task), source_project_manifest_sha256=plan['project_manifest_sha256'],
        status='completed', execution_evidence={'native_fitting_performed':True}, outputs=[file_record(original)],
        synthetic_lifecycle_fixture=True)
    receipt['receipt_sha256'] = canonical_hash(receipt)
    receipt_path = write_json(collection / 'evidence/native/task.receipt.json', receipt)
    target = collection / 'products/native.tsv'; target.parent.mkdir(); shutil.copyfile(original, target)
    transport = [{'source':file_record(original), 'collection_path':'products/native.tsv'}]
    write_json(collection / 'transport-provenance.json', {'products':transport})
    proof = dict(project_manifest='project.json', task_receipts=[dict(task_id='native',
        collection_path='evidence/native/task.receipt.json', sha256=file_sha256(receipt_path))])
    views = build_mapp_views(collection, plan, {'native':(receipt_path, receipt)}, transport)
    # A completed collection from the previous runtime has bare sample labels.
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    for view in views['extra_views']:
        for panel in view['panels']:
            table = pd.read_csv(next(x['path'] for x in panel['data_links'] if x['role'] == 'plot_data'), sep='\t')
            columns = [c for c in table if c.startswith('activities_')]
            fig = Figure(); FigureCanvasAgg(fig); axes = fig.subplots()
            axes.imshow(table[columns].to_numpy(), aspect='auto')
            axes.set_xticks(range(len(columns)), [c.removeprefix('activities_') for c in columns])
            fig.savefig(panel['path'])
            panel['metadata']['native_activity_image_sha256'] = file_sha256(panel['path'])
    out = root / 'report'
    report = build_report(dict(project_manifest=str(collection / 'project.json'), artifacts=[],
        extra_views=views['extra_views'], report_require_plot_data=True), out)
    def gaps():
        return [g for g in required_report_gaps(collection, plan, report, out, collection_receipt=proof)
                if g['component'] in {'mapp_native','mapp_parallel'}]
    assert not gaps(), gaps()
    panel = next(p for p in report['panels'] if p.get('metadata',{}).get('component') == 'mapp_native')
    assert panel['metadata']['activity_rendering'] == 'condition_sample_labels/v1'
    assert panel['metadata']['rendered_activity_image_sha256'] != panel['metadata']['native_activity_image_sha256']
    image = out / panel['image']; saved = image.read_bytes(); image.write_bytes(b'changed')
    assert gaps(), 'Modified activity image incorrectly accepted'
    image.write_bytes(saved)
    assert not gaps()
    assert file_sha256(original) == native_hash
    assert (out / 'REPORT.zip').stat().st_size > 0
    write_json(root / 'acceptance.json', dict(status='passed', generated_data_only=True,
        metadata_transition_verified=True, actual_report_gate=True, changed_image_rejected=True,
        native_fit_rerun=False, native_values_unchanged=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args(); args.outdir.mkdir(parents=True, exist_ok=False)
    verify(args.outdir)
