"""Fail closed before a preserved-plan terminal handoff can submit anything."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


MODULE = Path(__file__).resolve().parents[2] / 'pipeline/tools/check_completed_handoff.py'
SPEC = importlib.util.spec_from_file_location('completed_handoff_gate', MODULE)
handoff = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(handoff)


def save(path, value):
    path.write_text(json.dumps(value), encoding='utf-8')


def sign(receipt):
    body = {key: value for key, value in receipt.items() if key != 'receipt_sha256'}
    receipt['receipt_sha256'] = hashlib.sha256(json.dumps(
        body, sort_keys=True, separators=(',', ':'), ensure_ascii=False,
        default=str).encode()).hexdigest()
    return receipt


def fixture(tmp_path):
    model = tmp_path / 'models.tsv'
    model.write_text('motif\tsequence\nM1\tACGU\n', encoding='utf-8')
    product = tmp_path / 'completed.tsv'
    product.write_text('value\n1\n', encoding='utf-8')
    plan = {'plan_id': 'unchanged-original-plan', 'tasks': [{'id': 'first'}, {'id': 'second'}]}
    save(tmp_path / 'plan.json', plan)
    results = []
    for identity in ('first', 'second'):
        directory = tmp_path / identity
        directory.mkdir()
        receipt = sign({'task_id': identity, 'plan_id': plan['plan_id'],
                        'status': 'completed', 'inputs': [{'path': str(model)}],
                        'outputs': [{'path': str(product)}], 'producer_markers': []})
        save(directory / 'task.receipt.json', receipt)
        results.append(identity)
    inputs = tmp_path / 'inputs.json'
    save(inputs, {'plan': 'plan.json', 'results': results, 'models': ['models.tsv']})
    return inputs


def test_original_plan_coverage_and_relative_resource_identity(tmp_path):
    inputs = fixture(tmp_path)
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    data, plan, protected = handoff.bounded_inputs(inputs)
    assert len(plan['tasks']) == 2
    assert data['models'] == [str((tmp_path / 'models.tsv').resolve())]
    assert (tmp_path / 'completed.tsv').resolve() in protected
    assert handoff.inventory(protected)
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('change', ['tamper', 'missing', 'wrong_plan', 'failed', 'extra', 'duplicate'])
def test_rejects_unproven_task_before_output_or_submission(tmp_path, change):
    inputs = fixture(tmp_path)
    data = handoff.read(inputs)
    receipt_path = tmp_path / 'first/task.receipt.json'
    receipt = handoff.read(receipt_path)
    if change == 'tamper':
        receipt['outputs'].append({'path': 'unproven.tsv'})
    elif change == 'wrong_plan':
        receipt['plan_id'] = 'different-plan'
        sign(receipt)
    elif change == 'failed':
        receipt['status'] = 'failed'
        sign(receipt)
    elif change == 'extra':
        receipt['task_id'] = 'not-in-original-plan'
        sign(receipt)
    elif change == 'missing':
        data['results'].pop()
    else:
        plan = handoff.read(tmp_path / 'plan.json')
        plan['tasks'].append({'id': 'first'})
        save(tmp_path / 'plan.json', plan)
    save(receipt_path, receipt)
    save(inputs, data)
    output = tmp_path / 'attempt'
    with pytest.raises(ValueError):
        handoff.main(['--inputs', str(inputs), '--outdir', str(output), '--image', 'unused'])
    assert not output.exists()


def test_trace_rejects_any_scientific_or_extra_process(tmp_path):
    trace = tmp_path / 'trace.tsv'
    trace.write_text('name\tstatus\texit\nKATANA_FLOW_COLLECT\tCOMPLETED\t0\n'
                     'KATANA_FLOW_REPORT\tCOMPLETED\t0\n'
                     'KATANA_FLOW_NETWORK\tCOMPLETED\t0\n', encoding='utf-8')
    with pytest.raises(ValueError, match='exactly the two'):
        handoff.trace(trace)


def test_preserves_existing_attempt_even_when_inputs_are_valid(tmp_path):
    inputs = fixture(tmp_path)
    output = tmp_path / 'attempt'
    output.mkdir()
    marker = output / 'old-failure.json'
    marker.write_bytes(b'preserve this lifecycle')
    with pytest.raises(FileExistsError):
        handoff.main(['--inputs', str(inputs), '--outdir', str(output), '--image', 'unused'])
    assert marker.read_bytes() == b'preserve this lifecycle'


def test_readonly_source_cannot_contain_new_output(tmp_path):
    inputs = fixture(tmp_path)
    output = tmp_path / 'attempt'
    with pytest.raises(ValueError, match='independent'):
        handoff.main(['--inputs', str(inputs), '--outdir', str(output), '--image', 'unused',
                      '--bind-readonly', str(tmp_path)])
    assert not output.exists()


def test_inventory_detects_content_change_with_same_size_and_mtime(tmp_path):
    import os
    product = tmp_path / 'product.tsv'
    product.write_bytes(b'original')
    before = handoff.inventory({product})
    stamp = product.stat()
    product.write_bytes(b'changed!')
    os.utime(product, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
    assert handoff.inventory({product}) != before
