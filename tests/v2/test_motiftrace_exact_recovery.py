from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import resume_motiftrace_exact5712_20261016 as recovery


def model(model_id='p_model'):
    return SimpleNamespace(model_id=model_id, kind='PWM', width=4)


def tables():
    layout = pd.DataFrame({'within_method_row_id': [0, 1],
                           'status': ['estimable', 'insufficient_independent_pairs']})
    numbers = pd.DataFrame({'within_method_row_id': [0], 'method': ['PWM_mean_odds'],
                            'P_exact': [.5], 'paired_rank_Z': [0.], 'nonzero_gene_pairs': [1],
                            'minimum_attainable_P': [.5], 'mean_paired_score_difference': [.1]})
    return layout, numbers, pd.DataFrame({'event_id': ['event'], 'score': [.1]})


def completed(tmp_path, model_id='p_model'):
    (tmp_path / 'models').mkdir(exist_ok=True)
    layout, numbers, forms = tables()
    recovery.commit_model(tmp_path, model(model_id), 1, 'PBS', numbers, forms, .01, False, 'seal')
    return layout, numbers, forms


def test_completed_model_hashes_reload_and_missing_model_is_not_completed(tmp_path):
    layout, _, _ = completed(tmp_path)
    result = recovery.validate_model(tmp_path, model(), 1, 'PBS', layout)
    assert result['all_estimable_hypotheses'] == 1
    assert recovery.validate_model(tmp_path, model('missing'), 2, 'PBS', layout) is None


def test_completed_model_is_never_overwritten(tmp_path):
    _, numbers, forms = completed(tmp_path)
    before = recovery.artifact(tmp_path / 'models/p_model.COMPLETE.json')
    with pytest.raises(FileExistsError):
        recovery.commit_model(tmp_path, model(), 1, 'PBS', numbers, forms, 1, False, 'different')
    assert recovery.artifact(tmp_path / 'models/p_model.COMPLETE.json') == before


@pytest.mark.parametrize('field,value', [('comparison_id', 'KA'), ('model_id', '../escape'),
                                         ('ordinal', 2), ('method', 'different test')])
def test_changed_receipt_identity_refused(tmp_path, field, value):
    layout, _, _ = completed(tmp_path)
    path = tmp_path / 'models/p_model.COMPLETE.json'
    evidence = recovery.read(path)
    evidence[field] = value
    recovery.save(path, evidence)
    with pytest.raises(ValueError):
        recovery.validate_model(tmp_path, model(), 1, 'PBS', layout)


def test_corrupt_completed_file_refused_not_recomputed(tmp_path):
    layout, _, _ = completed(tmp_path)
    with (tmp_path / 'models/p_model.numbers.parquet').open('ab') as stream:
        stream.write(b'corruption')
    with pytest.raises(ValueError, match='hash/size'):
        recovery.validate_model(tmp_path, model(), 1, 'PBS', layout)


def test_completed_hypothesis_count_must_match_full_frozen_geometry(tmp_path):
    layout, _, _ = completed(tmp_path)
    layout['status'] = 'estimable'
    with pytest.raises(ValueError, match='hypothesis count'):
        recovery.validate_model(tmp_path, model(), 1, 'PBS', layout)


def test_uncommitted_partial_is_preserved_before_retry(tmp_path):
    folder, attempt = tmp_path / 'PBS', tmp_path / 'attempt'
    (folder / 'models').mkdir(parents=True)
    path = folder / 'models/p_model.numbers.parquet'
    path.write_bytes(b'interrupted-data')
    recovery.archive_uncommitted(folder, 'p_model', attempt)
    assert not path.exists()
    assert (attempt / 'uncommitted_preserved/PBS/p_model.numbers.parquet').read_bytes() == b'interrupted-data'


def test_interrupted_commit_does_not_make_completed_marker(tmp_path, monkeypatch):
    (tmp_path / 'models').mkdir()
    _, numbers, forms = tables()
    original = recovery.atomic_parquet
    calls = []

    def fail_second(frame, path):
        calls.append(path)
        if len(calls) == 2:
            raise OSError('simulated power loss before receipt')
        return original(frame, path)

    monkeypatch.setattr(recovery, 'atomic_parquet', fail_second)
    with pytest.raises(OSError):
        recovery.commit_model(tmp_path, model(), 1, 'PBS', numbers, forms, 1, False, 'seal')
    assert (tmp_path / 'models/p_model.numbers.parquet').exists()
    assert not (tmp_path / 'models/p_model.COMPLETE.json').exists()


def test_reused_forms_preserve_original_bytes(tmp_path):
    (tmp_path / 'models').mkdir()
    _, numbers, forms = tables()
    source = tmp_path / 'old_forms.parquet'
    forms.to_parquet(source, index=False)
    recovery.commit_model(tmp_path, model(), 1, 'PBS', numbers, source, 1, True, 'seal')
    assert recovery.artifact(source) == recovery.artifact(tmp_path / 'models/p_model.all_DAS_forms.parquet')


def test_kernel_lock_blocks_duplicates_and_releases_after_exit(tmp_path):
    path = tmp_path / 'lock'
    with recovery.exclusive_lock(path):
        with pytest.raises(BlockingIOError):
            with recovery.exclusive_lock(path):
                pytest.fail('second worker acquired existing lock')
    with recovery.exclusive_lock(path):
        pass


def test_publish_replays_without_overwrite_and_preserves_different_partial(tmp_path, monkeypatch):
    monkeypatch.setattr(recovery, 'ROOT', tmp_path / 'base/upstream')
    source = tmp_path / 'source'
    source.write_bytes(b'completed')
    target = tmp_path / 'base/upstream/result'
    target.parent.mkdir(parents=True)
    target.write_bytes(b'partial')
    attempt = tmp_path / 'attempt'
    recovery.publish_file(source, target, attempt)
    assert target.read_bytes() == b'completed'
    assert (attempt / 'replaced_incomplete_publications/upstream/result').read_bytes() == b'partial'
    before = target.stat().st_mtime_ns
    recovery.publish_file(source, target, attempt)
    assert target.stat().st_mtime_ns == before


def test_resume_skips_completed_compute_and_continues_next_model(tmp_path, monkeypatch):
    root, prior = tmp_path / 'root', tmp_path / 'prior'
    folder = root / 'PBS'
    (folder / 'models').mkdir(parents=True)
    (folder / 'layouts').mkdir()
    (prior / 'PBS').mkdir(parents=True)
    recovery.save(root / 'recovery_v3/REGISTRATION.json', {'test': True})
    recovery.save(root / 'FAMILY_FREEZE.json', {'unchanged': True})
    monkeypatch.setattr(recovery, 'ROOT', root)
    monkeypatch.setattr(recovery, 'RECOVERY', root / 'recovery_v3')
    layout, numbers, forms = tables()
    layout.to_parquet(folder / 'layouts/width_4.parquet', index=False)
    recovery.commit_model(folder, model('p_done'), 1, 'PBS', numbers, forms, .01, False, 'old')
    before = {p.name: recovery.artifact(p) for p in (folder / 'models').iterdir()}
    inventory = {'PBS': {'p_done': recovery.validate_model(folder, model('p_done'), 1, 'PBS', layout)}}
    pd.DataFrame({'event_id': ['event']}).to_parquet(prior / 'PBS/ALL_ORIGINAL_DAS.parquet', index=False)
    models = pd.DataFrame([{'model_id': value, 'width': 4, 'kind': 'PWM'} for value in ['p_done', 'p_new']])
    rna = pd.DataFrame({'sequence': ['AAAA'], 'role': ['foreground'], 'layer': ['mature_local_form']})
    calls = []

    class Packed:
        def __init__(self, sequences):
            pass

        def activate(self, matrix):
            pass

    def compute(plan, current, packed, engine, geometry):
        calls.append(current['model_id'])
        return numbers.copy()

    def finish(frozen, case, catalog, total, attempt):
        recovery.save(root / case / 'COMPLETE.json', {'verified_raw_rows': total})

    frozen = SimpleNamespace(CASES=['PBS'], PRIOR=prior, catalog=lambda: (models, {}, {'p_done': [], 'p_new': []}),
                             ExactSignedRank=lambda: object(), data_for=lambda case: (rna, None, object()),
                             PackedPWM=Packed, compute_model=compute, safety=lambda value: None,
                             form_scores=lambda *args: forms.copy())
    monkeypatch.setattr(recovery, 'finish_case_resumable', finish)
    recovery.analyze_resumable(frozen, {'families': {'PBS': {'all_estimable_hypotheses': 2}}}, inventory, tmp_path / 'attempt')
    assert calls == ['p_new']
    for name, expected in before.items():
        assert recovery.artifact(folder / 'models' / name) == expected
    assert recovery.validate_model(folder, model('p_new'), 2, 'PBS', layout)
    # A second restart with both committed receipts performs zero model computations.
    calls.clear()
    recovery.analyze_resumable(frozen, {'families': {'PBS': {'all_estimable_hypotheses': 2}}}, inventory, tmp_path / 'attempt2')
    assert calls == []
