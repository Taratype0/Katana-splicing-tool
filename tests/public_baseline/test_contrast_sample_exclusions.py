"""An incomplete pair is excluded only from its explicitly selected contrast."""
import csv
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


validator = load('subset_validator', 'pipeline/bin/katana_validate_inputs.py')
suppa = load('subset_suppa', 'pipeline/bin/katana_subset_suppa.py')


def samples():
    return [dict(sample=f'{arm}{i}', condition=arm, subject=f'mouse{i}',
                 single_end='false', strandedness='reverse', read_length='150')
            for arm, ids in [('A', [1, 2, 3]), ('B', [1, 2]), ('C', [4, 5, 6])]
            for i in ids]


def test_two_pairs_and_unrelated_contrast_keep_the_unmatched_sample():
    contrasts = [dict(contrast='B_vs_A', treatment='B', control='A', paired='true', exclude_samples='A3'),
                 dict(contrast='A_vs_C', treatment='A', control='C', paired='false')]
    _, audit = validator.validate_contrasts(contrasts, samples(), design='~ condition', analysis_purpose='analysis')
    assert audit[0]['rmats_group1_samples'] == ['B1', 'B2']
    assert audit[0]['rmats_group2_samples'] == ['A1', 'A2']
    assert audit[0]['pair_ids'] == ['mouse1', 'mouse2']
    assert audit[0]['excluded_samples'] == ['A3']
    assert audit[1]['rmats_group1_samples'] == ['A1', 'A2', 'A3']


@pytest.mark.parametrize('excluded', ['', 'unknown', 'A3;A3', 'C4', 'A1;A3'])
def test_bad_or_implicit_exclusions_are_rejected(excluded):
    row = dict(contrast='B_vs_A', treatment='B', control='A', paired='true', exclude_samples=excluded)
    with pytest.raises(validator.InputError):
        validator.validate_contrasts([row], samples(), design='~ condition', analysis_purpose='analysis')


@pytest.mark.parametrize('header', ['A1\tA2\tA3\n', 'event\tA1\tA2\tA3\n'])
def test_suppa_removes_only_the_selected_column_preserving_native_values(tmp_path, header):
    source, dest = tmp_path/'input.psi', tmp_path/'selected.psi'
    text = header + 'e1\t0.1\tNA\t0.9\ne2\t0.8\t0.3\t0.7\n'
    source.write_text(text)
    suppa.subset(source, dest, {'A3'})
    assert dest.read_text() == header.replace('\tA3', '') + 'e1\t0.1\tNA\ne2\t0.8\t0.3\n'
    assert source.read_text() == text


def test_import_keeps_whole_study_sample_but_excludes_it_from_contrast(tmp_path):
    fixture = load('subset_import_fixture', 'tests/public_baseline/test_flow_import.py')
    from katana_next.flow.importer import import_bundle
    from katana_next.flow.model import write_json
    root, manifest, template = fixture.bundle(tmp_path)
    manifest['samples'].append(dict(sample_id='extra', condition='A', subject='mouse-extra'))
    manifest['rmats_sample_order']['A'].append('extra')
    manifest['comparisons'][0]['excluded_samples'] = ['extra']
    write_json(root/'manifest.json', manifest)
    value = json.loads(import_bundle(root, template, tmp_path/'workspace').read_text())
    assert {s['sample_id'] for s in value['samples']} == {'a', 'b', 'extra'}
    assert value['contrasts'][0]['denominator']['sample_ids'] == ['a']
    assert value['assets'][0]['metadata']['denominator_sample_ids'] == ['a']


def test_whole_study_association_carries_explicit_subset_without_removing_discovery_sample():
    from katana_cosplice.study_design import project_study_design
    records = [SimpleNamespace(sample_id=r['sample'], biological_sample_id=r['subject'], condition=r['condition'])
               for r in samples()]
    contrast = SimpleNamespace(contrast_id='B_vs_A', numerator=SimpleNamespace(sample_ids=['B1','B2']),
                               denominator=SimpleNamespace(sample_ids=['A1','A2']))
    project = SimpleNamespace(samples=records, contrasts=[contrast], metadata={'comparison_designs': {'B_vs_A': {'paired': True}}})
    policy, rows = project_study_design(project)
    assert policy['samples'] == 8
    assert rows[0]['sample_ids'] == ['B1', 'B2', 'A1', 'A2']
    assert 'A3' in [s.sample_id for s in project.samples]


def test_bundle_expands_only_the_samples_present_in_native_vectors(tmp_path):
    exporter = load('subset_bundle', 'pipeline/bin/katana_bundle.py')
    root = tmp_path/'bundle'
    root.mkdir()
    native = root/'SE.MATS.JC.txt'
    coordinates = exporter.COORDINATES['SE']
    record = dict(GeneID='G', chr='chr1', strand='+', **{key: str(i+100) for i, key in enumerate(coordinates)})
    record.update(IncLevel1='0.1,0.2', IncLevel2='0.5,0.6', IJC_SAMPLE_1='10,11', IJC_SAMPLE_2='12,13', SJC_SAMPLE_1='20,21', SJC_SAMPLE_2='22,23')
    with native.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record), delimiter='\t')
        writer.writeheader(); writer.writerow(record)
    artifacts = [dict(kind='rmats', relative_path=native.name, status='available', comparison_id='B_vs_A', mode='JC', event_type='SE')]
    comparisons = [dict(comparison_id='B_vs_A', rmats_group1='B', rmats_group2='A', excluded_samples=['A3'])]
    exporter.make_matrices(root, artifacts, comparisons, {'B':['B1','B2'], 'A':['A1','A2','A3']})
    rows = list(csv.reader((root/'matrices/B_vs_A/JC/psi.tsv').open(), delimiter='\t'))
    assert rows[0] == ['event_key','B1','B2','A1','A2']
    assert rows[1][1:] == ['0.1','0.2','0.5','0.6']


def test_full_network_uses_full_condition_sources_and_retains_all_samples(tmp_path):
    import pandas as pd
    from katana_next.full_fixture import create_full_fixture
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.flow.network_inputs import prepare_project_network
    from katana_next.flow.model import file_sha256
    from katana_cosplice.matrix import build_sample_matrices
    path = create_full_fixture(tmp_path/'study')
    _, manifest = ManifestStore().load(path)
    original = manifest.contrasts[0]
    partial = original.model_copy(update={
        'contrast_id':'subset',
        'numerator':original.numerator.model_copy(update={'sample_ids':original.numerator.sample_ids[:2]}),
        'denominator':original.denominator.model_copy(update={'sample_ids':original.denominator.sample_ids[:2]})})
    duplicate = [a.model_copy(update={'asset_id':a.asset_id+'-subset', 'metadata':{**a.metadata,
        'contrast_id':'subset', 'numerator_sample_ids':list(partial.numerator.sample_ids),
        'denominator_sample_ids':list(partial.denominator.sample_ids)}})
        for a in manifest.assets if a.kind=='rmats_result' and a.metadata.get('contrast_id')==original.contrast_id]
    manifest = manifest.model_copy(update={'contrasts':tuple(manifest.contrasts)+(partial,), 'assets':tuple(manifest.assets)+tuple(duplicate)})
    path = ManifestStore().save(manifest, tmp_path/'subset.project.json')
    config, proof, _ = prepare_project_network(path, tmp_path/'projection', {'network_scope':'whole_study','rmats_mode':'JC'}, file_sha256)
    matrix = build_sample_matrices(config)
    assert set(matrix.psi.columns) == {s.sample_id for s in manifest.samples}
    mapping = json.loads(Path(proof['projection_mapping']).read_text())
    assert len(mapping['excluded_measurement_sources']) == 5
    assert all(row['contrast_id']=='subset' for row in mapping['excluded_measurement_sources'])
    assert 'subset' not in {row['contrast_id'] for row in mapping['projection']}
