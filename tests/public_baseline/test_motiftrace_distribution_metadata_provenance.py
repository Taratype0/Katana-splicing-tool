"""A distribution-version audit cannot replace a frozen scientific identity."""
from hashlib import sha256
import json
from pathlib import Path
import tomllib


ROOT=Path(__file__).resolve().parents[2]
PREVIOUS='cd24246b872984e6dffbe35cfcd28f6787aa3fe4c066c6c265006052898fd419'
CURRENT='8d441a2bec5274e14334b2ee129bea2ca900930dd91fc0e62af6ee4f8555a0e8'
SCIENTIFIC='a4bd76153a3975213c0994e004d65d7785113a2adef3ade86ecb46e7f289a7ed'
MODES='82f67e45f464c98276868e2b15b25fb9ec683f4bc0c67611d507eeb5284d105b'


def _record():
    # Preserve the exact historical metadata proof after the explicit execution
    # capacity revision; do not relabel new source bytes as an old release.
    value=json.loads((ROOT/'tests/fixtures/frozen_motiftrace_capacity_20261017/integration.before-capacity-removal.json').read_text('utf-8'))
    revision=next(row for row in value['interface_revisions']
                  if row.get('revision')=='distribution-version-metadata-20261025')
    return value,revision


def _canonical(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def test_metadata_revision_records_both_exact_interface_identities_and_history():
    record,revision=_record()
    assert record['integration_version']==revision['integration_version']=='1.0.0.post5'
    assert len(revision['changes'])==1
    change=revision['changes'][0]
    assert change['path']=='pyproject.toml'
    assert change['previous_sha256']==PREVIOUS
    assert change['current_sha256']==record['interface_source_sha256']['pyproject.toml']==CURRENT
    assert change['change_kind']=='distribution_version_only'
    assert change['changed_toml_fields']==['project.version']
    assert any(row.get('revision')=='explicit-network-threadpool-runtime-dependency-20261024'
               and row['changes'][0]['current_sha256']==PREVIOUS for row in record['interface_revisions'])


def test_only_declared_distribution_version_reconstructs_the_previous_exact_bytes():
    _,revision=_record()
    change=revision['changes'][0]
    current=(ROOT/'tests/fixtures/frozen_motiftrace_capacity_20261017/pyproject.before-capacity-removal.toml').read_bytes()
    old_token=('version = "'+change['previous_distribution_version']+'"').encode()
    new_token=('version = "'+change['current_distribution_version']+'"').encode()
    assert current.count(new_token)==1
    previous=current.replace(new_token,old_token,1)
    assert sha256(current).hexdigest()==CURRENT
    assert sha256(previous).hexdigest()==PREVIOUS
    old_doc=tomllib.loads(previous.decode('utf-8'))
    new_doc=tomllib.loads(current.decode('utf-8'))
    assert old_doc['project'].pop('version')==change['previous_distribution_version']
    assert new_doc['project'].pop('version')==change['current_distribution_version']
    assert old_doc==new_doc


def test_packaging_revision_preserves_all_scientific_source_and_mode_fingerprints():
    record,revision=_record()
    assert record['schema']=='katana.motiftrace_integration/v1'
    assert record['scientific_release']=='1.0.0-frozen.20261017'
    assert record['frozen_source_archive_sha256']=='fe3f0220d16940a7d90f06a22320d0dee7aa8ebc3d7186f76a131a6d351a5778'
    assert record['frozen_release_manifest_sha256']=='d8329931ecb1f453897cb2f8d157dff4fe046c7746d3e10805cfdc43db98357e'
    assert revision['scientific_release_unchanged'] is True
    assert revision['scientific_calculation_changes'] is revision['scientific_parameter_changes'] is False
    assert not record['scientific_calculation_changes'] and not record['scientific_parameter_changes']
    assert len(record['scientific_sources'])==23
    assert _canonical(record['scientific_sources'])==revision['scientific_sources_manifest_sha256']==SCIENTIFIC
    assert len(record['modes'])==6
    assert _canonical(record['modes'])==revision['workflow_modes_manifest_sha256']==MODES
    for name,values in record['scientific_sources'].items():
        snapshot=ROOT/'tests/fixtures/frozen_motiftrace_capacity_20261017'/Path(name).name
        source=snapshot if snapshot.is_file() else ROOT/name
        if name == 'katana_next/science/motiftrace_consequences.py':
            source=ROOT/'tests/fixtures/performance_reference_20261026/motiftrace_consequences.py'
        assert sha256(source.read_bytes()).hexdigest()==values['current_sha256']


def test_capacity_revision_preserves_historical_identity_and_records_current_sources():
    historical,_=_record()
    record=json.loads((ROOT/'docs/motiftrace/integration.v1.json').read_text('utf-8'))
    revision=next(r for r in record['interface_revisions']
                  if r.get('revision')=='remove-fixed-motiftrace-capacity-limits-20261026')
    assert record['integration_version']=='1.0.0.post7'
    assert revision['integration_version']=='1.0.0.post6'
    oldpath=ROOT/revision['previous_record_snapshot']
    assert sha256(oldpath.read_bytes()).hexdigest()==revision['previous_integration_record_sha256']
    assert json.loads(oldpath.read_text())==historical
    old_toml=ROOT/'tests/fixtures/frozen_motiftrace_capacity_20261017/pyproject.before-capacity-removal.toml'
    changes={r['path']:r for r in revision['changes']}
    assert changes['pyproject.toml']['previous_sha256']==sha256(old_toml.read_bytes()).hexdigest()
    packaging = next(r for r in record['interface_revisions']
                     if r.get('revision') == 'distribution-version-metadata-20261028')['changes'][0]
    current = (ROOT/'pyproject.toml').read_bytes()
    old_token = ('version = "'+packaging['previous_distribution_version']+'"').encode()
    new_token = ('version = "'+packaging['current_distribution_version']+'"').encode()
    assert current.count(new_token) == 1
    previous = current.replace(new_token, old_token, 1)
    assert changes['pyproject.toml']['current_sha256'] == packaging['previous_sha256'] == sha256(previous).hexdigest()
    assert packaging['current_sha256'] == record['interface_source_sha256']['pyproject.toml'] == sha256(current).hexdigest()
    old_doc, new_doc = tomllib.loads(previous.decode()), tomllib.loads(current.decode())
    assert old_doc['project'].pop('version') == packaging['previous_distribution_version']
    assert new_doc['project'].pop('version') == packaging['current_distribution_version']
    assert old_doc == new_doc
    assert tomllib.loads((ROOT/'pyproject.toml').read_text())['project']['version']=='2.0.0b5.dev20261026'
    assert record['execution_adapter']['fixed_capacity_limits_enforced'] is False
    assert record['execution_adapter']['max_tests_default'] is None
    assert record['execution_adapter']['max_matrix_cells_default'] is None
    for name,value in record['scientific_sources'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest()==value['current_sha256']
