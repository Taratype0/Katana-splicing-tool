"""The resource audit must cover current modules and preserve upstream identity."""
import importlib.util
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def test_inventory_covers_all_declarations_and_labels():
    spec = importlib.util.spec_from_file_location('resource_probe', ROOT / 'pipeline/tools/check_upstream_resources.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    declared = module.declared_processes()
    policy = json.loads((ROOT / 'pipeline/assets/resource-policy.json').read_text())
    records = {r['process']: r for r in policy['processes']}
    assert len(records) == len(policy['processes'])
    assert records.keys() == declared.keys()
    for name, (path, labels) in declared.items():
        assert records[name]['labels'] == labels
        assert (ROOT / records[name]['source']).resolve() == path.resolve()
        assert records[name]['rationale']


def test_upstream_links_are_immutable_and_granularity_exceptions_explicit():
    policy = json.loads((ROOT / 'pipeline/assets/resource-policy.json').read_text())
    commits = {source['commit'] for source in policy['sources']}
    for record in policy['processes']:
        if 'upstream' in record:
            source = record['upstream']
            assert source['commit'] in commits
            assert re.fullmatch('[0-9a-f]{64}', source['sha256'])
            assert f"/blob/{source['commit']}/" in source['url']
    prep = [r for r in policy['processes'] if r['process'].startswith('RMATS_PREP')]
    assert len(prep) == 2
    assert all(r['disposition'] == 'granularity exception retained' for r in prep)
    assert all('one BAM' in r['rationale'] for r in prep)


def test_allocation_policy_does_not_change_science_or_select_a_cluster():
    text = (ROOT / 'pipeline/conf/upstream_task_resources.config').read_text()
    assignments = re.findall(r'^\s*([\w.]+)\s*=', text, re.M)
    assert assignments and set(assignments) == {'cpus', 'memory', 'time'}
