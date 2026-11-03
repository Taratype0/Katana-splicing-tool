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


def test_upstream_links_are_immutable_and_per_bam_prep_inherits_single():
    policy = json.loads((ROOT / 'pipeline/assets/resource-policy.json').read_text())
    commits = {source['commit'] for source in policy['sources']}
    for record in policy['processes']:
        if 'upstream' in record:
            source = record['upstream']
            assert source['commit'] in commits
            assert re.fullmatch('[0-9a-f]{64}', source['sha256'])
            assert f"/blob/{source['commit']}/" in source['url']
    prep = next(r for r in policy['processes'] if r['process'] == 'RMATS_PREP_BAM')
    assert prep['initial'] == {'cpus':1, 'memory_gib':6, 'time_hours':4}
    workflow = (ROOT / 'pipeline/subworkflows/local/rmats.nf').read_text()
    assert 'RMATS_PREP_BAM(ch_sample_bams' in workflow
    assert 'RMATS_PREP (' not in workflow and 'RMATS_PREP_SINGLE (' not in workflow


def test_allocation_policy_does_not_change_science_or_select_a_cluster():
    text = (ROOT / 'pipeline/conf/upstream_task_resources.config').read_text()
    assignments = re.findall(r'^\s*([\w.]+)\s*=', text, re.M)
    assert assignments and set(assignments) == {'cpus', 'memory', 'time'}
