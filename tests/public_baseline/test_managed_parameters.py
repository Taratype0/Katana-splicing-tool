"""Unified public options and storage handoffs, using generated protocol data."""
import importlib.util
import json
from pathlib import Path
import re

import pytest

from katana_next.flow import managed
from katana_next.flow.managed_parameters import schema
from katana_next.flow.model import read_json, write_json
from test_flow_managed_run import lifecycle

ROOT = Path(__file__).resolve().parents[2]


def simple(lifecycle, monkeypatch):
    argv, inputs, calls, fail = lifecycle
    args = managed.parser().parse_args(argv)
    args.project_template = args.flow_config = None
    args.native_parameters = {'skip_' + key: True for key in (
        'network', 'mapp', 'motiftrace', 'motif_binding', 'enrichment', 'rbp_identity')}
    args.native_parameters.update(deg_fdr=.01, species='Mus musculus')
    spec = importlib.util.spec_from_file_location('managed_fixture_configure', ROOT / 'pipeline/bin/katana_configure.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    native = managed._native_phase
    commands = []
    def configured(command, work, scope, stage, **kwargs):
        commands.append(command)
        if stage == 'upstream':
            assert '--project_template' not in command and '--flow_config' not in command
            params = read_json(command[command.index('-params-file') + 1])
            output = Path(command[command.index('--outdir') + 1])
            module.generate(params, [], output / 'pipeline_info' / 'configuration')
        return native(command, work, scope, stage, **kwargs)
    monkeypatch.setattr(managed, '_native_phase', configured)
    return args, calls, fail, commands


def test_flat_entry_uses_same_generator_cleans_and_preserves_outputs(lifecycle, monkeypatch):
    args, calls, _, commands = simple(lifecycle, monkeypatch)
    managed.run(args)
    out = args.outdir
    scope = out.with_name(out.name + '.run')
    assert calls == [('upstream', False), ('downstream', False)]
    assert not (scope / 'work/upstream').exists()
    assert not (scope / 'work/downstream').exists()
    config = out / 'upstream-support/pipeline_info/configuration'
    assert read_json(config / 'flow.json')['parameters']['differential.deg']['fdr'] == .01
    assert read_json(config / 'project-template.json')['metadata']['species'] == 'Mus musculus'
    assert (out / 'upstream-support/models/fit.rds').is_file()
    assert (out / 'report/REPORT.zip').is_file()
    assert '--flow_config' in commands[1]
    args.resume = True
    assert managed.run(args)['status'] == 'completed'
    assert len(calls) == 2


def test_flat_failure_resume_reuses_upstream_and_configuration(lifecycle, monkeypatch):
    args, calls, fail, _ = simple(lifecycle, monkeypatch)
    fail['downstream'] = True
    with pytest.raises(RuntimeError, match='downstream failure'):
        managed.run(args)
    fail['downstream'] = False
    args.resume = True
    assert managed.run(args)['status'] == 'completed'
    assert calls == [('upstream', False), ('downstream', False), ('downstream', True)]


def test_changed_generated_configuration_is_not_adopted_on_resume(lifecycle, monkeypatch):
    args, calls, fail, _ = simple(lifecycle, monkeypatch)
    fail['downstream'] = True
    with pytest.raises(RuntimeError, match='downstream failure'):
        managed.run(args)
    config = args.outdir.with_name('results.run') / 'work/downstream/upstream-output/pipeline_info/configuration/flow.json'
    data = read_json(config)
    data['parameters']['differential.deg']['fdr'] = .5
    write_json(config, data)
    args.resume = True
    fail['downstream'] = False
    with pytest.raises(ValueError, match='configuration changed'):
        managed.run(args)
    assert calls == [('upstream', False), ('downstream', False)]


def test_flat_keep_work_then_explicit_cleanup(lifecycle, monkeypatch):
    args, _, _, _ = simple(lifecycle, monkeypatch)
    args.keep_work = True
    assert managed.run(args)['status'] == 'completed'
    scope = args.outdir.with_name('results.run')
    assert all((scope / 'work' / stage).is_dir() for stage in ('upstream', 'downstream'))
    assert managed.clean_completed_run(scope, apply=True)['status'] == 'completed'
    assert all(not (scope / 'work' / stage).exists() for stage in ('upstream', 'downstream'))
    assert (args.outdir / 'report/REPORT.zip').is_file()


def test_flat_failed_report_does_not_clean_downstream(lifecycle, monkeypatch):
    args, _, fail, _ = simple(lifecycle, monkeypatch)
    fail['acceptance'] = True
    with pytest.raises(ValueError, match='did not pass'):
        managed.run(args)
    assert (args.outdir.with_name('results.run') / 'work/downstream').is_dir()


def test_flat_downstream_only_generates_config_natively(lifecycle, monkeypatch):
    argv, inputs, calls, _ = lifecycle
    commands = []
    native = managed._native_phase
    def capture(command, *args, **kwargs):
        commands.append(command)
        return native(command, *args, **kwargs)
    monkeypatch.setattr(managed, '_native_phase', capture)
    args = managed.parser().parse_args(['--mode', 'downstream-only', '--project', str(inputs / 'project.json'),
        '--outdir', argv[1], '--skip_network', '--deg_fdr', '0.02'])
    assert managed.run(args)['status'] == 'completed'
    assert calls == [('downstream', False)]
    assert '--flow_config' not in commands[0]
    forwarded = read_json(commands[0][commands[0].index('-params-file') + 1])
    assert forwarded == {'skip_network': True, 'deg_fdr': .02}


def test_native_yaml_and_cli_share_settings_with_cli_precedence(tmp_path):
    params = tmp_path / 'params.yaml'
    params.write_text('mode: downstream-only\nproject: /data/project.json\noutdir: /data/results\n'
                      'deg_fdr: 0.05\nskip_network: true\ncontrast_ids: [B_vs_A, C_vs_A]\n'
                      'katana_stage_containers: {mapp: example:1}\nstorage_mode: retain\n')
    args = managed.parser().parse_args(['-params-file', str(params), '--deg_fdr', '.01',
        '--skip_network', 'false', '-profile', 'slurm,apptainer', '-c', 'site.config', '-resume'])
    assert args.mode == 'downstream-only' and args.outdir == Path('/data/results').resolve()
    assert args.profile == 'slurm,apptainer' and args.resume
    assert args.native_parameters['deg_fdr'] == .01 and args.native_parameters['skip_network'] is False
    assert args.native_parameters['contrast_ids'] == ['B_vs_A', 'C_vs_A']
    assert args.native_parameters['katana_stage_containers'] == {'mapp': 'example:1'}
    assert 'mode' not in args.native_parameters and 'storage_mode' not in args.native_parameters


@pytest.mark.parametrize('text', ['skip_network: maybe', 'deg_fdr: .nan', 'read_length: 2.5',
                                  'resource_metadata: []', 'managed_preflight: true', 'unknown_parameter: 1'])
def test_invalid_file_values_rejected_before_launch(tmp_path, text):
    path = tmp_path / 'parameters.yaml'
    path.write_text('outdir: results\n' + text + '\n')
    with pytest.raises(SystemExit):
        managed.parser().parse_args(['-params-file', str(path)])


def test_all_native_public_parameter_names_are_supported():
    config = (ROOT / 'pipeline/nextflow.config').read_text().split('params {', 1)[1].split('\n}', 1)[0]
    names = set(re.findall(r'^\s*(\w+)\s*=', config, re.M)) - {'help', 'managed_preflight'}
    assert set(schema()) == names
    parser = managed.parser()
    assert names - {'storage_mode', 'publish_dir_mode'} <= {a.dest for a in parser._actions}


def test_missing_native_configuration_never_authorizes_cleanup(lifecycle):
    argv, _, _, _ = lifecycle
    args = managed.parser().parse_args(argv)
    args.project_template = args.flow_config = None
    with pytest.raises(ValueError, match='Missing published native configuration'):
        managed.run(args)
    assert (args.outdir.with_name('results.run') / 'work/upstream').is_dir()
