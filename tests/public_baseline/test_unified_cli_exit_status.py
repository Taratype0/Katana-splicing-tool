"""A failed scientific job must fail the public shell command as well."""
from types import SimpleNamespace
import json
import pytest


@pytest.mark.parametrize('engine', ['cosplice', 'motiftrace'])
@pytest.mark.parametrize('status,exit_code', [('complete', 0), ('cache_hit', 0), ('failed', 1), ('cancelled', 1)])
def test_public_commands_preserve_job_exit_status(monkeypatch, capsys, engine, status, exit_code):
    class Client:
        def __init__(self, *args): pass
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def run(self, *args, **kwargs):
            return SimpleNamespace(model_dump=lambda **kwargs: {'status': status})
    if engine == 'cosplice':
        from katana_cosplice import study_cli as cli
        monkeypatch.setattr(cli, 'CoSpliceStudy', Client)
        command = ['--project', 'generated-project.json', 'run']
    else:
        from katana_motiftrace import cli
        monkeypatch.setattr(cli, 'MotifTrace', Client)
        command = ['--project', 'generated-project.json', 'run', 'consequences']
    assert cli.main(command) == exit_code
    assert json.loads(capsys.readouterr().out)['status'] == status
