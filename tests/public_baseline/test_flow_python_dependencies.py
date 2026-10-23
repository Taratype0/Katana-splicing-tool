"""Import-only readiness catches broken interpreters before scientific work."""
import os
from pathlib import Path
import sys
import venv

import pytest

from katana_next.flow import dependency_checks, stages


def test_actual_import_probe_reports_missing_bio_before_mapp_inputs(tmp_path, monkeypatch):
    # A separate interpreter loads this deliberate broken-package fixture. It
    # cannot borrow the already imported Bio modules from the pytest process.
    shadow = tmp_path / "incomplete-environment"
    package = shadow / "Bio"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        "raise ModuleNotFoundError('Bio unavailable in dependency-test fixture')\n", encoding="utf-8")
    monkeypatch.setenv("PYTHONPATH", str(shadow) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    monkeypatch.setattr(stages, "_run", lambda *_a, **_k: pytest.fail("No scientific command may run"))
    spec = {"adapter": "mapp.prepare", "inputs": {}, "parameters": {"python": sys.executable}}
    issues = stages.readiness(spec)
    assert any("Bio.Seq" in issue and "biopython" in issue and "Katana runtime" in issue for issue in issues)
    assert all("No analysis was started" in issue for issue in issues)
    output = tmp_path / "science-output"
    with pytest.raises(ValueError, match="Required Python imports failed"):
        stages.execute_stage(spec, output)
    assert not output.exists()
    assert not (package / "__pycache__").exists()


def test_selected_worker_is_probed_not_only_host_python(tmp_path):
    worker = tmp_path / "empty-worker"
    venv.EnvBuilder(with_pip=False, system_site_packages=False).create(worker)
    executable = worker / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    assert stages._executable(str(executable), sys.executable) == str(executable.absolute())
    issues = dependency_checks.python_dependency_issues(
        {"adapter": "mapp.maei", "parameters": {"python": str(executable)}},
        resolve_executable=stages._executable,
    )
    assert any("configured worker runtime" in issue and "numpy" in issue and "pandas" in issue for issue in issues)
    assert not any("Katana runtime" in issue for issue in issues)
    assert not list(worker.glob("**/stage.receipt.json"))


def test_actual_import_failure_detects_broken_binary_dependency(tmp_path, monkeypatch):
    (tmp_path / "broken_binary_dependency.py").write_text("raise OSError('shared library missing')\n")
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    with pytest.raises(ValueError, match="OSError: shared library missing"):
        dependency_checks.require_python_imports(sys.executable, ["broken_binary_dependency"])


def test_import_probe_handles_import_noise_and_success_without_scientific_files(tmp_path, monkeypatch):
    (tmp_path / "import_noise_fixture.py").write_text("print('unrelated library startup message')\n")
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    actual = dependency_checks.require_python_imports(sys.executable, ["import_noise_fixture", "json"])
    assert Path(actual).resolve() == Path(sys.executable).resolve()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["import_noise_fixture.py"]


def test_historical_adoption_does_not_require_unrelated_compute_dependencies(monkeypatch):
    monkeypatch.setattr(dependency_checks, "require_python_imports",
                        lambda *_a, **_k: pytest.fail("Adoption must not probe compute dependencies"))
    assert dependency_checks.python_dependency_issues(
        {"adapter": "historical.adopt"}, resolve_executable=stages._executable) == []


def test_dependency_names_cannot_inject_a_script(monkeypatch):
    monkeypatch.setattr(dependency_checks.subprocess, "run",
                        lambda *_a, **_k: pytest.fail("Invalid module names must never launch a command"))
    with pytest.raises(ValueError, match="module names only"):
        dependency_checks.require_python_imports(sys.executable, ["json; open('file','w')"])
