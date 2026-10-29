"""The installed wheel exports the same self-contained workflow sources."""
import json
import os
import re
from pathlib import Path
import subprocess
import sys

import pytest


def test_installed_managed_entry_has_executable_scripts_without_export():
    from katana_next.flow.managed import _pipeline_source
    main = _pipeline_source()
    if os.environ.get('KATANA_WHEEL_ONLY') == '1':
        assert main.parent.name == '_pipeline'
    if os.name == 'posix':
        for name in ['katana_differential.R', 'tximport.R', 'run_dexseq_exon.R']:
            assert os.access(main.parent / 'bin' / name, os.X_OK), name


def test_exported_pipeline_has_verified_sources_and_preserves_existing_exports(tmp_path):
    from katana_next.flow.distribution import export_pipeline
    pipeline = export_pipeline(tmp_path / "pipeline")
    assert (pipeline / "main.nf").is_file()
    assert (pipeline / "workflows/katana_downstream.nf").is_file()
    assert (pipeline / "bin/katana_bundle.py").is_file()
    result = subprocess.run([sys.executable, "-B", str(pipeline / "tools/verify_sources.py")],
                            cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["verified_source_files"] == 89
    # BioContainers publishes these fixed builds at Quay, not Docker Hub.
    for module in (pipeline / "modules").rglob("*.nf"):
        assert not re.search(r"['\"]biocontainers/", module.read_text(encoding="utf-8")), module
    with pytest.raises(FileExistsError, match="preserved"):
        export_pipeline(pipeline)


def test_flow_cli_is_batch_safe_and_available_outside_checkout(tmp_path):
    import katana_next.flow
    command = [sys.executable, *(["-I"] if os.environ.get("KATANA_WHEEL_ONLY") == "1" else []),
               "-m", "katana_next.flow.cli", "--help"]
    environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    if os.environ.get("KATANA_WHEEL_ONLY") != "1":
        environment["PYTHONPATH"] = str(Path(katana_next.flow.__file__).resolve().parents[2])
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, env=environment, timeout=30)
    assert result.returncode == 0, result.stderr
    assert "export-pipeline" in result.stdout and "collect" in result.stdout
    for action in ("run", "retain-project", "verify-retained", "clean-run"):
        checked = subprocess.run([*command[:-1], action, "--help"], cwd=tmp_path,
            capture_output=True, text=True, env=environment, timeout=30)
        assert checked.returncode == 0, checked.stderr
        assert "usage:" in checked.stdout
        if action == "run":
            assert "--keep-work" in checked.stdout and "--run-dir" in checked.stdout
