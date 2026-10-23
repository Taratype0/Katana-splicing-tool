"""The installed wheel exports the same self-contained workflow sources."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


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
