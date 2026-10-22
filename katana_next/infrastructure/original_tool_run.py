"""Auditable subprocess execution; completion requires verified output artifacts."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Mapping, Sequence


def file_record(path: Path) -> dict:
    path = path.resolve(strict=True)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Expected a non-empty file: {path}")
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def run_original_tool(
    *, output: Path, repository: Path, commands: Sequence[Sequence[str]],
    inputs: Sequence[Path], expected_outputs: Sequence[str], boundary: str,
    cwd: Path | None = None, environment: Mapping[str, str] | None = None,
    timeout_seconds: float = 3600,
) -> dict:
    """Run a clean original checkout without a shell, retaining failures and logs.

    This verifies execution provenance, not numerical or biological validity.
    Callers must separately validate the scientific contents of the artifacts.
    """
    if not commands or not expected_outputs or timeout_seconds <= 0:
        raise ValueError("Commands, expected outputs and positive timeout required")
    output, repository = output.resolve(), repository.resolve(strict=True)
    outputs = [(output / name).resolve() for name in expected_outputs]
    if any(not path.is_relative_to(output) for path in outputs):
        raise ValueError("Expected outputs must stay inside the run directory")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use a fresh output directory; existing runs are immutable")
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=repository, text=True).strip()
    if dirty:
        raise ValueError("Original tracked source has modifications")
    input_records = [file_record(Path(path)) for path in inputs]
    output.mkdir(parents=True, exist_ok=True)
    parameters = {
        "original_repository": str(repository), "original_commit": commit,
        "commands": [list(command) for command in commands],
        "cwd": str(cwd or repository), "inputs": input_records,
        "expected_outputs": list(expected_outputs), "boundary": boundary,
        "timeout_seconds_per_command": timeout_seconds,
        "environment_overrides": dict(environment or {}),
        "executor_python": sys.version,
        "executor_source": file_record(Path(__file__)),
    }
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    started = time.perf_counter()
    steps = []
    try:
        for index, command in enumerate(commands):
            tick = time.perf_counter()
            with (output / f"step_{index + 1:02d}.log").open("wb") as log:
                result = subprocess.run(
                    list(command), cwd=cwd or repository,
                    env={**os.environ, **dict(environment or {})},
                    stdout=log, stderr=subprocess.STDOUT,
                    timeout=timeout_seconds, check=True,
                )
            steps.append({"step": index + 1, "returncode": result.returncode,
                          "seconds": time.perf_counter() - tick})
        records = [file_record(path) for path in outputs]
        if [file_record(Path(path)) for path in inputs] != input_records:
            raise RuntimeError("Input changed during execution")
        if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip() != commit:
            raise RuntimeError("Original checkout changed during execution")
        if subprocess.check_output(["git", "diff", "HEAD", "--"], cwd=repository):
            raise RuntimeError("Original tracked source changed during execution")
        status = {"status": "execution_complete", "steps": steps,
                  "seconds": time.perf_counter() - started, "outputs": records,
                  "scientific_validation": "not_asserted_by_executor"}
        (output / "completed.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
        return status
    except Exception as exc:
        (output / "failed.json").write_text(json.dumps({
            "status": "failed", "error": str(exc), "steps": steps,
            "seconds": time.perf_counter() - started}, indent=2), encoding="utf-8")
        raise
