#!/usr/bin/env python3
"""Run real Nextflow with synthetic transport fixtures, never biological fits."""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import tempfile
from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[1]


def run(command: list[str], cwd: Path, label: str, *, successful: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=300)
    (cwd / f"{label}.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (cwd / f"{label}.stderr.txt").write_text(result.stderr, encoding="utf-8")
    if successful != (result.returncode == 0):
        raise AssertionError(f"{label}: unexpected exit {result.returncode}\n{result.stdout[-5000:]}\n{result.stderr[-5000:]}")
    return result


def tasks() -> list[dict]:
    def task(identifier: str, group: str, dependencies: list[str], delay: int = 1) -> dict:
        return {"id": identifier, "execution_group": group, "workflow_id": "fixture",
                "dependencies": dependencies, "delay": delay}
    return [task("alpha", "audit", []), task("beta", "audit", [], 8),
            task("gamma", "quality", ["alpha"]), task("delta", "as", ["gamma"]),
            task("epsilon", "quality", ["beta"]), task("zeta", "as", ["epsilon"]),
            task("eta", "motiftrace", ["delta"]), task("theta", "motiftrace", ["delta"]), task("omega", "network_fit", ["alpha"]),
            task("network-summary", "network_followup", ["omega"])]


def make_case(root: Path, *, failure: bool = False, runtime_failure: bool = False) -> Path:
    root.mkdir(parents=True)
    (root / "project").mkdir()
    (root / "configuration").mkdir()
    (root / "project/project.json").write_text("{}\n", encoding="utf-8")
    planned = tasks()
    models = root / 'external-models'
    models.mkdir()
    model = models / 'custom.tsv'
    model.write_text('model\tconsensus\nfixture\tACGU\n')
    next(task for task in planned if task['id'] == 'theta')['parameters'] = {'custom_motif_path': str(model.resolve())}
    if failure:
        planned[1]["fail"] = True
    if runtime_failure:
        planned[1]["runtime_fail"] = True
    (root / "configuration/flow.json").write_text(json.dumps({"tasks": planned}), encoding="utf-8")
    return root


def command(nextflow: list[str]) -> list[str]:
    return [*nextflow, "run", str(PIPELINE / "main.nf"), "--mode", "downstream-only",
            "--project", "project/project.json", "--flow_config", "configuration/flow.json",
            "--flow_executable", str(PIPELINE / "tests/protocol_fixture.py"), "--outdir", "output",
            "-c", str(PIPELINE / "tests/transport.config"), "-ansi-log", "false", "-work-dir", "work"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_mutually_exclusive_group()
    commands.add_argument("--nextflow-command-json", help="JSON argv; no shell evaluation")
    commands.add_argument("--nextflow-command-file", type=Path, help="UTF-8 JSON argv file")
    parser.add_argument("--outdir", type=Path, help="New isolated evidence directory (default temporary)")
    args = parser.parse_args()
    encoded = args.nextflow_command_file.read_text(encoding="utf-8") if args.nextflow_command_file else args.nextflow_command_json or '["nextflow"]'
    nextflow = json.loads(encoded)
    if not isinstance(nextflow, list) or not nextflow or not all(isinstance(item, str) for item in nextflow):
        parser.error("--nextflow-command-json must be a nonempty string array")
    root = args.outdir.resolve() if args.outdir else Path(tempfile.mkdtemp(prefix="katana-nextflow-transport-"))
    if args.outdir:
        root.mkdir(parents=True, exist_ok=False)
    fixture = PIPELINE / "tests/protocol_fixture.py"
    fixture.chmod(fixture.stat().st_mode | 0o111)

    case = make_case(root / "success")
    run(command(nextflow), case, "first")
    receipts = {path.parent.name: json.loads(path.read_text())
                for path in (case / "output/downstream").rglob("task.receipt.json")}
    assert len(receipts) == len(tasks()), receipts.keys()
    probes = {json.loads(path.read_text())["task_id"]: json.loads(path.read_text())
              for path in (case / "output/downstream").glob("runtime-*.json")}
    assert set(probes) == set(receipts)
    assert probes["omega"]["runtime_label"] == "fixture-network-image"
    assert probes["network-summary"]["runtime_label"] == "fixture-followup-image"
    assert probes["alpha"]["runtime_label"] == "configured-host"
    assert receipts["gamma"]["finished"] < receipts["beta"]["finished"], "A false whole-group barrier delayed gamma"
    assert receipts["omega"]["finished"] < receipts["beta"]["finished"], "Network was blocked by unrelated beta"
    assert receipts["delta"]["started"] >= receipts["gamma"]["finished"]
    assert (case / "output/report/REPORT.html").is_file()
    before = {path.relative_to(case).as_posix(): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
              for path in (case / "work").rglob("task.receipt.json")}
    old_traces = set((case / "output/pipeline_info").glob("trace-*.tsv"))
    run([*command(nextflow), "-resume"], case, "resume")
    traces = set((case / "output/pipeline_info").glob("trace-*.tsv")) - old_traces
    assert len(traces) == 1
    with next(iter(traces)).open(encoding="utf-8") as handle:
        statuses = [row["status"] for row in csv.DictReader(handle, delimiter="\t")]
    assert statuses and all(status == "CACHED" for status in statuses), statuses
    after = {path.relative_to(case).as_posix(): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
             for path in (case / "work").rglob("task.receipt.json")}
    assert after == before, "A completed worker ran again during resume"

    failed = make_case(root / "failure", failure=True)
    run(command(nextflow), failed, "required-task-failure", successful=False)
    assert not (failed / "output/report/REPORT.html").exists(), "A failed required task produced a report"

    failed_runtime = make_case(root / "runtime-failure", runtime_failure=True)
    run(command(nextflow), failed_runtime, "task-runtime-failure", successful=False)
    assert not list((failed_runtime / "work").glob("**/task-results/beta/task.receipt.json")), "A failed runtime started its calculation"
    assert not (failed_runtime / "output/report/REPORT.html").exists()
    absent_image = make_case(root / "missing-image")
    rejected = run([*command(nextflow), "-profile", "docker"], absent_image, "missing-container", successful=False)
    assert "host fallback is disabled" in rejected.stdout + rejected.stderr
    assert not list((absent_image / "work").glob("**/task.receipt.json"))

    raw_without_image = root / "raw-missing-dexseq-image"
    raw_without_image.mkdir()
    rejected = run([*nextflow, "run", str(PIPELINE / "main.nf"), "--mode", "upstream-only",
                    "--input", "absent.csv", "--contrasts", "absent-contrasts.csv",
                    "--fasta", "absent.fa", "--gtf", "absent.gtf", "-profile", "docker",
                    "-ansi-log", "false"], raw_without_image, "missing-dexseq-container", successful=False)
    assert "requires --dexseq_tools_container" in rejected.stdout + rejected.stderr
    assert not list((raw_without_image / "work").glob("**/.command.sh"))

    # Actual stdlib GTF filtering in its own native Nextflow process.
    reference = root / "reference"
    reference.mkdir()
    spec = importlib.util.spec_from_file_location("reference_fixture", PIPELINE / "tests/test_reference_filter.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.fixture_gtf(reference / "fixture.gtf")
    run([*nextflow, "run", str(PIPELINE / "tests/reference_filter.nf"),
         "--fixture_gtf", "fixture.gtf", "--filter_script", str(PIPELINE / "bin/kanata_filter_analysis_gtf.py"),
         "-ansi-log", "false", "-work-dir", "work"], reference, "native-reference-filter")
    filtering = [json.loads(path.read_text()) for path in (reference / "work").rglob("filtering.receipt.json")]
    assert len(filtering) == 1 and filtering[0]["excluded_genes"] == 3
    receipt = {"schema": "katana.native_nextflow_transport_check/v1", "status": "passed",
               "scientific_model_execution": False, "native_transport_tasks": len(tasks()),
               "dependency_exact_parallelism": True, "required_failure_blocks_report": True,
               "runtime_probe_precedes_calculation": True, "stage_container_dispatch": True,
               "container_execution_tested": False, "missing_container_blocks_host_fallback": True,
               "missing_dexseq_image_blocks_raw_entry_before_submission": True,
               "resume_all_cached": True, "completed_worker_bytes_and_mtimes_unchanged": True,
               "external_custom_motif_file_staged": True,
               "native_reference_filter": filtering[0]}
    (root / "acceptance.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "passed_native_transport_only", "receipt": str(root / "acceptance.receipt.json")}))


if __name__ == "__main__":
    main()
