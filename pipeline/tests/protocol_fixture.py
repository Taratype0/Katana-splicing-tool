#!/usr/bin/env python3
"""Nextflow transport fixture only: does not calculate biological results.

Used to verify real native task fan-out, dependency staging, empty levels,
failure propagation and -resume. Never selected by a production profile.
"""
from __future__ import annotations
import argparse
import json
import os
import time
from pathlib import Path


def write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    for command_name in ('preflight', 'import'):
        entry = sub.add_parser(command_name)
        entry.add_argument('--project-template', type=Path, required=True)
        entry.add_argument('--config', type=Path, required=True)
        if command_name == 'preflight':
            entry.add_argument('--out', type=Path, required=True)
        else:
            entry.add_argument('--bundle', type=Path, required=True)
            entry.add_argument('--outdir', type=Path, required=True)
    plan = sub.add_parser("plan")
    plan.add_argument("--project", type=Path, required=True)
    plan.add_argument("--config", type=Path, required=True)
    plan.add_argument("--out", type=Path, required=True)
    plan.add_argument("--require-ready", action="store_true")
    task = sub.add_parser("task")
    task.add_argument("--plan", type=Path, required=True)
    task.add_argument("--project", type=Path, required=True)
    task.add_argument("--task-id", required=True)
    task.add_argument("--outdir", type=Path, required=True)
    task.add_argument("--dependency-receipts", nargs="*", type=Path, default=[])
    runtime = sub.add_parser("runtime-check")
    runtime.add_argument("--plan", type=Path, required=True)
    runtime.add_argument("--project", type=Path, required=True)
    runtime.add_argument("--task-id", required=True)
    runtime.add_argument("--runtime-label")
    runtime.add_argument("--out", type=Path, required=True)
    collect = sub.add_parser("collect")
    collect.add_argument("--plan", type=Path, required=True)
    collect.add_argument("--task-results", nargs="+", type=Path, required=True)
    collect.add_argument("--outdir", type=Path, required=True)
    report = sub.add_parser("report")
    report.add_argument("--collection", type=Path, required=True)
    report.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'preflight':
        assert args.project_template.is_file() and args.config.is_file()
        write(args.out, {'fixture': True, 'ready': True})
    elif args.command == 'import':
        assert (args.bundle / 'fixture.txt').read_text() == 'generated transport only'
        write(args.outdir / 'project.json', {'fixture': True})
    elif args.command == "plan":
        document = json.loads(args.config.read_text())
        assert args.project.is_file()
        if document.get("schema") == "katana.flow.config/v1":
            # Native flat-entry transport test: parse real generated settings,
            # but perform no biological computation in this explicit test CLI.
            assert document["network"]["scope"] == "contrast"
            assert document["parameters"]["differential.deg"]["fdr"] > 0
            document["tasks"] = [
                {"id": "fixture-audit", "execution_group": "audit", "workflow_id": "fixture", "dependencies": [], "delay": 0},
                {"id": "fixture-quality", "execution_group": "quality", "workflow_id": "fixture", "dependencies": ["fixture-audit"], "delay": 0},
                {"id": "fixture-as", "execution_group": "as", "workflow_id": "fixture", "dependencies": ["fixture-quality"], "delay": 0}]
        planned = {"schema": "test.nextflow_transport/v1", "tasks": document["tasks"]}
        if document.get("task_cache_scope") == "ancestry-v1":
            import importlib.util
            source = Path(__file__).resolve().parents[2] / "katana_next/flow/model.py"
            spec = importlib.util.spec_from_file_location("fixture_plan_model", source)
            model = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(model)
            planned = model.seal_plan({**planned, "task_cache_scope": "ancestry-v1",
                "project_manifest": "project.json", "project_manifest_sha256": "generated-fixture",
                "report_options": document.get("report_options", {})})
        write(args.out, planned)
    elif args.command == "runtime-check":
        spec = next(item for item in json.loads(args.plan.read_text())["tasks"] if item["id"] == args.task_id)
        write(args.out, {"task_id": args.task_id, "runtime_label": args.runtime_label,
                        "status": "failed" if spec.get("runtime_fail") else "ready", "scientific_execution": False,
                        "thread_budget": {key: os.environ.get(key) for key in
                            ("KATANA_CPUS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                             "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")}})
        if spec.get("runtime_fail"):
            raise SystemExit("intentional transport-fixture runtime failure")
    elif args.command == "task":
        document = json.loads(args.plan.read_text())
        spec = next(item for item in document["tasks"] if item["id"] == args.task_id)
        probe = json.loads(Path("runtime-" + args.task_id + ".json").read_text())
        assert probe["task_id"] == args.task_id and probe["status"] == "ready"
        seen = {json.loads(path.read_text())["task_id"]: path for path in args.dependency_receipts}
        assert set(spec["dependencies"]).issubset(seen)
        # A staged dependency needs its output payload, not only a detached receipt.
        for dep in spec["dependencies"]:
            assert (seen[dep].parent / "payload.txt").read_text() == dep
        assert args.project.is_file()
        if spec.get("fail"):
            raise SystemExit("intentional transport-fixture failure")
        custom = spec.get('parameters', {}).get('custom_motif_path')
        if custom:
            staged = list(Path('runtime-data').glob('*/*'))
            assert len(staged) == 1 and staged[0].read_bytes() == Path(custom).read_bytes()
        started = time.time()
        time.sleep(float(spec.get("delay", 1.0)))
        args.outdir.mkdir(parents=True)
        (args.outdir / "payload.txt").write_text(args.task_id)
        write(args.outdir / "task.receipt.json", {"task_id": args.task_id, "status": "completed",
            "started": started, "finished": time.time(), "pid": os.getpid(),
            "scientific_execution": False})
    elif args.command == "collect":
        planned = {task["id"] for task in json.loads(args.plan.read_text())["tasks"]}
        observed = {json.loads((path / "task.receipt.json").read_text())["task_id"] for path in args.task_results}
        assert observed == planned
        write(args.outdir / "collection.receipt.json", {"planned": sorted(planned), "observed": sorted(observed),
            "scientific_execution": False})
    else:
        assert (args.collection / "collection.receipt.json").is_file()
        args.outdir.mkdir(parents=True)
        (args.outdir / "REPORT.html").write_text("<html><body>Transport fixture. No scientific report.</body></html>")
        write(args.outdir / "report-manifest.json", {"fixture": True, "scientific_execution": False})
        write(args.outdir / "flow.acceptance.json", {
            "schema": "test.nextflow_transport_acceptance/v1",
            "status": "passed", "fixture": True, "scientific_execution": False,
        })


if __name__ == "__main__":
    main()
