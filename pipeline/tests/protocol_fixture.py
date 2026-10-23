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
    collect = sub.add_parser("collect")
    collect.add_argument("--plan", type=Path, required=True)
    collect.add_argument("--task-results", nargs="+", type=Path, required=True)
    collect.add_argument("--outdir", type=Path, required=True)
    report = sub.add_parser("report")
    report.add_argument("--collection", type=Path, required=True)
    report.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "plan":
        document = json.loads(args.config.read_text())
        assert args.project.is_file()
        write(args.out, {"schema": "test.nextflow_transport/v1", "tasks": document["tasks"]})
    elif args.command == "task":
        document = json.loads(args.plan.read_text())
        spec = next(item for item in document["tasks"] if item["id"] == args.task_id)
        seen = {json.loads(path.read_text())["task_id"]: path for path in args.dependency_receipts}
        assert set(spec["dependencies"]).issubset(seen)
        # A staged dependency needs its output payload, not only a detached receipt.
        for dep in spec["dependencies"]:
            assert (seen[dep].parent / "payload.txt").read_text() == dep
        assert args.project.is_file()
        if spec.get("fail"):
            raise SystemExit("intentional transport-fixture failure")
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


if __name__ == "__main__":
    main()
