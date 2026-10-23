#!/usr/bin/env python3
"""Check the actual downstream failure policy with tiny shell-only tasks.

Uses the repository's Nextflow config, no sequencing or study inputs. One branch
fails while two independent branches run; a second explicit -resume repairs
only the failed branch and must reuse both successful branches unchanged.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys


PIPELINE = Path(__file__).resolve().parents[1]
MAIN = r'''
nextflow.enable.dsl = 2
params.fail_code = 23
process FAIL_BRANCH {
    label 'katana_downstream'
    input:
    val spec
    output:
    path 'failed-branch.txt'
    script:
    """
    date +%s%N > started.time
    sleep 3
    date +%s%N > exit.time
    if [ ${params.fail_code} -ne 0 ]; then exit ${params.fail_code}; fi
    printf 'repaired fixture branch\n' > failed-branch.txt
    """
}
process SURVIVE_BRANCH {
    tag "${spec.name}"
    label 'katana_downstream'
    input:
    val spec
    output:
    path "${spec.name}.txt"
    script:
    """
    date +%s%N > started.time
    sleep 12
    printf '${spec.name} completed independently\n' > ${spec.name}.txt
    date +%s%N > finished.time
    """
}
process BLOCKED_REPORT {
    label 'katana_report'
    input:
    path failed_branch
    path successful_branches
    output:
    path 'REPORT.ok'
    script:
    """
    printf 'all generated prerequisites completed\n' > REPORT.ok
    """
}
workflow {
    resources = [cpus: 1, memory: '64 MB', time: '1 min']
    failed = FAIL_BRANCH([name: 'failure', resources: resources])
    passed = SURVIVE_BRANCH(Channel.of([name: 'alpha', resources: resources], [name: 'beta', resources: resources]))
    BLOCKED_REPORT(failed, passed.collect())
}
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def trace(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def execute(command, directory, label, timeout):
    with (directory / f"{label}.stdout.txt").open("wb") as stdout, (directory / f"{label}.stderr.txt").open("wb") as stderr:
        result = subprocess.run(command, cwd=directory, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--nextflow", default="nextflow", help="Executable/path to the Nextflow launcher")
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()
    destination = args.outdir.expanduser().resolve()
    if destination == PIPELINE.parent or destination.is_relative_to(PIPELINE.parent):
        parser.error("Choose a new evidence directory outside the checkout")
    destination.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "katana.actual_nextflow_finish_acceptance/v1", "status": "running",
               "started_utc": datetime.now(timezone.utc).isoformat(),
               "scientific_computation": False, "study_inputs_used": False,
               "production_config_sha256": sha(PIPELINE / "nextflow.config"),
               "documentation": "https://docs.seqera.io/nextflow/reference/process/directives/error-strategy"}
    def save():
        (destination / "acceptance.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    try:
        save()
        (destination / "main.nf").write_text(MAIN, encoding="utf-8")
        config = str(PIPELINE / "nextflow.config").replace("\\", "/").replace("'", "\\'")
        # Import the actual policy; only the synthetic report's resource request
        # is reduced. No strategy is duplicated or overridden in this fixture.
        (destination / "nextflow.config").write_text(
            f"includeConfig '{config}'\n"
            "process { withName: BLOCKED_REPORT { cpus = 1; memory = '64 MB'; time = '1 min' } }\n"
            "trace.fields = 'task_id,hash,name,status,exit,workdir'\n"
            "report.enabled = false\ntimeline.enabled = false\ndag.enabled = false\n", encoding="utf-8")
        version = subprocess.run([args.nextflow, "-version"], cwd=destination, capture_output=True, text=True, check=True, timeout=30)
        receipt["nextflow_version"] = version.stdout.strip()
        first = [args.nextflow, "run", "main.nf", "--fail_code", "23", "-ansi-log", "false", "-with-trace", "trace-first.tsv"]
        receipt["first_command"] = first
        receipt["first_exit_code"] = execute(first, destination, "first", args.timeout)
        initial = trace(destination / "trace-first.tsv")
        receipt["first_tasks"] = initial
        assert receipt["first_exit_code"] != 0, "Task failure must make Nextflow fail"
        assert len(initial) == 3, "The failed run must not launch its dependent report or retry a task"
        failed = next(row for row in initial if row["name"] == "FAIL_BRANCH")
        assert failed["status"] == "FAILED" and failed["exit"] == "23"
        survivors = [row for row in initial if row["name"].startswith("SURVIVE_BRANCH (")]
        assert len(survivors) == 2 and all(row["status"] == "COMPLETED" and row["exit"] == "0" for row in survivors)
        failure_time = int((Path(failed["workdir"]) / "exit.time").read_text())
        before = {}
        for row in survivors:
            work = Path(row["workdir"])
            assert int((work / "started.time").read_text()) < failure_time < int((work / "finished.time").read_text()), "Independent task must actually be running when failure occurs"
            output = next(work.glob("*.txt"))
            before[row["name"]] = {"workdir": str(work), "sha256": sha(output), "output": output.name,
                                   "mtime_ns": output.stat().st_mtime_ns, "finished_time": (work / "finished.time").read_text()}
        receipt.update(successful_branches_before_resume=before, report_blocked_on_failure=True, running_independent_tasks_not_killed=True)
        save()
        resumed = [args.nextflow, "run", "main.nf", "--fail_code", "0", "-ansi-log", "false", "-with-trace", "trace-resume.tsv", "-resume"]
        receipt["resume_command"] = resumed
        receipt["resume_exit_code"] = execute(resumed, destination, "resume", args.timeout)
        final = trace(destination / "trace-resume.tsv")
        receipt["resume_tasks"] = final
        assert receipt["resume_exit_code"] == 0 and len(final) == 4
        for row in final:
            if row["name"] in before:
                old = before[row["name"]]
                assert row["status"] == "CACHED" and row["exit"] == "0" and row["workdir"] == old["workdir"]
                output = Path(row["workdir"]) / old["output"]
                assert sha(output) == old["sha256"] and output.stat().st_mtime_ns == old["mtime_ns"]
                assert (output.parent / "finished.time").read_text() == old["finished_time"]
            else:
                assert row["name"] in {"FAIL_BRANCH", "BLOCKED_REPORT"} and row["status"] == "COMPLETED" and row["exit"] == "0"
        receipt.update(status="passed", resume_cached_successful_branches=2, failed_branch_recomputed=1,
                       dependent_report_started_only_after_repair=True, automatic_retries=0,
                       successful_output_bytes_and_mtimes_preserved=True)
    except Exception as error:
        receipt.update(status="failed", error=f"{type(error).__name__}: {error}")
    finally:
        receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save()
    print(json.dumps({"status": receipt["status"], "receipt": str(destination / "acceptance.receipt.json")}))
    return 0 if receipt["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
