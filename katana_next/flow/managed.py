"""Run the existing native workflow with verified storage handoffs.

The scientific DAG and report renderer are unchanged. Separate native sessions
let a completed upstream bundle release its work before downstream analysis.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
from uuid import uuid4

from .model import canonical_hash, file_sha256, read_json, write_json


def verify_upstream_bundle(bundle):
    root = Path(bundle).resolve(strict=True)
    manifest = read_json(root / "manifest.json")
    if manifest.get("schema") != "katana-pipeline-bundle/v1" or manifest.get("status") != "complete":
        raise ValueError("Upstream checkpoint is not a completed native bundle")
    available, seen = [], set()
    for item in manifest.get("artifacts", []):
        if not str(item.get("status", "")).startswith("available"):
            continue
        rel = PurePosixPath(item["relative_path"])
        if rel.is_absolute() or ".." in rel.parts or ":" in str(rel) or str(rel) in seen:
            raise ValueError("Unsafe or duplicate upstream checkpoint member")
        seen.add(str(rel))
        path = root / str(rel)
        if (path.is_symlink() or not path.resolve().is_relative_to(root) or not path.is_file()
                or path.stat().st_size != item.get("bytes") or file_sha256(path) != item.get("sha256")):
            raise ValueError("Upstream checkpoint failed member verification: " + str(rel))
        available.append(str(rel))
    if not available or not manifest.get("samples") or not manifest.get("comparisons"):
        raise ValueError("Upstream checkpoint has no complete study or native products")
    return {"status": "completed", "bundle": str(root), "member_count": len(available),
            "manifest_sha256": file_sha256(root / "manifest.json")}


def verify_report(report, collection):
    from .collect import validate_collection
    root = Path(report).resolve(strict=True)
    proof = read_json(root / "flow.acceptance.json")
    collected = validate_collection(collection)
    if (proof.get("schema") != "katana.flow.report_acceptance/v1" or proof.get("status") != "completed"
            or proof.get("plan_id") != collected["plan_id"]):
        raise ValueError("The fixed report did not pass complete calculation acceptance")
    verify_delivered_report(root, proof)
    return {"status": "completed", "plan_id": collected["plan_id"],
            "archive_sha256": proof["archive_sha256"]}


def verify_delivered_report(root, proof):
    root = Path(root).resolve(strict=True)
    for name in ("REPORT.html", "REPORT.zip", "report-manifest.json"):
        path = root / name
        if not path.is_file() or path.is_symlink() or not path.stat().st_size:
            raise ValueError("Missing accepted final report product: " + name)
    if file_sha256(root / "REPORT.zip") != proof.get("archive_sha256"):
        raise ValueError("Accepted report archive has changed")
    return {"status": "completed", "archive_sha256": proof["archive_sha256"]}


def _copy_verified_tree(source, destination):
    """Copy a new final tree without following links into work or user data."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if destination.exists():
        raise FileExistsError("Refusing to overwrite a final delivery: " + str(destination))
    members = list(source.rglob("*"))
    if any(path.is_symlink() for path in members):
        raise ValueError("Final delivery must not depend on symbolic links")
    pending = destination.with_name("." + destination.name + ".pending-" + uuid4().hex)
    shutil.copytree(source, pending)
    for path in members:
        if path.is_file() and file_sha256(path) != file_sha256(pending / path.relative_to(source)):
            raise ValueError("Final delivery copy failed content verification")
    pending.rename(destination)


def _delivery_record(out):
    rows = []
    for name in ("report", "upstream-support", "pipeline_info", "runner-logs"):
        root = out / name
        for path in sorted(root.rglob("*")) if root.exists() else []:
            if path.is_symlink() or not path.resolve().is_relative_to(out):
                raise ValueError("Final supporting products cannot depend on symbolic links")
            if path.is_file():
                rows.append({"path": path.relative_to(out).as_posix(), "bytes": path.stat().st_size,
                             "sha256": file_sha256(path)})
    data = {"schema": "katana.flow.delivery/v1", "status": "completed", "files": rows,
            "retained_project_receipt_sha256": file_sha256(out / "project" / "retention.receipt.json")}
    return {**data, "receipt_sha256": canonical_hash(data)}


def _verify_delivery(out):
    from .storage_project import validate_retained_project
    proof = read_json(out / "delivery.receipt.json")
    if (proof.get("schema") != "katana.flow.delivery/v1" or proof.get("status") != "completed" or
            proof.get("receipt_sha256") != canonical_hash({k: v for k, v in proof.items() if k != "receipt_sha256"})):
        raise ValueError("Final delivery has no valid completed receipt")
    if file_sha256(out / "project" / "retention.receipt.json") != proof["retained_project_receipt_sha256"]:
        raise ValueError("Final retained-project receipt changed")
    for row in proof["files"]:
        relative = Path(row["path"])
        path = out / relative
        if (relative.is_absolute() or ".." in relative.parts or path.is_symlink() or
                not path.resolve().is_relative_to(out) or not path.is_file() or
                path.stat().st_size != row["bytes"] or file_sha256(path) != row["sha256"]):
            raise ValueError("Final delivery member is missing or changed: " + row["path"])
    validate_retained_project(out / "project")
    return verify_report(out / "report", out / "project" / "collection")


def _software_identity(pipeline):
    package = Path(__file__).resolve().parent.parent
    source = {"package/" + p.relative_to(package).as_posix(): file_sha256(p)
              for p in sorted(package.rglob("*.py")) if "__pycache__" not in p.parts}
    native = [pipeline]
    if (pipeline.parent / "nextflow.config").is_file():
        native.append(pipeline.parent / "nextflow.config")
    # Visit only code directories, never a neighbouring study's work/results.
    for name in ("modules", "workflows", "subworkflows", "conf", "bin", "lib"):
        folder = pipeline.parent / name
        if folder.is_dir():
            native.extend(p for p in folder.rglob("*") if p.is_file() and
                          p.suffix in {".nf", ".config", ".py", ".R", ".awk"})
    source.update({"pipeline/" + p.relative_to(pipeline.parent).as_posix(): file_sha256(p) for p in sorted(native)})
    return canonical_hash(source)


def _check_tool_cache_locations(work):
    for name in ("conda", "singularity", "apptainer"):
        path = work / "tasks" / name
        if path.exists() and not path.is_symlink():
            raise ValueError("Tool cache is inside reclaimable work; configure it outside this run before cleanup: " + str(path))
    for path in work.rglob("*") if work.exists() else []:
        if not path.is_symlink() and (path.name == "conda-meta" or path.suffix.lower() in {".sif", ".img"}):
            raise ValueError("A tool environment/image is inside reclaimable work; retain it outside this run: " + str(path))


def _native_phase(command, work_root, scope, stage, *, resume=False):
    launch = work_root / "launch"
    launch.mkdir(parents=True, exist_ok=True)
    log_dir = scope / "logs"
    log_dir.mkdir(exist_ok=True)
    attempt = len(list(log_dir.glob(stage + "-*.native.json"))) + 1
    previous = sorted(log_dir.glob(stage + "-*.native.json"))
    if previous and resume and not read_json(previous[-1]).get("engine_stopped", False):
        raise RuntimeError("Previous engine termination is unverified; inspect native jobs before resuming. Work retained.")
    log = log_dir / f"{stage}-{attempt:03d}.log"
    # Never let a site-wide cleanup=true destroy a checkpoint before validation.
    overlay = launch / "retention.config"
    overlay.write_text("cleanup = false\n", encoding="utf-8")
    actual = [*command, "-c", str(overlay), "-work-dir", str(work_root / "tasks")]
    if resume:
        actual.append("-resume")
    started = datetime.now(timezone.utc).isoformat()
    process = None
    try:
        with log.open("ab") as handle:
            environment = os.environ.copy()
            shared_cache = scope.parent / ".katana-runtime-cache"
            for engine in ("CONDA", "SINGULARITY", "APPTAINER"):
                # Explicit site config has priority over these Nextflow default
                # environment settings. Existing user settings stay untouched.
                environment.setdefault("NXF_" + engine + "_CACHEDIR", str(shared_cache / engine.lower()))
                configured = Path(environment["NXF_" + engine + "_CACHEDIR"]).expanduser().resolve()
                if configured.is_relative_to(work_root):
                    raise ValueError("Configured tool cache must be outside the reclaimable stage work")
            process = subprocess.Popen(actual, cwd=launch, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            try:
                for line in iter(process.stdout.readline, b""):
                    handle.write(line)
                    handle.flush()
                    sys.stderr.write(line.decode("utf-8", errors="replace"))
                result_code = process.wait()
            finally:
                process.stdout.close()
    except BaseException as error:
        not_started = process is None and isinstance(error, (FileNotFoundError, ValueError))
        write_json(log.with_suffix(".native.json"), {"state": "NOT_STARTED" if not_started else "INTERRUPTED",
            "engine_stopped": not_started, "stage": stage,
            "command": actual, "started_at": started, "scope": "No cleanup authorized"})
        raise
    engine_log = launch / ".nextflow.log"
    orderly_end = (engine_log.is_file() and
                   "Execution complete -- Goodbye" in engine_log.read_text(encoding="utf-8", errors="replace"))
    receipt = {"engine": "nextflow", "state": "COMPLETED" if result_code == 0 else "FAILED",
               "exit_code": result_code, "engine_stopped": result_code == 0 or orderly_end,
               "tool_caches": {engine.lower(): environment["NXF_" + engine + "_CACHEDIR"]
                               for engine in ("CONDA", "SINGULARITY", "APPTAINER")},
               "stage": stage, "command": actual,
               "started_at": started, "ended_at": datetime.now(timezone.utc).isoformat(), "log": str(log)}
    write_json(log.with_suffix(".native.json"), receipt)
    # Preserve the small native diagnostics outside the reclaimable work root.
    diagnostic = log_dir / f"{stage}-{attempt:03d}.nextflow.log"
    if (launch / ".nextflow.log").is_file():
        shutil.copyfile(launch / ".nextflow.log", diagnostic)
    if result_code:
        raise RuntimeError(f"Native {stage} failed ({result_code}); work retained. See {log}")
    return receipt


def _pipeline_source(explicit=None):
    if explicit:
        candidate = Path(explicit).expanduser().resolve()
    else:
        package = Path(__file__).resolve().parent
        candidate = package / "_pipeline" / "main.nf"
        if not candidate.is_file():
            candidate = package.parents[1] / "pipeline" / "main.nf"
    if not candidate.is_file():
        raise FileNotFoundError("Installed Nextflow pipeline is missing")
    return candidate


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("complete", "downstream-only"), default="complete")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--run-dir", type=Path, help="New private run directory; defaults to a sibling of outdir")
    p.add_argument("--resume", action="store_true", help="Resume this exact managed run; completed handoffs are reused")
    p.add_argument("--keep-work", action="store_true", help="Retain work even after success for native -resume")
    p.add_argument("--pipeline", type=Path, help="Optional explicit main.nf; default uses the installed pipeline")
    p.add_argument("--nextflow-command-json", default='["nextflow"]', help="Native executable prefix as a JSON string list")
    p.add_argument("--profile", default="local")
    p.add_argument("--site-config", action="append", default=[], type=Path)
    for key in ("input", "contrasts", "fasta", "gtf", "project_template", "flow_config", "project", "flow_plan"):
        p.add_argument("--" + key, type=Path)
    p.add_argument("--katana_container")
    p.add_argument("--dexseq_tools_container")
    p.add_argument("--flow_executable", default="katana-flow")
    p.add_argument("--params-file", type=Path, help="Native Nextflow scientific/resource parameter file; managed routing overrides are reserved")
    return p


def run(args):
    from .storage import create_managed_scope, managed_run_lock, register_stage, cleanup_stage
    out = args.outdir.expanduser().resolve()
    args.outdir = out
    scope = (args.run_dir or out.with_name(out.name + ".run")).expanduser().resolve()
    prefix = json.loads(args.nextflow_command_json)
    if not isinstance(prefix, list) or not prefix or any(not isinstance(v, str) or not v for v in prefix):
        raise ValueError("--nextflow-command-json must be a nonempty list of executable arguments")
    needed = ["project"] if args.mode == "downstream-only" else ["input", "contrasts", "fasta", "gtf", "project_template"]
    if not args.flow_plan:
        needed.append("flow_config")
    if args.mode == "complete" and args.flow_plan:
        raise ValueError("--flow_plan is only valid with --mode downstream-only")
    for key in needed:
        if getattr(args, key) is None or not getattr(args, key).expanduser().is_file():
            raise ValueError("Missing required file: --" + key)
    path_fields = ("input", "contrasts", "fasta", "gtf", "project_template", "flow_config", "project", "flow_plan", "params_file")
    for key in path_fields:
        if getattr(args, key) is not None:
            setattr(args, key, getattr(args, key).expanduser().resolve(strict=True))
    args.site_config = [v.expanduser().resolve(strict=True) for v in args.site_config]
    pipeline = _pipeline_source(args.pipeline)
    request = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()
               if key not in {"resume", "run_dir", "site_config"}}
    request.update(pipeline=str(pipeline), site_config=[str(v) for v in args.site_config])
    small = [getattr(args, key) for key in path_fields if key not in {"fasta", "gtf"} and getattr(args, key)]
    request["configuration_sha256"] = {str(p): file_sha256(p) for p in [*small, *args.site_config]}
    request["reference_metadata"] = {str(p): [p.stat().st_size, p.stat().st_mtime_ns]
                                     for p in (args.fasta, args.gtf) if p}
    request["software_identity"] = _software_identity(pipeline)
    request_id = canonical_hash(request)
    state_file = scope / "managed-run.json"
    if args.resume:
        state = read_json(state_file)
        if state.get("request_id") != request_id:
            raise ValueError("Resume request differs; preserve this run and use a new --run-dir for new inputs/settings")
    else:
        if out.exists() and any(out.iterdir()):
            raise FileExistsError("Choose a new final output directory")
        create_managed_scope(scope, protected_roots=[out, pipeline, *[getattr(args, k) for k in path_fields if getattr(args, k)]])
        state = {"schema": "katana.flow.managed_run/v1", "request_id": request_id, "request": request,
                 "status": "running", "stages": {}}
        write_json(state_file, state)
    with managed_run_lock(scope):
        def finish():
            _verify_delivery(out)
            state.update(status="completed", final_project=str(out / "project" / "project.json"),
                         report=str(out / "report" / "REPORT.html"), archive=str(out / "report" / "REPORT.zip"),
                         native_resume_available=args.keep_work,
                         completed_at=datetime.now(timezone.utc).isoformat())
            write_json(state_file, state)
            write_json(out / "managed-run.receipt.json", state)
            return state
        if state.get("status") == "completed":
            _verify_delivery(out)
            return state
        if (scope / "cleanup-downstream.receipt.json").exists():
            # Recover the small terminal-commit window without re-creating work
            # or requiring the already deleted upstream checkpoint again.
            _verify_delivery(out)
            cleanup_stage(scope, "downstream", native_status=state["stages"]["downstream"]["native"],
                          checkpoint_validator=lambda _: _verify_delivery(out))
            return finish()
        register_stage(scope, "downstream", scope / "work" / "downstream")
        down_work = scope / "work" / "downstream"
        common = [*prefix, "run", str(pipeline), "-profile", args.profile, "-ansi-log", "false"]
        for config in args.site_config:
            common.extend(["-c", str(config)])
        if args.params_file:
            common.extend(["-params-file", str(args.params_file)])
        for key in ("katana_container", "dexseq_tools_container", "flow_executable"):
            if getattr(args, key):
                common.extend(["--" + key, getattr(args, key)])
        common += ["--storage_mode", "lean"]
        try:
            if args.mode == "complete":
                bundle = down_work / "upstream-output" / "katana_bundle"
                project = Path(state.get("working_project", down_work / "project" / "project.json"))
                if "upstream" not in state["stages"]:
                    register_stage(scope, "upstream", scope / "work" / "upstream")
                    up_work = scope / "work" / "upstream"
                    command = [*common, "--mode", "upstream-only", "--outdir", str(bundle.parent), "--managed_preflight", "true",
                               "--project_template", str(args.project_template), "--flow_config", str(args.flow_config)]
                    for key in ("input", "contrasts", "fasta", "gtf"):
                        command.extend(["--" + key, str(getattr(args, key))])
                    native = _native_phase(command, up_work, scope, "upstream", resume=args.resume)
                    verified = verify_upstream_bundle(bundle)
                    state["stages"]["upstream"] = {"native": native, "checkpoint": verified}
                    write_json(state_file, state)
                else:
                    verified = verify_upstream_bundle(bundle)
                    if verified != state["stages"]["upstream"]["checkpoint"]:
                        raise ValueError("The successful upstream handoff has changed")
                if not args.keep_work:
                    _check_tool_cache_locations(scope / "work" / "upstream")
                    cleanup_stage(scope, "upstream", native_status=state["stages"]["upstream"]["native"],
                                  checkpoint_validator=lambda _: verify_upstream_bundle(bundle))
                if not project.exists() or not (project.parent / "import.receipt.json").exists():
                    from .importer import import_bundle
                    if project.parent.exists():
                        if list((scope / "logs").glob("downstream-*.native.json")):
                            raise ValueError("The imported project changed after downstream execution; work retained")
                        project = down_work / ("project-import-" + uuid4().hex) / "project.json"
                    import_bundle(bundle, args.project_template, project.parent,
                                  owned_link_root=down_work)
                    state["working_project"] = str(project)
                    write_json(state_file, state)
            else:
                project = args.project
            phase_outputs = down_work / "outputs"
            if "downstream" not in state["stages"]:
                command = [*common, "--mode", "downstream-only", "--project", str(project), "--outdir", str(phase_outputs)]
                for key in ("flow_config", "flow_plan"):
                    if getattr(args, key):
                        command.extend(["--" + key, str(getattr(args, key))])
                native = _native_phase(command, down_work, scope, "downstream", resume=args.resume)
                proof = verify_report(phase_outputs / "report", phase_outputs / "collection")
                state["stages"]["downstream"] = {"native": native, "checkpoint": proof}
                write_json(state_file, state)
            verify_report(phase_outputs / "report", phase_outputs / "collection")
            from .storage_project import retain_collection_project, validate_retained_project
            out.mkdir(parents=True, exist_ok=True)
            if not (out / "project").exists():
                staging_project = out / (".project-pending-" + uuid4().hex)
                retain_collection_project(phase_outputs / "collection", project, staging_project, owned_scope=scope)
                staging_project.rename(out / "project")
            if args.mode == "complete" and not (out / "upstream-support").exists():
                support = down_work / "upstream-support"
                support.mkdir(exist_ok=True)
                for path in (down_work / "upstream-output").iterdir():
                    if path.name in {"katana_bundle", "reference_index"}:
                        continue
                    if path.is_symlink():
                        raise ValueError("Upstream supporting results must not be symbolic links")
                    if path.is_dir():
                        _copy_verified_tree(path, support / path.name)
                    else:
                        shutil.copyfile(path, support / path.name)
                        if file_sha256(path) != file_sha256(support / path.name):
                            raise ValueError("Upstream provenance copy changed")
                _copy_verified_tree(support, out / "upstream-support")
            if not (out / "report").exists():
                _copy_verified_tree(phase_outputs / "report", out / "report")
            if (phase_outputs / "pipeline_info").exists() and not (out / "pipeline_info").exists():
                _copy_verified_tree(phase_outputs / "pipeline_info", out / "pipeline_info")
            if not (out / "runner-logs").exists():
                _copy_verified_tree(scope / "logs", out / "runner-logs")
            if not (out / "delivery.receipt.json").exists():
                write_json(out / "delivery.receipt.json", _delivery_record(out))
            final_report = state["stages"]["downstream"]["checkpoint"]
            def final_checkpoint(_):
                _verify_delivery(out)
                return {"status": "completed", "project": str(out / "project"), **final_report}
            if not args.keep_work:
                _check_tool_cache_locations(down_work)
                cleanup_stage(scope, "downstream", native_status=state["stages"]["downstream"]["native"],
                              checkpoint_validator=final_checkpoint)
            else:
                final_checkpoint(down_work)
            # Reload after the source work has actually gone; do not just trust
            # that a renderer or pre-deletion copy appeared to succeed.
            return finish()
        except BaseException:
            state["status"] = "failed_or_interrupted"
            write_json(state_file, state)
            raise


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = run(args)
        print(json.dumps(result, indent=2, ensure_ascii=True))
        return 0
    except (ValueError, OSError, KeyError, RuntimeError, ImportError, TypeError) as exc:
        print(json.dumps({"status": "failed", "error": str(exc), "work_retained_on_failure": True}), file=sys.stderr)
        return 1


def clean_completed_run(scope, *, apply=False):
    """Preview or clean only a verified, already completed managed run."""
    from .storage import managed_run_lock, cleanup_stage
    scope = Path(scope).expanduser().resolve(strict=True)
    state_file = scope / "managed-run.json"
    state = read_json(state_file)
    if state.get("schema") != "katana.flow.managed_run/v1" or state.get("status") != "completed":
        raise ValueError("Only a completed managed run can be cleaned; failed/running/unknown work stays intact")
    out = Path(state["request"]["outdir"]).resolve(strict=True)
    with managed_run_lock(scope):
        checkpoint = _verify_delivery(out)
        results = []
        for stage in ("upstream", "downstream"):
            if stage not in state["stages"]:
                continue
            _check_tool_cache_locations(scope / "work" / stage)
            results.append(cleanup_stage(scope, stage, native_status=state["stages"][stage]["native"],
                checkpoint_validator=lambda _: {**_verify_delivery(out), "status": "completed"}, dry_run=not apply))
        if apply:
            _verify_delivery(out)
            state["native_resume_available"] = False
            write_json(state_file, state)
            write_json(out / "managed-run.receipt.json", state)
    return {"schema": "katana.flow.clean_managed_run/v1", "status": "completed" if apply else "planned",
            "dry_run": not apply, "run_dir": str(scope), "stages": results}
