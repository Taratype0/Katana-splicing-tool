"""One execution entrypoint for Nextflow tasks and the final report gate."""
from __future__ import annotations

from .container_runtime import sanitize_container_environment
sanitize_container_environment()

import argparse
from contextlib import redirect_stdout
import json
from pathlib import Path
import sys

from .model import load_plan, read_json, write_json


def template_preflight(template_path, config_path):
    """Check configured resources before upstream creates native result files."""
    from .planner import create_plan, preflight_plan
    # A template is a resource contract, not fake experimental data or fitted
    # results. The final imported project receives its own complete preflight.
    from katana_next.workflows.registry import WorkflowRegistry
    template, config = read_json(template_path), read_json(config_path)
    base = Path(template_path).resolve().parent
    assets = template.get("assets", [])
    errors = []
    required = {"gene_sets_gmt", "rbp_gene_sets", "rbp_motif_database", "fimo_runtime", "rbp_binding_peaks"}
    for kind in required - {a.get("kind") for a in assets}:
        errors.append(f"Missing fixed-report resource kind: {kind}")
    for asset in assets:
        path = Path(asset["path"]).expanduser()
        path = path if path.is_absolute() else base / path
        if not path.is_file():
            errors.append(f"Missing configured resource: {asset['asset_id']}")
    if not template.get("metadata", {}).get("species"):
        errors.append("Explicit project species is required")
    config_base = Path(config_path).resolve().parent
    from .stages import readiness
    for adapter, options in (("network.fit", config.get("network", {})),
                             (config.get("mapp", {}).get("adapter", "mapp.maei"), config.get("mapp", {}))):
        inputs = options.get("inputs", {})
        if adapter in {"mapp.prepare", "mapp.maei"}:
            from .bundled_runtime import mapp_input_defaults
            inputs = mapp_input_defaults(inputs)
        if adapter == "network.fit":
            # Raw-study network inputs will be generated from the verified
            # imported manifest rather than guessed future work-directory paths.
            inputs = {"config": options["config"]} if options.get("config") else {"project": str(Path(template_path).resolve())}
        path_keys = {"native_root", "inclusion_table", "sitecount_matrix", "design", "config"}
        inputs = {k: str((Path(v) if Path(v).is_absolute() else config_base / v).resolve())
                  if isinstance(v, str) and k in path_keys else v for k, v in inputs.items()}
        future = ["project"] if adapter == "network.fit" and not inputs.get("config") else []
        if adapter == "mapp.maei" and not inputs.get("inclusion_table"):
            future = ["inclusion_table", "sitecount_matrix", "design", "quantity_provenance"]
        issues = readiness({"adapter": adapter, "inputs": inputs, "parameters": options.get("parameters", {}),
                            "readiness_deferred_inputs": future})
        if isinstance(issues, dict):
            issues = issues.get("errors", issues.get("issues", []))
        errors.extend(str(issue) for issue in issues)
        if adapter == "mapp.maei" and future:
            issues = readiness({"adapter": "mapp.prepare", "inputs": inputs,
                                "parameters": options.get("prepare_parameters", {}),
                                "readiness_deferred_parameters": ["sample_ids", "site", "window_start", "window_end"],
                                "readiness_deferred_inputs": ["events", "gtf", "genome_fasta", "quantification"]})
            if isinstance(issues, dict):
                issues = issues.get("errors", issues.get("issues", []))
            errors.extend(str(issue) for issue in issues)
    qc = config.get("qc", {})
    bed = Path(qc.get("reference_bed", ""))
    bed = bed if bed.is_absolute() else config_base / bed
    if qc.get("reference_bed") and not bed.is_file():
        errors.append("Gene-body QC requires a reference BED matching the declared annotation")
    import shutil
    executable = qc.get("rseqc_executable", "geneBody_coverage.py")
    if not shutil.which(executable) and not Path(executable).is_file():
        errors.append("Missing RSeQC executable: " + executable)
    return {"schema": "katana.flow.template_preflight/v1", "ready": not errors, "errors": errors,
            "scope": "Resource/runtime readiness before upstream. Native project and calculation acceptance follow upstream export."}


def build_complete_report(collection, outdir):
    """No complete report unless every planned calculation has actual evidence."""
    from .collect import validate_collection, materialize_report_descriptor
    from .acceptance import required_report_gaps
    collection, outdir = Path(collection).resolve(), Path(outdir).resolve()
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise FileExistsError("Choose a new report output directory; existing report and acceptance evidence are preserved")
    receipt = validate_collection(collection)
    from katana_next.reports import build_report
    descriptor = collection / "report-source.json"
    if not descriptor.exists():
        descriptor = collection / "report-descriptor.json"
    if not descriptor.exists():
        descriptor = collection / "report-project.json"
    if not descriptor.exists():
        raise ValueError("Verified collection has no report descriptor")
    plan = load_plan(collection / "plan.json")
    if plan.get("report_profile") != "complete":
        raise ValueError("A test subgraph or partial collection cannot be published as a complete fixed report")
    # Rendering success and calculation success are distinct. Preserve a failed
    # render for diagnosis, and withhold its complete-report receipt/archive.
    outdir.mkdir(parents=True, exist_ok=True)
    required_gaps = []
    try:
        manifest = build_report(materialize_report_descriptor(collection), outdir, event_set=plan.get("event_set", "stable_significant"))
        required_gaps = required_report_gaps(collection, plan, manifest, outdir, collection_receipt=receipt)
        if required_gaps:
            raise ValueError("Required fixed-report analyses are still absent from the rendered result; see flow.acceptance.failed.json")
        chapters = {int(panel["chapter"]) for panel in manifest.get("panels", []) if panel.get("image")}
        if not set(range(1, 8)).issubset(chapters):
            raise ValueError("Fixed seven-chapter report is missing rendered figures; collection remains diagnostic")
        for name in ("REPORT.html", "REPORT.zip", "report-manifest.json"):
            if not (outdir / name).is_file() or not (outdir / name).stat().st_size:
                raise ValueError("Incomplete final report artifact: " + name)
        from .model import file_sha256
        import os
        result = {"schema": "katana.flow.report_acceptance/v1", "status": "completed", "plan_id": plan["plan_id"],
                  "verified_calculation_tasks": len(plan["tasks"]), "seven_chapter_specification": plan["report_specification"],
                  "report": str(outdir / "REPORT.html"), "archive": str(outdir / "REPORT.zip"),
                  "archive_sha256": file_sha256(outdir / "REPORT.zip")}
        pending = outdir / ".flow.acceptance.pending.json"
        with pending.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2, ensure_ascii=True, allow_nan=False) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        pending.rename(outdir / "flow.acceptance.json")
    except Exception as error:
        archive = outdir / "REPORT.zip"
        diagnostic = outdir / "REPORT.diagnostic.zip"
        if archive.exists():
            if diagnostic.exists():
                raise FileExistsError("Refusing to overwrite existing diagnostic archive") from error
            archive.rename(diagnostic)
        failure = {"schema": "katana.flow.report_acceptance/v1", "status": "failed",
                   "plan_id": plan["plan_id"], "required_gaps": required_gaps,
                   "error": str(error),
                   "message": "The renderer produced a diagnostic artifact; fixed-report calculation acceptance failed",
                   "diagnostic_archive": str(diagnostic) if diagnostic.is_file() else None}
        with (outdir / "flow.acceptance.failed.json").open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(failure, indent=2, ensure_ascii=True, allow_nan=False) + "\n")
        raise
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    command = sub.add_parser("export-pipeline")
    command.add_argument("--outdir", required=True, type=Path)
    command = sub.add_parser("seal-existing", help="Seal explicitly bounded historical results without running analysis")
    command.add_argument("--descriptor", required=True, type=Path)
    command.add_argument("--out", required=True, type=Path)
    command.add_argument("--allowed-root", required=True, action="append", type=Path)
    command.add_argument("--completion-evidence", action="append", default=[], type=Path)
    command.add_argument("--family-root", action="append", default=[], type=Path)
    command.add_argument("--include-file", action="append", default=[], type=Path)
    command.add_argument("--purpose", choices=("report-results", "analysis-inputs"), default="report-results")
    command.add_argument("--exclude-exports", action="store_true", help="Analysis-inputs only; retain original export pointers as provenance")
    command = sub.add_parser("adopt-existing", help="Copy a sealed historical snapshot and preserve its native evidence")
    command.add_argument("--source", required=True, type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    command = sub.add_parser("bind-existing", help="Verify an explicit compute overlay retains the original study/design")
    command.add_argument("--adoption", required=True, type=Path)
    command.add_argument("--project", required=True, type=Path)
    command.add_argument("--contrast-id", required=True, action="append")
    command.add_argument("--out", required=True, type=Path)
    for action in ("plan", "preflight"):
        command = sub.add_parser(action)
        command.add_argument("--project", type=Path)
        command.add_argument("--project-template", type=Path)
        command.add_argument("--config", required=True, type=Path)
        command.add_argument("--out", type=Path)
        if action == "plan":
            command.add_argument("--require-ready", action="store_true", help="Fail before scheduling when fixed-report resources are missing")
    command = sub.add_parser("import")
    command.add_argument("--bundle", required=True, type=Path)
    command.add_argument("--project-template", required=True, type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    command = sub.add_parser("validate-plan", help="Validate an existing sealed plan without replanning or changing its path")
    command.add_argument("--plan", required=True, type=Path)
    command.add_argument("--project", required=True, type=Path)
    command.add_argument("--out", required=True, type=Path)
    command = sub.add_parser("task")
    command.add_argument("--plan", required=True, type=Path)
    command.add_argument("--task-id", required=True)
    command.add_argument("--project", type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    command.add_argument("--dependency-receipts", nargs="*", default=[], type=Path)
    command.add_argument("--presentation-recovery", type=Path,
                         help="Explicit unchanged MotifTrace 64 MiB render failure receipt; reuse its full statistics in a new task directory")
    command = sub.add_parser("runtime-check", help="Probe the selected task runtime without starting analysis")
    command.add_argument("--plan", required=True, type=Path)
    command.add_argument("--task-id", required=True)
    command.add_argument("--project", required=True, type=Path)
    command.add_argument("--runtime-label")
    command.add_argument("--out", required=True, type=Path)
    command = sub.add_parser("collect")
    command.add_argument("--plan", required=True, type=Path)
    command.add_argument("--project", type=Path)
    command.add_argument("--task-results", required=True, nargs="+", type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    command = sub.add_parser("report")
    command.add_argument("--collection", required=True, type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            if args.action == "export-pipeline":
                from .distribution import export_pipeline
                result = {"pipeline": str(export_pipeline(args.outdir))}
            elif args.action == "seal-existing":
                from .adoption import seal_existing_results
                result = seal_existing_results(args.descriptor, args.out, allowed_roots=args.allowed_root,
                    completion_evidence=args.completion_evidence, family_roots=args.family_root,
                    include_files=args.include_file, purpose=args.purpose, include_exports=not args.exclude_exports)
            elif args.action == "adopt-existing":
                from .adoption import adopt_existing_results
                result = adopt_existing_results(args.source, args.outdir)
            elif args.action == "bind-existing":
                from .adoption import create_identity_projection
                result = create_identity_projection(args.adoption, args.project, args.out, contrast_ids=args.contrast_id)
            elif args.action in {"plan", "preflight"}:
                if args.out and args.out.exists():
                    raise ValueError("Refusing to overwrite an existing plan or readiness receipt")
                from .planner import create_plan, preflight_plan
                if args.action == "preflight" and args.project_template:
                    result = template_preflight(args.project_template, args.config)
                else:
                    if not args.project:
                        raise ValueError("--project is required for downstream planning")
                    plan = create_plan(args.project, args.config)
                    result = plan if args.action == "plan" else preflight_plan(plan)
                    if args.action == "plan" and args.require_ready:
                        checked = preflight_plan(plan)
                        if args.out:
                            write_json(args.out.with_name(args.out.stem + ".preflight.json"), checked)
                        if not checked["ready"]:
                            raise ValueError("Fixed-report readiness failed: " + json.dumps(checked["errors"], ensure_ascii=True))
                if args.out:
                    if args.out.exists():
                        raise ValueError("Refusing to overwrite an existing plan or readiness receipt")
                    write_json(args.out, result)
            elif args.action == "import":
                from .importer import import_bundle
                result = {"project": str(import_bundle(args.bundle, args.project_template, args.outdir))}
            elif args.action == "validate-plan":
                from .plan_reuse import validate_existing_plan
                if args.out.exists():
                    raise ValueError("Refusing to overwrite an existing plan validation receipt")
                result = validate_existing_plan(args.plan, args.project)
                write_json(args.out, result)
            elif args.action == "task":
                from .task import execute_task
                options = {"project_path": args.project} if args.project else {}
                if args.presentation_recovery:
                    options["presentation_recovery"] = args.presentation_recovery
                result = execute_task(args.plan, args.task_id, args.outdir, dependency_receipts=args.dependency_receipts, **options)
            elif args.action == "runtime-check":
                from .runtime_checks import check_task_runtime
                if args.out.exists():
                    raise ValueError("Refusing to overwrite an existing runtime probe receipt")
                result = check_task_runtime(args.plan, args.task_id, project_path=args.project, runtime_label=args.runtime_label)
                write_json(args.out, result)
            elif args.action == "collect":
                from .collect import collect_results
                options = {"project_path": args.project} if args.project else {}
                result = collect_results(args.plan, args.task_results, args.outdir, **options)
            else:
                result = build_complete_report(args.collection, args.outdir)
        if isinstance(result, Path):
            result = {"receipt": str(result)}
        print(json.dumps(result, ensure_ascii=True, indent=2, default=str))
        return 2 if (args.action in {"preflight", "validate-plan"} and not result["ready"]) or (args.action == "runtime-check" and result["status"] != "ready") else 0
    except (ValueError, OSError, KeyError, RuntimeError, ImportError, TypeError) as exc:
        print(json.dumps({"schema": "katana.flow.error/v1", "status": "failed", "error": str(exc)}), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
