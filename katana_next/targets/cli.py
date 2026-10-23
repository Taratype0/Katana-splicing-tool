"""Export and execute portable selected-gene/event evidence requests."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import json
from pathlib import Path
import sys


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    export = sub.add_parser("export", help="Export selected source rows and completed evidence.")
    export.add_argument("--project", required=True, type=Path, help="Existing Katana project JSON.")
    export.add_argument("--selections", required=True, type=Path,
                        help="JSON array (or {selections: [...]}) of source_artifact_id, row_ids, target_kind, filters.")
    export.add_argument("--outdir", required=True, type=Path, help="New independent bundle directory.")
    export.add_argument("--module", action="append", dest="modules", choices=("gene_card", "event_card", "sashimi"))
    export.add_argument("--contrast-id", action="append", dest="contrast_ids",
                        help="Limit to declared comparisons; repeat for multiple comparisons.")
    for name, help_text in (("inspect", "Verify the request and show its targets/resources."),
                            ("tasks", "List target IDs for independent scheduler tasks."),
                            ("plan", "Check evidence, bound resources and output readiness."),
                            ("run", "Build independent target reports, reusing verified products.")):
        command = sub.add_parser(name, help=help_text)
        command.add_argument("--request", required=True, type=Path, help="Request ZIP, directory or request.json.")
        if name in ("plan", "run"):
            command.add_argument("--bindings", "--resources", dest="bindings", type=Path,
                                 help="Resource-ID to execution-machine path bindings JSON.")
            command.add_argument("--outdir", type=Path, required=(name == "run"))
        if name == "run":
            command.add_argument("--workers", type=int, help="Local worker budget; default derives from the allocation.")
            command.add_argument("--target-id", action="append", dest="target_ids",
                                 help="Execute this target only; repeat for a subset. Use a separate outdir per scheduler task.")
        if name == "tasks":
            command.add_argument("--format", choices=("json", "tsv"), default="json")
    inspect = sub.add_parser("results", help="Verify a returned result directory or manifest.")
    inspect.add_argument("--outdir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        # Scientific-library diagnostics belong on stderr; stdout is machine-readable.
        with redirect_stdout(sys.stderr):
            if args.action == "export":
                from katana_next.services.project_service import ProjectService
                from .bundle import export_targets
                selected = json.loads(args.selections.read_text(encoding="utf-8-sig"))
                if isinstance(selected, dict):
                    selected = selected["selections"]
                if not isinstance(selected, list) or not selected:
                    raise ValueError("Provide at least one explicit source selection.")
                context = ProjectService().open_project(args.project)
                result = export_targets(context, selected, args.outdir,
                                        modules=args.modules, contrast_ids=args.contrast_ids)
            elif args.action in ("inspect", "tasks"):
                from .bundle import open_bundle
                with open_bundle(args.request) as bundle:
                    request = bundle.request
                    if args.action == "inspect":
                        result = {"schema": "katana.target_inspection/v1", "verified": True,
                                  "request_id": request["request_id"], "project": request["project"],
                                  "targets": request["targets"], "resources": request.get("resources", []),
                                  "modules": request.get("modules", []),
                                  "evidence_count": len(request.get("evidence", []))}
                    else:
                        result = {"schema": "katana.target_tasks/v1", "request_id": request["request_id"],
                                  "tasks": [{"target_id": row["target_id"], "target_kind": row["target_kind"],
                                             "gene_id": row.get("gene_id"), "event_id": row.get("event_id"),
                                             "comparison_ids": row.get("comparison_ids", [])}
                                            for row in request["targets"]]}
            elif args.action == "plan":
                from .runner import plan_targets
                result = plan_targets(args.request, resources=args.bindings, outdir=args.outdir)
            elif args.action == "run":
                from .runner import run_targets
                if args.workers is not None and args.workers < 1:
                    raise ValueError("--workers must be positive.")
                result = run_targets(args.request, args.outdir, resources=args.bindings,
                                     workers=args.workers, target_ids=args.target_ids)
            else:
                from .runner import inspect_results
                result = inspect_results(args.outdir)
        if args.action == "tasks" and args.format == "tsv":
            print("target_id\ttarget_kind\tgene_id\tevent_id\tcomparison_ids")
            for task in result["tasks"]:
                values = [task["target_id"], task["target_kind"], task.get("gene_id") or "",
                          task.get("event_id") or "", ",".join(task["comparison_ids"])]
                if any("\t" in str(value) or "\n" in str(value) for value in values):
                    raise ValueError("Task identity cannot be encoded safely as TSV; use --format json.")
                print("\t".join(map(str, values)))
        else:
            print(json.dumps(result, ensure_ascii=True, allow_nan=False, indent=2))
        if args.action == "plan" and result.get("ready") is False:
            return 2
        if args.action == "run" and result.get("execution_status", result.get("status")) in ("failed", "partial"):
            return 2
        return 0
    except (OSError, ValueError, KeyError, RuntimeError, ImportError, TypeError) as error:
        print(json.dumps({"schema": "katana.target_cli_result/v1", "status": "failed",
                          "error": str(error), "no_upstream_rerun": True}, ensure_ascii=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
