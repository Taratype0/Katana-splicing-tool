"""Local CLI, not a separate scientific implementation or public hosting service."""
import argparse
import json
import sys
from pathlib import Path

from katana_next.client import KatanaProject


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path, help="Existing v2 project manifest")
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("workflows")
    sub.add_parser("design")
    sub.add_parser("design-export")
    sub.add_parser("jobs")
    retry = sub.add_parser("retry")
    retry.add_argument("job_id")
    contrasts = sub.add_parser("design-contrasts")
    contrasts.add_argument("--factors", nargs="+", default=["condition"])
    contrasts.add_argument("--include-batch", action="store_true")
    contrasts.add_argument("--references", default="{}", help="JSON object of explicit factor reference levels")
    contrasts.add_argument("--pairing", choices=["auto", "paired", "unpaired"], default="auto")
    edit = sub.add_parser("design-edit")
    edit.add_argument("edit_file", type=Path, help="JSON with revision, samples and reason")
    edit.add_argument("--apply", action="store_true", help="Explicitly save after validation; default is preview only")
    artifacts = sub.add_parser("artifacts")
    artifacts.add_argument("--workflow")
    module = sub.add_parser("module")
    module.add_argument("artifact_id", help="Source artifact from the desired frozen discovery run")
    module.add_argument("module", help="Exact module label")
    module.add_argument("--limit", type=int, default=30)
    for action in ("plan", "run"):
        command = sub.add_parser(action)
        command.add_argument("workflow")
        command.add_argument("--parameters", type=Path, help="JSON object using the same parameters as the web UI")
        if action == "run":
            command.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        with KatanaProject(args.project) as project:
            parameters = json.loads(args.parameters.read_text(encoding="utf-8")) if getattr(args, "parameters", None) else {}
            if not isinstance(parameters, dict):
                raise ValueError("Parameters must be a JSON object")
            if args.action == "workflows":
                result = project.inventory()
            elif args.action == "design":
                result = project.design()
            elif args.action == "design-contrasts":
                result = project.contrasts(args.factors, include_batch=args.include_batch, pairing=args.pairing, reference_levels=json.loads(args.references))
            elif args.action == "jobs":
                result = [j.model_dump(mode="json") for j in project.jobs.list(project.project_id)]
            elif args.action == "retry":
                result = project.retry(args.job_id).model_dump(mode="json")
            elif args.action == "design-export":
                design = project.editable_design()
                result = {"revision": design["revision"], "samples": design["samples"], "reason": ""}
            elif args.action == "design-edit":
                result = project.revise_design(json.loads(args.edit_file.read_text(encoding="utf-8")), apply=args.apply)
            elif args.action == "artifacts":
                result = [a.model_dump(mode="json") for a in project.artifacts(args.workflow)]
            elif args.action == "module":
                result = project.module(args.artifact_id, args.module, limit=args.limit)
            elif args.action == "plan":
                result = project.plan(args.workflow, parameters)
            else:
                result = project.run(args.workflow, parameters, force=args.force).model_dump(mode="json")
            print(json.dumps(result, ensure_ascii=True, indent=2))
            if args.action == "plan" and not result["ready"]:
                raise SystemExit(2)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(json.dumps({"error": str(exc)}), file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
