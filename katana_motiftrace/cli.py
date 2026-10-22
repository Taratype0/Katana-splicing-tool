"""Batch command-line entry to the same MotifTrace API and workflow services."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from . import MotifTrace, __version__, scientific_release
from .catalog import mode_catalog, resolve_mode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, help="Existing Katana project manifest; required except for modes")
    parser.add_argument("--version", action="version", version=f"MotifTrace integration {__version__}; science {scientific_release}")
    commands = parser.add_subparsers(dest="action", required=True)
    commands.add_parser("modes", help="List established modes without opening a project")
    listing = commands.add_parser("artifacts", help="List existing result artifacts; no analysis")
    listing.add_argument("--mode")
    listing.add_argument("--run-fingerprint")
    selection = commands.add_parser("select", help="Save exact event row IDs from an existing source artifact")
    selection.add_argument("source_artifact_id")
    selection.add_argument("--rows", required=True, type=Path, help="JSON array of source row IDs")
    selection.add_argument("--name", default="MotifTrace events")
    for action in ("plan", "run"):
        command = commands.add_parser(action)
        command.add_argument("mode")
        command.add_argument("--parameters", type=Path, help="JSON object; same scientific parameters as the interface")
        command.add_argument("--selection-id")
        command.add_argument("--contrast", action="append", dest="contrast_ids", help="Exact declared comparison ID; repeat only when the selected mode supports it")
        if action == "run":
            command.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    if args.action == "modes":
        print(json.dumps(mode_catalog(), indent=2))
        return 0
    if not args.project:
        parser.error("--project is required for this action")
    try:
        if getattr(args, "mode", None):
            resolve_mode(args.mode)
        parameters = json.loads(args.parameters.read_text(encoding="utf-8")) if getattr(args, "parameters", None) else {}
        if not isinstance(parameters, dict):
            raise ValueError("Parameters must be a JSON object")
        with MotifTrace(args.project) as motiftrace:
            if args.action == "artifacts":
                result = [a.model_dump(mode="json") for a in motiftrace.artifacts(args.mode, run_fingerprint=args.run_fingerprint)]
            elif args.action == "select":
                row_ids = json.loads(args.rows.read_text(encoding="utf-8"))
                if not isinstance(row_ids, list):
                    raise ValueError("--rows must contain a JSON array")
                result = motiftrace.select(args.source_artifact_id, row_ids, name=args.name).model_dump(mode="json")
            elif args.action == "plan":
                result = motiftrace.plan(args.mode, parameters, selection_id=args.selection_id, contrast_ids=args.contrast_ids)
            else:
                result = motiftrace.run(args.mode, parameters, selection_id=args.selection_id,
                    contrast_ids=args.contrast_ids, force=args.force).model_dump(mode="json")
            print(json.dumps(result, indent=2, ensure_ascii=True))
            return 2 if args.action == "plan" and not result["ready"] else 0
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, ImportError) as error:
        print(json.dumps({"error": str(error)}), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
