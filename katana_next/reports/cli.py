"""Build the seven-section Katana report from completed results, read only."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import json
from pathlib import Path
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path,
                        help="Existing katana.project.json or completed-report source descriptor JSON.")
    parser.add_argument("--outdir", required=True, type=Path,
                        help="New independent report directory; completed reports are never overwritten.")
    parser.add_argument("--event-set", default="stable_significant",
                        choices=("stable_significant", "significant"),
                        help="Default RS-DAS; significant selects ordinary DAS. Original source statistics remain unchanged.")
    parser.add_argument("--catalog", "--asset-catalog", dest="catalog", type=Path,
                        help="Optional JSON analysis-artifact export. --asset-catalog is an alias for --catalog.")
    args = parser.parse_args(argv)
    try:
        from katana_next.reports import build_report
        # Rendering-library diagnostics belong on stderr; stdout is one JSON
        # result, suitable for Nextflow and other command-line callers.
        with redirect_stdout(sys.stderr):
            manifest = build_report(args.project, args.outdir, catalog=args.catalog,
                                    event_set=args.event_set)
        output = args.outdir.expanduser().resolve()
        result = {
            "schema": "katana.report_cli_result/v1", "status": "completed",
            "project_id": manifest.get("project_id"), "event_set": manifest.get("event_set", args.event_set),
            "report": str(output / "REPORT.html"), "archive": str(output / "REPORT.zip"),
            "manifest": str(output / "report-manifest.json"), "provenance": str(output / "provenance.json"),
            "figures": len(manifest.get("panels", [])), "complete_results": len(manifest.get("complete_results", [])),
            "unavailable": len(manifest.get("unavailable", [])),
            "no_upstream_rerun": True, "no_source_statistics_refit": True,
        }
        print(json.dumps(result, ensure_ascii=True, indent=2))
        return 0
    except (OSError, ValueError, KeyError, RuntimeError, ImportError) as error:
        print(json.dumps({"schema": "katana.report_cli_result/v1", "status": "failed",
                          "error": str(error), "no_upstream_rerun": True}, ensure_ascii=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
