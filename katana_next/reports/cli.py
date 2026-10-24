"""Build the seven-section Katana report from completed results, read only."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import json
from pathlib import Path
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path,
                        help="Existing katana.project.json or completed-report source descriptor JSON.")
    parser.add_argument("--outdir", required=True, type=Path,
                        help="New independent report directory; completed reports are never overwritten.")
    parser.add_argument("--event-set", default="stable_significant",
                        choices=("stable_significant", "significant"),
                        help="Default RS-DAS; significant selects ordinary DAS. Original source statistics remain unchanged.")
    parser.add_argument("--catalog", "--asset-catalog", dest="catalog", type=Path,
                        help="Optional JSON analysis-artifact export. --asset-catalog is an alias for --catalog.")
    parser.add_argument("--from-report", type=Path,
                        help="Refresh presentation of an existing portable report; completed analyses are never rerun.")
    parser.add_argument("--compact", action="store_true", help="Losslessly consolidate portable report storage.")
    archive_options=parser.add_mutually_exclusive_group()
    archive_options.add_argument("--no-archive", action="store_true", help="Generate the complete report directory without a duplicate ZIP.")
    archive_options.add_argument("--archive", type=Path, help="Write the portable ZIP at this new path instead of REPORT.zip.")
    parser.add_argument("--result-store", type=Path, help="Optional separate immutable payload cache; report files remain independently portable.")
    parser.add_argument("--path-map", type=Path, help="JSON list of explicit source_prefix/target_prefix bindings for this execution machine.")
    parser.add_argument("--network-display-parameters", type=Path,
                        help="Optional JSON mapping from comparison IDs to their original network parameter files, for display-only edge export.")
    args = parser.parse_args(argv)
    if args.result_store and args.from_report:
        parser.error('--result-store applies to fresh completed-result builds; --from-report already reuses its payloads')
    if args.path_map and args.from_report:
        parser.error('--path-map applies to completed-result inputs; portable reports already use relative payload paths')
    if args.network_display_parameters and not args.from_report:
        parser.error('--network-display-parameters requires --from-report')
    try:
        from katana_next.reports import build_report
        mappings=None
        if args.path_map:
            mappings=json.loads(args.path_map.read_text('utf-8-sig'))
            if isinstance(mappings,dict): mappings=mappings['path_mappings']
            mappings=[{**row,'target_prefix':str((args.path_map.resolve().parent/row['target_prefix']).resolve())}
                      for row in mappings]
        # Rendering-library diagnostics belong on stderr; stdout is one JSON
        # result, suitable for Nextflow and other command-line callers.
        with redirect_stdout(sys.stderr):
            if args.from_report:
                from .refresh import refresh_report
                parameters = None
                if args.network_display_parameters:
                    parameters = json.loads(args.network_display_parameters.read_text('utf-8-sig'))
                    parameters = {cid: str((args.network_display_parameters.resolve().parent / path).resolve())
                                  if not Path(path).is_absolute() else path for cid, path in parameters.items()}
                manifest = refresh_report(args.from_report, args.outdir, project=args.project,
                                          compact=args.compact, network_parameters=parameters,create_archive=False)
                if not args.no_archive:
                    from .refresh import write_archive
                    write_archive(args.outdir,archive_path=args.archive)
            else:
                if not args.project: parser.error('--project or --from-report is required')
                manifest = build_report(args.project, args.outdir, catalog=args.catalog,
                                        event_set=args.event_set,compact=args.compact,
                                        create_archive=not args.no_archive,archive_path=args.archive,
                                        result_store=args.result_store,path_mappings=mappings)
        output = args.outdir.expanduser().resolve()
        result = {
            "schema": "katana.report_cli_result/v1", "status": "completed",
            "project_id": manifest.get("project_id"), "event_set": manifest.get("event_set", args.event_set),
            "report": str(output / "REPORT.html"),
            "archive":None if args.no_archive else str(args.archive.resolve() if args.archive else output/"REPORT.zip"),
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
