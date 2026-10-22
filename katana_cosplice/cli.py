from __future__ import annotations

import argparse
import json

from katana_cosplice.config import load_config
from katana_cosplice.branding import MODULE_NAME
from katana_cosplice.discovery import discover_config
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.stage_registry import stage_specs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="katana-splice-network", description=f"{MODULE_NAME} sample-level co-splicing network analysis")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init", help="Discover project inputs and write a starter YAML config")
    init.add_argument("--project-root", required=True)
    init.add_argument("--output")
    init.add_argument("--config", required=True)

    stage = subparsers.add_parser("stage", help="Run exactly one registered stage")
    stage.add_argument("--config", required=True)
    stage.add_argument("--name", required=True, choices=[spec.key for spec in stage_specs(visible_only=True)])
    stage.add_argument("--force", action="store_true")

    legacy_commands = list(dict.fromkeys(spec.cli_command for spec in stage_specs() if spec.cli_command))
    for name in [*legacy_commands, "status", "run"]:
        command = subparsers.add_parser(name)
        command.add_argument("--config", required=True)
        if name in set(legacy_commands) | {
            "build-matrix",
            "qc",
            "fit-network",
            "consensus-modules",
            "associate",
            "validate-network",
            "annotate-modules",
            "benchmark-methods",
            "preserve-modules",
            "contrast-factors",
            "run",
        }:
            command.add_argument("--force", action="store_true")
        if name in {"fit-network", "associate"}:
            command.add_argument("--network", choices=["observed", "residual"], default="observed")
        if name == "run":
            command.add_argument("--no-residual", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "init":
        config = discover_config(args.project_root, args.output)
        path = config.write(args.config)
        print(path)
        return 0

    config = load_config(args.config)
    pipeline = CoSplicePipeline(config)
    if args.command == "stage":
        result = pipeline.run_stage(args.name, force=args.force)
        if args.name == "audit":
            print(result.to_string(index=False))
    elif args.command == "audit":
        table = pipeline.run_stage("audit")
        print(table.to_string(index=False))
    elif args.command == "build-matrix":
        pipeline.run_stage("matrix", force=args.force)
    elif args.command == "qc":
        pipeline.run_stage("stable_qc", force=args.force)
    elif args.command == "fit-network":
        pipeline.run_stage(f"network_{args.network}", force=args.force)
    elif args.command == "consensus-modules":
        pipeline.run_stage("consensus_modules", force=args.force)
    elif args.command == "associate":
        pipeline.run_stage(f"associate_{args.network}", force=args.force)
    elif args.command == "validate-network":
        pipeline.run_stage("validate_modules", force=args.force)
    elif args.command == "annotate-modules":
        pipeline.run_stage("annotate_modules", force=args.force)
    elif args.command == "benchmark-methods":
        pipeline.run_stage("benchmark_methods", force=args.force)
    elif args.command == "preserve-modules":
        pipeline.run_stage("preserve_modules", force=args.force)
    elif args.command == "report":
        print(pipeline.run_stage("report"))
    elif args.command == "status":
        print(json.dumps(pipeline.status(), indent=2))
    elif args.command == "run":
        pipeline.run_all(include_residual=not args.no_residual, force=args.force)
    else:
        spec = next((item for item in stage_specs() if item.cli_command == args.command), None)
        if spec is None:
            raise ValueError(f"No stage is registered for CLI command: {args.command}")
        pipeline.run_stage(spec.key, force=getattr(args, "force", False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
