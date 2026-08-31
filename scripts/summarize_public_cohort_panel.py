from __future__ import annotations

import argparse
from pathlib import Path

from katana_cosplice.cohort_panel_validation import (
    evaluate_cohort_panel,
    write_cohort_panel_validation,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize a locked multi-cohort CoSplice benchmark."
    )
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--panel-config", type=Path, required=True)
    parser.add_argument("--panel-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate_cohort_panel(
        args.benchmark,
        args.panel_config,
        args.panel_manifest,
    )
    files = write_cohort_panel_validation(result, args.output.resolve())
    print(result.claim_audit.to_string(index=False))
    print("\nGenerated files:")
    for path in files:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
