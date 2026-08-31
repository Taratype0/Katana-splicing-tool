from __future__ import annotations

import argparse
from pathlib import Path

from katana_cosplice.method_extensions_benchmark import (
    run_method_extension_benchmark,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run known-truth, real read-count and public-cohort ablations for "
            "the CoSplice method extensions."
        )
    )
    parser.add_argument("--cosplice-output", required=True, type=Path)
    parser.add_argument("--public-preservation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = run_method_extension_benchmark(
        cosplice_output=args.cosplice_output.resolve(),
        public_preservation_path=args.public_preservation.resolve(),
        output_dir=args.output.resolve(),
    )
    print(result.recommendation.to_string(index=False))
    print(f"\nOutputs: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
