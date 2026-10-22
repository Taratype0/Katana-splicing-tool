from __future__ import annotations

import argparse
import csv
from pathlib import Path

from katana_next.domain.models import AnnotationReference, InputAsset, ProjectManifest
from katana_next.infrastructure.manifest_store import ManifestStore


FIXTURE_ROWS = (
    ("SE|chr1|+|100|200|300|400", "SE", -0.31, 0.002, 84, True, "Afg2a"),
    ("RI|chr2|-|500|700", "RI", 0.18, 0.011, 42, True, "Setd2"),
    ("A3SS|chr3|+|900|930|1100", "A3SS", 0.07, 0.24, 33, False, "Srsf7"),
    ("A5SS|chr4|-|1200|1270|1450", "A5SS", -0.22, 0.018, 12, True, "Clk1"),
    ("MXE|chr5|+|1500|1600|1700|1800", "MXE", 0.11, 0.19, 8, False, "Mbnl1"),
    ("SE|chr6|+|2100|2200|2300|2400", "SE", 0.27, 0.004, 61, True, "Ptbp1"),
    ("RI|chr7|-|2500|2680", "RI", -0.05, 0.44, 27, False, "Ythdc1"),
    ("SE|chr8|+|3000|3100|3200|3300", "SE", 0.14, 0.031, 19, True, "Srsf3"),
)


def create_phase1_fixture(root: str | Path) -> Path:
    base = Path(root).expanduser().resolve()
    inputs = base / "inputs"
    outputs = base / "results"
    inputs.mkdir(parents=True, exist_ok=True)
    outputs.mkdir(parents=True, exist_ok=True)
    table = inputs / "event_profile.tsv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "event_id",
                "event_type",
                "effect",
                "fdr",
                "read_support",
                "significant",
                "gene_symbol",
            ]
        )
        writer.writerows(FIXTURE_ROWS)

    manifest = ProjectManifest(
        project_id="katana-next-phase1-demo",
        name="Katana Next platform fixture",
        input_root=inputs,
        output_root=outputs,
        annotation=AnnotationReference(genome_build="GRCm39"),
        assets=(
            InputAsset(
                asset_id="phase1-events",
                kind="phase1_event_table",
                path=table,
            ),
        ),
        metadata={
            "purpose": "Platform acceptance only; not a biological analysis.",
        },
    )
    return ManifestStore().save(manifest, base / "katana.project.json")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create the Katana Next Phase 1 fixture."
    )
    parser.add_argument(
        "root", help="External destination for fixture inputs and outputs."
    )
    args = parser.parse_args()
    print(create_phase1_fixture(args.root))


if __name__ == "__main__":
    main()
