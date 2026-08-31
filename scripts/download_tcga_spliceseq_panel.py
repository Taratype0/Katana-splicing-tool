from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

try:
    from scripts.prepare_tcga_spliceseq import prepare_tcga_spliceseq
except ModuleNotFoundError:  # Support direct execution from the scripts directory.
    from prepare_tcga_spliceseq import prepare_tcga_spliceseq


DEFAULT_COHORTS = [
    "KIRC",
    "KIRP",
    "KICH",
    "LUAD",
    "LUSC",
    "COAD",
    "READ",
    "BRCA",
    "PRAD",
    "LIHC",
]
DEFAULT_SPLICE_TYPES = ["ES"]
SUPPORTED_SPLICE_TYPES = {"ES", "RI", "AA", "AD", "ME", "AP", "AT"}
DOWNLOAD_URL = "https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIDownload"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download_cohort(
    cohort: str,
    event_type: str,
    output: Path,
    *,
    force: bool,
) -> Path:
    cohort_dir = output / cohort
    cohort_dir.mkdir(parents=True, exist_ok=True)
    event_type = event_type.upper()
    legacy = cohort_dir / f"PSI_download_{cohort}.txt"
    raw_path = (
        legacy
        if event_type == "ES" and legacy.exists() and not force
        else cohort_dir / f"PSI_download_{cohort}_{event_type}.txt"
    )
    if raw_path.exists() and not force:
        return raw_path

    archive_path = cohort_dir / f"PSI_download_{cohort}_{event_type}.zip"
    payload = urllib.parse.urlencode(
        [
            ("tissue", cohort),
            ("genes", ""),
            ("samples", ""),
            ("splicetype", event_type),
            ("pctwithval", "75"),
            ("avgexppct", "50"),
            ("psirange", "0"),
            ("psistd", "0"),
            ("annotation", "true"),
            ("clinicalData", "true"),
        ]
    ).encode("utf-8")
    request = urllib.request.Request(DOWNLOAD_URL, data=payload, method="POST")
    with urllib.request.urlopen(request, timeout=300) as response:
        with archive_path.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    with zipfile.ZipFile(archive_path) as archive:
        members = [name for name in archive.namelist() if not name.endswith("/")]
        if len(members) != 1:
            raise RuntimeError(
                f"Expected one TCGA SpliceSeq file for {cohort}; found {members}."
            )
        with archive.open(members[0]) as source, raw_path.open("wb") as target:
            shutil.copyfileobj(source, target)
    return raw_path


def _prepared_rows(output: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for matrix in sorted(output.glob("*/prepared/*_PSI_tumor_only.tsv")):
        cohort = matrix.parent.parent.name.upper()
        prefix = f"{cohort}_"
        suffix = "_PSI_tumor_only.tsv"
        if not matrix.name.startswith(prefix) or not matrix.name.endswith(suffix):
            continue
        event_type = matrix.name[len(prefix) : -len(suffix)].upper()
        parameter_path = (
            matrix.parent / "parameters.json"
            if event_type == "ES"
            else matrix.parent / f"parameters_{event_type}.json"
        )
        if not parameter_path.exists():
            continue
        parameters = json.loads(parameter_path.read_text(encoding="utf-8"))
        legacy_raw = matrix.parent.parent / f"PSI_download_{cohort}.txt"
        typed_raw = matrix.parent.parent / f"PSI_download_{cohort}_{event_type}.txt"
        raw_path = typed_raw if typed_raw.exists() else legacy_raw
        rows.append(
            {
                "cohort": cohort,
                "event_type": event_type,
                "source_url": DOWNLOAD_URL,
                "raw_file": str(raw_path.resolve()) if raw_path.exists() else "",
                "raw_sha256": _sha256(raw_path) if raw_path.exists() else "",
                "matrix": str(matrix.resolve()),
                "n_events": int(parameters.get("n_events", 0)),
                "n_samples": int(parameters.get("n_tumor_samples", 0)),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Download and prepare a reproducible TCGA SpliceSeq exon-skipping "
            "cohort panel."
        )
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--cohort",
        action="append",
        default=[],
        help="TCGA cohort code. Repeat for multiple cohorts; defaults to the locked panel.",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--splice-type",
        action="append",
        default=[],
        help=(
            "TCGA SpliceSeq event type. Repeat for multiple types; defaults to ES. "
            "Supported: ES, RI, AA, AD, ME, AP, AT."
        ),
    )
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cohorts = [value.upper() for value in (args.cohort or DEFAULT_COHORTS)]
    splice_types = [
        value.upper() for value in (args.splice_type or DEFAULT_SPLICE_TYPES)
    ]
    unsupported = sorted(set(splice_types).difference(SUPPORTED_SPLICE_TYPES))
    if unsupported:
        raise ValueError(f"Unsupported TCGA SpliceSeq event types: {unsupported}")
    rows: list[dict[str, object]] = []
    for cohort in cohorts:
        for event_type in splice_types:
            print(f"{cohort}/{event_type}: downloading/preparing")
            raw_path = _download_cohort(
                cohort, event_type, output, force=args.force
            )
            prepared_dir = output / cohort / "prepared"
            parameters = prepare_tcga_spliceseq(
                raw_path,
                prepared_dir,
                cohort=cohort,
                event_type=event_type,
            )
            rows.append(
                {
                    "cohort": cohort,
                    "event_type": event_type,
                    "source_url": DOWNLOAD_URL,
                    "raw_file": str(raw_path),
                    "raw_sha256": _sha256(raw_path),
                    "matrix": str(
                        prepared_dir
                        / f"{cohort}_{event_type}_PSI_tumor_only.tsv"
                    ),
                    "n_events": parameters["n_events"],
                    "n_samples": parameters["n_tumor_samples"],
                }
            )

    combined = {
        (str(row["cohort"]), str(row["event_type"])): row
        for row in [*_prepared_rows(output), *rows]
    }
    manifest = {
        "panel_id": "TCGA_SpliceSeq_multi_event_cohort_panel_v2",
        "cohorts": [combined[key] for key in sorted(combined)],
        "download_filters": {
            "event_types": sorted(
                {str(row["event_type"]) for row in combined.values()}
            ),
            "pct_with_values": 75,
            "minimum_average_expression_percentage": 50,
            "minimum_psi_range": 0,
            "minimum_psi_standard_deviation": 0,
        },
        "independence_boundary": (
            "Each TCGA disease type is an independent patient cohort, but all cohorts "
            "share the TCGA/SpliceSeq processing ecosystem and are not ten independent "
            "studies."
        ),
    }
    (output / "panel_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
