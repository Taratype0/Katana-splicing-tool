from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from katana_cosplice.config import CoSpliceConfig
from katana_cosplice.utils import normalize_sample_id


@dataclass
class AuditResult:
    checks: pd.DataFrame
    samples: pd.DataFrame

    @property
    def passed(self) -> bool:
        blocking = self.checks.loc[self.checks["severity"].eq("error")]
        return bool(blocking.empty or blocking["passed"].all())


def resolve_rmats_root(path: Path) -> Path:
    candidates = [path, path / "rmats", path / "star_salmon" / "rmats"]
    for candidate in candidates:
        if (candidate / "bamlist").is_dir() and any(candidate.glob("*/rmats_post")):
            return candidate
    for candidate in path.glob("**/star_salmon/rmats"):
        if (candidate / "bamlist").is_dir() and any(candidate.glob("*/rmats_post")):
            return candidate
    raise FileNotFoundError(f"Could not locate an rMATS root with bamlist/ and comparison/rmats_post folders under {path}")


def load_sample_metadata(path: Path, selected_conditions: list[str] | None = None) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    frame = pd.read_csv(path, sep=sep, low_memory=False)
    rename: dict[str, str] = {}
    if "sample_id" not in frame.columns and "sample" in frame.columns:
        rename["sample"] = "sample_id"
    if "condition" not in frame.columns and "group" in frame.columns:
        rename["group"] = "condition"
    frame = frame.rename(columns=rename)
    required = {"sample_id", "condition"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Sample metadata is missing required columns: {', '.join(missing)}")
    frame["sample_id"] = frame["sample_id"].map(normalize_sample_id)
    frame["condition"] = frame["condition"].astype(str).str.strip()
    invalid_samples = frame["sample_id"].isin({"", "nan", "None"})
    invalid_conditions = frame["condition"].isin({"", "nan", "None"})
    if invalid_samples.any() or invalid_conditions.any():
        raise ValueError(
            "Sample metadata contains blank sample_id or condition values: "
            f"invalid_sample_rows={frame.index[invalid_samples].tolist()}, "
            f"invalid_condition_rows={frame.index[invalid_conditions].tolist()}"
        )
    duplicated = frame["sample_id"].duplicated(keep=False)
    if duplicated.any():
        values = sorted(frame.loc[duplicated, "sample_id"].unique())
        raise ValueError(f"Sample metadata contains duplicate normalized sample_id values: {values}")
    if selected_conditions:
        missing_conditions = sorted(set(selected_conditions).difference(frame["condition"].unique()))
        if missing_conditions:
            raise ValueError(f"Configured selected_conditions are absent from sample metadata: {missing_conditions}")
        frame = frame.loc[frame["condition"].isin(selected_conditions)].copy()
    return frame.reset_index(drop=True)


def _bamlist_samples(path: Path) -> list[str]:
    if not path.exists():
        return []
    raw = path.read_text(encoding="utf-8").replace("\r", "\n").replace("\n", ",")
    return [normalize_sample_id(value) for value in raw.split(",") if value.strip()]


def audit_inputs(config: CoSpliceConfig) -> AuditResult:
    rows: list[dict[str, object]] = []

    def add(name: str, passed: bool, severity: str, detail: str) -> None:
        rows.append({"check": name, "passed": bool(passed), "severity": severity, "detail": detail})

    add("sample_metadata_exists", config.sample_metadata.exists(), "error", str(config.sample_metadata))
    add("rmats_root_exists", config.rmats_root.exists(), "error", str(config.rmats_root))
    if not config.sample_metadata.exists() or not config.rmats_root.exists():
        return AuditResult(pd.DataFrame(rows), pd.DataFrame())

    try:
        samples = load_sample_metadata(config.sample_metadata, config.selected_conditions)
    except Exception as exc:
        add("sample_metadata_schema", False, "error", str(exc))
        return AuditResult(pd.DataFrame(rows), pd.DataFrame())
    add("sample_metadata_schema", True, "error", f"rows={len(samples)}, columns={len(samples.columns)}")
    add("unique_sample_ids", samples["sample_id"].is_unique, "error", f"n={len(samples)}")
    condition_counts = samples.groupby("condition", observed=True)["sample_id"].nunique()
    add("at_least_two_conditions", len(condition_counts) >= 2, "error", f"conditions={len(condition_counts)}")
    add(
        "replicates_per_condition",
        bool((condition_counts >= 3).all()),
        "warning",
        "; ".join(f"{condition}={count}" for condition, count in condition_counts.items()),
    )

    try:
        rmats_root = resolve_rmats_root(config.rmats_root)
    except FileNotFoundError as exc:
        add("rmats_layout", False, "error", str(exc))
        return AuditResult(pd.DataFrame(rows), samples)

    bamlist_dir = rmats_root / "bamlist"
    bamlist_samples: set[str] = set()
    for path in bamlist_dir.glob("*_bamlist.txt"):
        bamlist_samples.update(_bamlist_samples(path))
    metadata_samples = set(samples["sample_id"])
    add(
        "metadata_matches_rmats_bamlists",
        metadata_samples == bamlist_samples if not config.selected_conditions else metadata_samples.issubset(bamlist_samples),
        "error",
        f"metadata={len(metadata_samples)}, bamlist={len(bamlist_samples)}, missing={sorted(metadata_samples-bamlist_samples)}",
    )

    comparisons = [path for path in rmats_root.glob("*/rmats_post") if path.is_dir()]
    add("rmats_comparison_outputs", bool(comparisons), "error", f"comparison_dirs={len(comparisons)}")
    expected = [f"{event}.MATS.{config.rmats_mode}.txt" for event in ("SE", "RI", "A3SS", "A5SS", "MXE")]
    complete = sum(all((path / name).exists() for name in expected) for path in comparisons)
    add("complete_rmats_event_sets", complete == len(comparisons), "error", f"complete={complete}/{len(comparisons)}")

    if config.quant_root:
        quant_files = list(config.quant_root.glob("**/quant.sf")) if config.quant_root.exists() else []
        quant_ids = {
            sample_id
            for sample_id in metadata_samples
            if any(
                normalize_sample_id(parent.name) == sample_id
                for path in quant_files
                for parent in list(path.parents)[:4]
            )
        }
        add(
            "quant_sf_coverage",
            metadata_samples.issubset(quant_ids),
            "warning",
            f"quant={len(quant_ids)}, missing={sorted(metadata_samples-quant_ids)}",
        )
    else:
        add("quant_sf_coverage", False, "info", "quant_root not configured; transcript abundance and isoform support will be unavailable")

    candidate_tables = (
        list(config.support_root.glob("*/gene_level_integrated_candidates.tsv"))
        if config.support_root and config.support_root.exists()
        else []
    )
    add(
        "module_annotation_support",
        bool(candidate_tables),
        "info",
        f"support_root={config.support_root or 'not configured'}; candidate_tables={len(candidate_tables)}",
    )
    add(
        "custom_rbp_gene_set",
        bool(config.rbp_gene_set and config.rbp_gene_set.exists()),
        "info",
        str(config.rbp_gene_set or "not configured; built-in generic RBP seed remains available"),
    )

    for key, path in [("gtf", config.gtf), ("genome_fasta", config.genome_fasta), ("transcript_fasta", config.transcript_fasta)]:
        add(key, bool(path and path.exists()), "info", str(path or "not configured"))

    local_bams: list[Path] = []
    if config.bam_root and config.bam_root.exists():
        local_bams = list(config.bam_root.glob("**/*.bam"))
    local_bam_ids = {normalize_sample_id(path.name) for path in local_bams}
    missing_local_bams = sorted(metadata_samples.difference(local_bam_ids))
    bam_mode = config.matrix_input.input_mode == "bam_junctions"
    add(
        "local_bam_coverage",
        not missing_local_bams,
        "error" if bam_mode else "info",
        (
            f"input_mode={config.matrix_input.input_mode}; local_bams={len(local_bams)}; "
            f"missing={missing_local_bams}"
        ),
    )
    add(
        "bam_junction_gtf",
        bool(config.gtf and config.gtf.exists()),
        "error" if bam_mode else "info",
        str(config.gtf or "not configured"),
    )

    missing_covariates = [column for column in config.residual_covariates if column not in samples.columns]
    add("residual_covariates", not missing_covariates, "error", f"configured={config.residual_covariates}, missing={missing_covariates}")
    model_columns = list(dict.fromkeys([*config.association_factors, *config.factorial_factors]))
    if config.batch_column:
        model_columns.append(config.batch_column)
    missing_model_columns = sorted(set(model_columns).difference(samples.columns))
    add(
        "association_metadata_columns",
        not missing_model_columns,
        "error",
        f"configured={model_columns}, missing={missing_model_columns}",
    )
    return AuditResult(pd.DataFrame(rows), samples)
