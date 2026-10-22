from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import html
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas.errors import EmptyDataError


DIRECT_METHODS = {
    "cosplice_tight": "Katana CoSplice tight branch",
    "wgcna_exact": "Exact R WGCNA",
    "pearson_hierarchical": "Pearson hierarchical baseline",
    "graphical_lasso": "Graphical-lasso baseline",
}

METHOD_ALIASES = {
    "CoSplice": "cosplice_tight",
    "Katana CoSplice": "cosplice_tight",
    "Katana CoSplice tight branch": "cosplice_tight",
    "cosplice": "cosplice_tight",
    "katana_sparse": "cosplice_tight",
    "Exact WGCNA": "wgcna_exact",
    "WGCNA": "wgcna_exact",
}


@dataclass(frozen=True)
class PublicationDataset:
    dataset_id: str
    accession: str
    result_dir: Path
    data_type: str
    design: str
    validation_role: str
    citation_url: str
    notes: str = ""


@dataclass
class PublicationBenchmarkResult:
    datasets: pd.DataFrame
    method_metrics: pd.DataFrame
    paired_contrasts: pd.DataFrame
    hybrid_architecture: pd.DataFrame
    hybrid_incremental_programs: pd.DataFrame
    preservation: pd.DataFrame
    method_scope: pd.DataFrame
    claim_audit: pd.DataFrame
    files: list[Path]


def _resolve_path(value: str, manifest_path: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path.resolve()


def load_publication_manifest(path: Path) -> tuple[dict[str, Any], list[PublicationDataset]]:
    path = path.resolve()
    # Windows PowerShell commonly writes a UTF-8 BOM; accepting it keeps manifests portable.
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if int(payload.get("schema_version", 0)) != 1:
        raise ValueError("Publication benchmark manifest schema_version must be 1.")
    raw_datasets = payload.get("datasets")
    if not isinstance(raw_datasets, list) or not raw_datasets:
        raise ValueError("Publication benchmark manifest must define at least one dataset.")
    datasets: list[PublicationDataset] = []
    seen: set[str] = set()
    for raw in raw_datasets:
        dataset_id = str(raw.get("dataset_id", "")).strip()
        if not dataset_id:
            raise ValueError("Every publication benchmark dataset requires dataset_id.")
        if dataset_id in seen:
            raise ValueError(f"Duplicate publication benchmark dataset_id: {dataset_id}")
        seen.add(dataset_id)
        datasets.append(
            PublicationDataset(
                dataset_id=dataset_id,
                accession=str(raw.get("accession", dataset_id)),
                result_dir=_resolve_path(str(raw.get("result_dir", "")), path),
                data_type=str(raw.get("data_type", "event_by_sample_PSI")),
                design=str(raw.get("design", "not documented")),
                validation_role=str(raw.get("validation_role", "internal stability")),
                citation_url=str(raw.get("citation_url", "")),
                notes=str(raw.get("notes", "")),
            )
        )
    return payload, datasets


def _hash_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalise_method(value: object) -> str:
    text = str(value)
    return METHOD_ALIASES.get(text, text)


def _coerce_bool(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y"}:
            return True
        if normalized in {"false", "0", "no", "n", ""}:
            return False
    return bool(value)


def _read_method_metrics(dataset: PublicationDataset) -> tuple[pd.DataFrame, list[Path]]:
    candidates = [
        dataset.result_dir / "method_benchmark_summary.tsv",
        dataset.result_dir / "internal_matched_bootstrap_benchmark.tsv",
    ]
    path = next((candidate for candidate in candidates if candidate.exists()), None)
    if path is None:
        return pd.DataFrame(), []
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    if "method" not in frame:
        raise ValueError(f"Missing method column in {path}")
    frame["method"] = frame["method"].map(_normalise_method)
    frame.insert(0, "dataset_id", dataset.dataset_id)
    frame.insert(1, "accession", dataset.accession)
    if "valid_partition" not in frame:
        modules = pd.to_numeric(frame.get("n_modules"), errors="coerce")
        frame["valid_partition"] = modules.ge(2)
    else:
        frame["valid_partition"] = frame["valid_partition"].map(_coerce_bool)
    if "status" not in frame:
        frame["status"] = np.where(frame["valid_partition"], "finished", "degenerate")
    return frame, [path]


def _read_truth(dataset: PublicationDataset) -> tuple[pd.DataFrame, list[Path]]:
    path = dataset.result_dir / "known_truth_recovery.tsv"
    if not path.exists():
        return pd.DataFrame(), []
    frame = pd.read_csv(path, sep="\t")
    frame["method"] = frame["method"].map(_normalise_method)
    columns = ["method", "truth_adjusted_rand_index"]
    return frame[columns], [path]


def _read_preservation(dataset: PublicationDataset) -> tuple[pd.DataFrame, list[Path]]:
    candidates = [
        dataset.result_dir / "external_preservation_summary.tsv",
        dataset.result_dir / "reciprocal_external_preservation_summary.tsv",
        dataset.result_dir / "tumor_only_external_preservation_summary.tsv",
    ]
    path = next((candidate for candidate in candidates if candidate.exists()), None)
    if path is None:
        return pd.DataFrame(), []
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    if "method" in frame:
        frame["method"] = frame["method"].map(_normalise_method)
    frame.insert(0, "dataset_id", dataset.dataset_id)
    return frame, [path]


def _read_hybrid_audit(
    dataset: PublicationDataset,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    architecture_path = dataset.result_dir / "hybrid_architecture_summary.tsv"
    incremental_path = dataset.result_dir / "hybrid_tight_incremental_programs.tsv"
    try:
        architecture = (
            pd.read_csv(architecture_path, sep="\t", low_memory=False)
            if architecture_path.exists()
            else pd.DataFrame()
        )
    except EmptyDataError:
        architecture = pd.DataFrame()
    try:
        incremental = (
            pd.read_csv(incremental_path, sep="\t", low_memory=False)
            if incremental_path.exists()
            else pd.DataFrame()
        )
    except EmptyDataError:
        incremental = pd.DataFrame()
    files: list[Path] = []
    if not architecture.empty:
        architecture.insert(0, "dataset_id", dataset.dataset_id)
        files.append(architecture_path)
    if not incremental.empty:
        incremental.insert(0, "dataset_id", dataset.dataset_id)
        files.append(incremental_path)
    return architecture, incremental, files


def _method_scope_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "method": "cosplice_tight",
                "comparison_class": "primary tight branch",
                "input_contract": "event-by-sample PSI",
                "head_to_head_status": "direct",
                "scientific_role": (
                    "Sparse signed module discovery for compact local programs. This row is not "
                    "the complete hybrid architecture."
                ),
            },
            {
                "method": "cosplice_hybrid",
                "comparison_class": "production architecture",
                "input_contract": "shared event-by-sample PSI",
                "head_to_head_status": "multi-view audit",
                "scientific_role": (
                    "Parallel tight sparse and diffuse exact-WGCNA/TOM programs. It is evaluated "
                    "through branch metrics, incremental coverage, overlap and preservation rather "
                    "than an artificial single-partition ARI."
                ),
            },
            {
                "method": "wgcna_exact",
                "comparison_class": "published network comparator",
                "input_contract": "same event-by-sample PSI",
                "head_to_head_status": "direct",
                "scientific_role": "Exact R WGCNA::blockwiseModules on the identical event universe.",
            },
            {
                "method": "pearson_hierarchical",
                "comparison_class": "classical baseline",
                "input_contract": "same event-by-sample PSI",
                "head_to_head_status": "direct baseline",
                "scientific_role": "Fixed-cut Pearson hierarchical clustering baseline.",
            },
            {
                "method": "graphical_lasso",
                "comparison_class": "sparse inverse-covariance baseline",
                "input_contract": "same event-by-sample PSI",
                "head_to_head_status": "direct baseline",
                "scientific_role": "Shrinkage partial-correlation clustering baseline.",
            },
            {
                "method": "WGCNA/Mantel co-splicing",
                "comparison_class": "published domain comparator",
                "input_contract": "multi-exon gene expression profiles",
                "head_to_head_status": "separate reproduction track",
                "scientific_role": (
                    "Gene-level exon-profile similarity; not pooled with event-PSI module scores."
                ),
            },
            {
                "method": "CoSpliceNet",
                "comparison_class": "published domain comparator",
                "input_contract": "time-course splicing plus RBP/motif information",
                "head_to_head_status": "domain-specific only",
                "scientific_role": (
                    "RBP-centred regulatory network reference; not an interchangeable PSI "
                    "module partitioner."
                ),
            },
            {
                "method": "Spycone",
                "comparison_class": "published domain comparator",
                "input_contract": "time-course transcript/isoform usage",
                "head_to_head_status": "domain-specific only",
                "scientific_role": (
                    "Time-course isoform-cluster reference; requires a dedicated longitudinal track."
                ),
            },
            {
                "method": "MAJIQ HET",
                "comparison_class": "published event-robustness reference",
                "input_contract": "junction reads and local splicing variations",
                "head_to_head_status": "not a module competitor",
                "scientific_role": (
                    "Heterogeneity/event robustness reference, not a co-splicing module method."
                ),
            },
        ]
    )


def _paired_contrasts(metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    numeric_metrics = [
        ("mean_bootstrap_ARI", "higher"),
        ("min_bootstrap_ARI", "higher"),
        ("silhouette_cosine", "higher"),
        ("neutral_graph_modularity", "higher"),
        ("within_between_separation", "higher"),
        ("grey_fraction", "lower"),
        ("runtime_seconds", "lower"),
        ("runtime_seconds_initial_fit", "lower"),
        ("truth_adjusted_rand_index", "higher"),
    ]
    for dataset_id, frame in metrics.groupby("dataset_id", sort=False):
        indexed = frame.drop_duplicates("method").set_index("method")
        if not {"cosplice_tight", "wgcna_exact"}.issubset(indexed.index):
            continue
        cosplice_valid = bool(indexed.at["cosplice_tight", "valid_partition"])
        wgcna_valid = bool(indexed.at["wgcna_exact", "valid_partition"])
        for metric, direction in numeric_metrics:
            if metric not in indexed:
                continue
            left = pd.to_numeric(
                pd.Series([indexed.at["cosplice_tight", metric]]),
                errors="coerce",
            ).iloc[0]
            right = pd.to_numeric(
                pd.Series([indexed.at["wgcna_exact", metric]]), errors="coerce"
            ).iloc[0]
            interpretable = bool(
                pd.notna(left)
                and pd.notna(right)
                and (metric.startswith("runtime_") or (cosplice_valid and wgcna_valid))
            )
            if not interpretable:
                winner = "not_interpretable"
            elif np.isclose(float(left), float(right)):
                winner = "tie"
            elif direction == "higher":
                winner = "cosplice_tight" if left > right else "wgcna_exact"
            else:
                winner = "cosplice_tight" if left < right else "wgcna_exact"
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "metric": metric,
                    "preferred_direction": direction,
                    "cosplice_tight": left,
                    "wgcna_exact": right,
                    "cosplice_tight_minus_wgcna": (
                        float(left) - float(right) if pd.notna(left) and pd.notna(right) else np.nan
                    ),
                    "interpretable": interpretable,
                    "winner": winner,
                }
            )
    return pd.DataFrame(rows)


def _claim_audit(
    metrics: pd.DataFrame,
    preservation: pd.DataFrame,
    hybrid_architecture: pd.DataFrame,
) -> pd.DataFrame:
    direct = metrics.loc[
        metrics["method"].isin(["cosplice_tight", "wgcna_exact"])
    ].copy()
    valid = direct.loc[direct["valid_partition"].astype(bool)]
    paired_valid = 0
    for _, frame in valid.groupby("dataset_id"):
        if {"cosplice_tight", "wgcna_exact"}.issubset(set(frame["method"])):
            paired_valid += 1
    has_external = not preservation.empty
    has_hybrid = not hybrid_architecture.empty
    return pd.DataFrame(
        [
            {
                "claim": "The implementation executes exact R WGCNA.",
                "status": "supported"
                if "wgcna_exact" in set(metrics["method"])
                else "not evaluated",
                "evidence": "method metrics and method availability audit",
            },
            {
                "claim": "The CoSplice tight branch and WGCNA use the same event universe.",
                "status": "supported by benchmark contract",
                "evidence": "common_event_universe.tsv / discovery-only event selection",
            },
            {
                "claim": "At least one dataset has two valid, directly comparable partitions.",
                "status": "supported" if paired_valid else "not supported",
                "evidence": (
                    f"{paired_valid} dataset(s) with valid CoSplice-tight and WGCNA partitions"
                ),
            },
            {
                "claim": "The complete tight/diffuse hybrid architecture was audited.",
                "status": "supported" if has_hybrid else "not evaluated",
                "evidence": (
                    "hybrid_architecture_summary.tsv and "
                    "hybrid_tight_incremental_programs.tsv"
                ),
            },
            {
                "claim": "The diffuse branch is an independent new method superior to WGCNA.",
                "status": "rejected",
                "evidence": (
                    "The exact diffuse branch uses WGCNA labels by construction; its value is "
                    "preserving broad programs alongside the complementary tight branch."
                ),
            },
            {
                "claim": "Independent-cohort module preservation was evaluated.",
                "status": "supported" if has_external else "not evaluated",
                "evidence": "external preservation summary",
            },
            {
                "claim": "The complete hybrid CoSplice architecture is universally superior.",
                "status": "not supported",
                "evidence": (
                    "A universal claim is invalid when metrics trade off, partitions degenerate, "
                    "or tools have different input domains."
                ),
            },
            {
                "claim": "Domain-specific tools were scored as if they were interchangeable.",
                "status": "rejected by design",
                "evidence": "method_scope.tsv",
            },
        ]
    )


def _render_scorecard(
    metrics: pd.DataFrame,
    contrasts: pd.DataFrame,
    output: Path,
) -> None:
    direct = metrics.loc[
        metrics["method"].isin(["cosplice_tight", "wgcna_exact"])
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    status = (
        direct.assign(valid=lambda frame: frame["valid_partition"].astype(bool).map({True: 1, False: 0}))
        .pivot_table(index="dataset_id", columns="method", values="valid", aggfunc="first")
        .reindex(columns=["cosplice_tight", "wgcna_exact"])
    )
    if not status.empty:
        sns.heatmap(
            status,
            cmap=sns.color_palette(["#e5e7eb", "#0f766e"], as_cmap=True),
            vmin=0,
            vmax=1,
            annot=True,
            fmt=".0f",
            cbar=False,
            linewidths=1,
            linecolor="white",
            ax=axes[0],
        )
        axes[0].set_title("Valid non-degenerate partition", fontweight="bold")
        axes[0].set_xlabel("")
        axes[0].set_ylabel("")
    else:
        axes[0].axis("off")

    stability = contrasts.loc[
        contrasts["metric"].eq("mean_bootstrap_ARI") & contrasts["interpretable"].astype(bool)
    ].copy()
    if stability.empty:
        axes[1].axis("off")
        axes[1].text(
            0.5,
            0.5,
            "No dataset had two valid partitions\nwith interpretable bootstrap ARI.",
            ha="center",
            va="center",
            fontsize=12,
        )
    else:
        long = stability.melt(
            id_vars="dataset_id",
            value_vars=["cosplice_tight", "wgcna_exact"],
            var_name="method",
            value_name="mean_bootstrap_ARI",
        )
        sns.barplot(
            data=long,
            x="dataset_id",
            y="mean_bootstrap_ARI",
            hue="method",
            palette={"cosplice_tight": "#0f766e", "wgcna_exact": "#d97706"},
            ax=axes[1],
        )
        axes[1].set_ylim(0, 1)
        axes[1].set_title("Matched-bootstrap partition stability", fontweight="bold")
        axes[1].set_xlabel("")
        axes[1].tick_params(axis="x", rotation=25)
        axes[1].legend(frameon=False)
    fig.suptitle(
        "Publication head-to-head audit: CoSplice tight branch vs exact WGCNA",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _markdown_table(frame: pd.DataFrame, max_rows: int = 50) -> str:
    if frame.empty:
        return "No rows."
    display = frame.head(max_rows).copy()
    for column in display.select_dtypes(include=["float"]).columns:
        display[column] = display[column].map(
            lambda value: "" if pd.isna(value) else f"{float(value):.4g}"
        )
    headers = [str(column).replace("|", "\\|") for column in display.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in display.itertuples(index=False, name=None):
        lines.append(
            "| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |"
        )
    if len(frame) > max_rows:
        lines.append(f"\nPreview limited to {max_rows} of {len(frame)} rows.")
    return "\n".join(lines)


def _build_report(
    result: PublicationBenchmarkResult,
    output: Path,
    title: str,
) -> tuple[Path, Path]:
    direct = result.method_metrics.loc[
        result.method_metrics["method"].isin(["cosplice_tight", "wgcna_exact"])
    ]
    valid_counts = (
        direct.groupby("method")["valid_partition"].sum().astype(int).to_dict()
        if not direct.empty
        else {}
    )
    scientific = result.paired_contrasts.loc[
        result.paired_contrasts["interpretable"].astype(bool)
        & ~result.paired_contrasts["metric"].str.startswith("runtime_")
    ]
    wins = (
        scientific["winner"]
        .value_counts()
        .to_dict()
        if not scientific.empty
        else {}
    )
    public_valid = result.method_metrics.loc[
        result.method_metrics["valid_partition"].astype(bool)
        & result.method_metrics["method"].isin(["cosplice_tight", "wgcna_exact"])
        & ~result.method_metrics["accession"].eq("simulation")
    ]
    public_valid_pairs = [
        dataset_id
        for dataset_id, frame in public_valid.groupby("dataset_id")
        if {"cosplice_tight", "wgcna_exact"}.issubset(set(frame["method"]))
    ]
    public_detail = ""
    public_interpretation = ""
    if public_valid_pairs:
        lines: list[str] = []
        interpretation_lines: list[str] = []
        for dataset_id in public_valid_pairs:
            frame = public_valid.loc[public_valid["dataset_id"].eq(dataset_id)].set_index("method")
            lines.append(
                f"- `{dataset_id}`: CoSplice-tight modules="
                f"{int(frame.at['cosplice_tight', 'n_modules'])}, "
                f"mean bootstrap ARI="
                f"{float(frame.at['cosplice_tight', 'mean_bootstrap_ARI']):.3f}; "
                f"exact WGCNA modules={int(frame.at['wgcna_exact', 'n_modules'])}, "
                f"mean bootstrap ARI={float(frame.at['wgcna_exact', 'mean_bootstrap_ARI']):.3f}."
            )
            cosplice = frame.loc["cosplice_tight"]
            wgcna = frame.loc["wgcna_exact"]
            interpretation = (
                f"- `{dataset_id}` CoSplice-tight versus exact WGCNA: "
                f"bootstrap ARI {float(cosplice.get('mean_bootstrap_ARI', float('nan'))):.3f} "
                f"versus {float(wgcna.get('mean_bootstrap_ARI', float('nan'))):.3f}; "
                f"cosine silhouette {float(cosplice.get('silhouette_cosine', float('nan'))):.3f} "
                f"versus {float(wgcna.get('silhouette_cosine', float('nan'))):.3f}; "
                f"within/between separation "
                f"{float(cosplice.get('within_between_separation', float('nan'))):.3f} versus "
                f"{float(wgcna.get('within_between_separation', float('nan'))):.3f}; "
                f"grey fraction {float(cosplice.get('grey_fraction', float('nan'))):.3f} versus "
                f"{float(wgcna.get('grey_fraction', float('nan'))):.3f}."
            )
            all_dataset_methods = result.method_metrics.loc[
                result.method_metrics["dataset_id"].eq(dataset_id)
            ].set_index("method")
            if (
                "pearson_hierarchical" in all_dataset_methods.index
                and bool(all_dataset_methods.at["pearson_hierarchical", "valid_partition"])
            ):
                interpretation += (
                    " Fixed Pearson hierarchical bootstrap ARI="
                    f"{float(all_dataset_methods.at['pearson_hierarchical', 'mean_bootstrap_ARI']):.3f}."
                )
            interpretation_lines.append(interpretation)
        public_detail = "\n".join(lines)
        public_interpretation = "\n".join(interpretation_lines)
    else:
        public_detail = "- No public biological dataset produced two valid partitions."
        public_interpretation = (
            "No direct public biological winner can be interpreted until both primary methods "
            "return non-degenerate partitions on the same event universe."
        )
    report = f"""# {title}

## Executive conclusion

This report is a publication-readiness audit, not a marketing score. Single-partition metrics compare
the CoSplice **tight sparse branch** with exact R WGCNA only where both receive the same event-by-sample
PSI universe and return non-degenerate partitions. The complete CoSplice hybrid is a parallel
tight/diffuse architecture, so it is audited through branch quality, incremental coverage, overlap
and preservation rather than an invented fused-partition ARI.

- Valid CoSplice-tight partitions: {valid_counts.get('cosplice_tight', 0)}
- Valid exact-WGCNA partitions: {valid_counts.get('wgcna_exact', 0)}
- Interpretable non-runtime metric wins: CoSplice-tight {wins.get('cosplice_tight', 0)}, exact WGCNA {wins.get('wgcna_exact', 0)}, ties {wins.get('tie', 0)}
- Universal superiority claim: **not supported**

The scientifically defensible claim is narrower: CoSplice preserves exact-WGCNA diffuse programs
while adding a complementary tight sparse view, event-level PSI inputs, bootstrap confidence,
replicate-aware filtering and external preservation. The added tight programs become biological
evidence only when they are reproducible and externally preserved.

## Current public biological head-to-head

{public_detail}

{public_interpretation}

These measurements support dataset-specific coverage/stability/compactness trade-offs, not
universal superiority. Other public tracks include deliberate failure boundaries: some datasets
yielded an all-grey/one-module result for one or both methods and are therefore not used as winner
evidence.

## Public datasets

{_markdown_table(result.datasets)}

## Direct method metrics

{_markdown_table(result.method_metrics)}

## Paired CoSplice-tight-minus-WGCNA contrasts

{_markdown_table(result.paired_contrasts)}

## Hybrid tight/diffuse architecture

{_markdown_table(result.hybrid_architecture)}

## Tight-branch incremental programs

{_markdown_table(result.hybrid_incremental_programs)}

## Independent-cohort preservation

{_markdown_table(result.preservation)}

## Comparator scope

{_markdown_table(result.method_scope)}

## Claim audit

{_markdown_table(result.claim_audit)}

## Interpretation rules

1. A zero- or one-module result is degenerate; its apparent stability is not interpreted.
2. Event selection must use the discovery cohort only. External samples cannot tune event selection,
   module labels or network parameters.
3. Exact WGCNA means the official R package. A Python soft-TOM approximation is never labelled WGCNA.
4. WGCNA/Mantel, CoSpliceNet, Spycone and MAJIQ HET answer adjacent but non-identical questions.
   Their published-domain reproductions are reported separately from the same-PSI head-to-head.
5. External module preservation and biological enrichment support generalization; neither proves
   causal regulation.
"""
    markdown = output / "PUBLICATION_HEAD_TO_HEAD_REPORT.md"
    markdown.write_text(report, encoding="utf-8")
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1180px;margin:36px auto;padding:0 24px;color:#172033}}
h1,h2{{color:#0f4c5c}} body{{line-height:1.55}} code{{background:#eef2f4;padding:2px 5px}}
pre{{white-space:pre-wrap;overflow-wrap:anywhere;font:15px/1.55 Arial,sans-serif}}
</style></head><body><pre>{html.escape(report)}</pre></body></html>"""
    html_path = output / "PUBLICATION_HEAD_TO_HEAD_REPORT.html"
    html_path.write_text(document, encoding="utf-8")
    return markdown, html_path


def build_publication_benchmark(manifest_path: Path, output: Path) -> PublicationBenchmarkResult:
    payload, datasets = load_publication_manifest(manifest_path)
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    dataset_rows: list[dict[str, object]] = []
    metric_frames: list[pd.DataFrame] = []
    preservation_frames: list[pd.DataFrame] = []
    hybrid_architecture_frames: list[pd.DataFrame] = []
    hybrid_incremental_frames: list[pd.DataFrame] = []
    source_files: list[Path] = [manifest_path.resolve()]
    for dataset in datasets:
        if not dataset.result_dir.exists():
            raise FileNotFoundError(
                f"Benchmark result directory does not exist for {dataset.dataset_id}: "
                f"{dataset.result_dir}"
            )
        metrics, metric_files = _read_method_metrics(dataset)
        truth, truth_files = _read_truth(dataset)
        if not truth.empty and not metrics.empty:
            metrics = metrics.merge(truth, on="method", how="left")
        preservation, preservation_files = _read_preservation(dataset)
        hybrid_architecture, hybrid_incremental, hybrid_files = _read_hybrid_audit(
            dataset
        )
        metric_frames.append(metrics)
        preservation_frames.append(preservation)
        hybrid_architecture_frames.append(hybrid_architecture)
        hybrid_incremental_frames.append(hybrid_incremental)
        source_files.extend(
            metric_files + truth_files + preservation_files + hybrid_files
        )
        dataset_rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "accession": dataset.accession,
                "data_type": dataset.data_type,
                "design": dataset.design,
                "validation_role": dataset.validation_role,
                "citation_url": dataset.citation_url,
                "result_dir": str(dataset.result_dir),
                "n_methods_loaded": int(metrics["method"].nunique()) if not metrics.empty else 0,
                "has_external_preservation": not preservation.empty,
                "has_hybrid_architecture_audit": not hybrid_architecture.empty,
                "notes": dataset.notes,
            }
        )
    dataset_frame = pd.DataFrame(dataset_rows)
    method_metrics = (
        pd.concat([frame for frame in metric_frames if not frame.empty], ignore_index=True)
        if any(not frame.empty for frame in metric_frames)
        else pd.DataFrame()
    )
    preservation = (
        pd.concat([frame for frame in preservation_frames if not frame.empty], ignore_index=True)
        if any(not frame.empty for frame in preservation_frames)
        else pd.DataFrame()
    )
    hybrid_architecture = (
        pd.concat(
            [frame for frame in hybrid_architecture_frames if not frame.empty],
            ignore_index=True,
        )
        if any(not frame.empty for frame in hybrid_architecture_frames)
        else pd.DataFrame()
    )
    hybrid_incremental = (
        pd.concat(
            [frame for frame in hybrid_incremental_frames if not frame.empty],
            ignore_index=True,
        )
        if any(not frame.empty for frame in hybrid_incremental_frames)
        else pd.DataFrame()
    )
    if method_metrics.empty:
        raise ValueError("No method benchmark result table was found in the configured datasets.")
    paired = _paired_contrasts(method_metrics)
    method_scope = _method_scope_table()
    claim_audit = _claim_audit(
        method_metrics,
        preservation,
        hybrid_architecture,
    )
    result = PublicationBenchmarkResult(
        datasets=dataset_frame,
        method_metrics=method_metrics,
        paired_contrasts=paired,
        hybrid_architecture=hybrid_architecture,
        hybrid_incremental_programs=hybrid_incremental,
        preservation=preservation,
        method_scope=method_scope,
        claim_audit=claim_audit,
        files=[],
    )
    tables = {
        "publication_dataset_manifest.tsv": dataset_frame,
        "dataset_method_metrics.tsv": method_metrics,
        "paired_cosplice_tight_vs_wgcna.tsv": paired,
        "hybrid_architecture_summary.tsv": hybrid_architecture,
        "hybrid_tight_incremental_programs.tsv": hybrid_incremental,
        "external_preservation.tsv": preservation,
        "method_scope.tsv": method_scope,
        "publication_claim_audit.tsv": claim_audit,
    }
    for name, frame in tables.items():
        path = output / name
        frame.to_csv(path, sep="\t", index=False)
        result.files.append(path)
    plot = output / "publication_head_to_head_scorecard.png"
    _render_scorecard(method_metrics, paired, plot)
    result.files.append(plot)
    markdown, html_report = _build_report(
        result,
        output,
        str(payload.get("benchmark_title", "Katana CoSplice publication head-to-head benchmark")),
    )
    result.files.extend([markdown, html_report])
    fingerprints = [
        {
            "path": str(path),
            "sha256": _hash_file(path),
            "bytes": path.stat().st_size,
        }
        for path in sorted(set(source_files))
        if path.exists() and path.is_file()
    ]
    parameters = output / "parameters.json"
    parameters.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "manifest": str(manifest_path.resolve()),
                "selection_fairness_rule": (
                    "Identical discovery-selected event universe for direct method comparisons."
                ),
                "primary_endpoints": [
                    "matched-bootstrap partition stability",
                    "independent-cohort module preservation",
                    "known-truth module recovery where available",
                ],
                "generated_files": [str(path) for path in result.files],
                "source_fingerprints": fingerprints,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    result.files.append(parameters)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build a publication-grade, traceable CoSplice head-to-head report from public "
            "benchmark result directories."
        )
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = build_publication_benchmark(args.manifest, args.output)
    print(f"Publication benchmark finished: {args.output.resolve()}")
    for path in result.files:
        print(path)
    return 0
