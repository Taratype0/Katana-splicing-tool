from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom

from katana_cosplice.utils import bh_adjust, normalize_gene_id


BUILTIN_RBP_SYMBOLS: tuple[str, ...] = (
    "Srsf1", "Srsf2", "Srsf3", "Srsf4", "Srsf5", "Srsf6", "Srsf7", "Srsf9", "Srsf10",
    "Hnrnpa1", "Hnrnpa2b1", "Hnrnpc", "Hnrnpf", "Hnrnph1", "Hnrnpk", "Hnrnpm", "Hnrnpu",
    "Rbfox1", "Rbfox2", "Rbfox3", "Mbnl1", "Mbnl2", "Mbnl3", "Celf1", "Celf2", "Ptbp1", "Ptbp2",
    "Qk", "Nova1", "Nova2", "Tardbp", "Fus", "Sfpq", "Tra2a", "Tra2b", "Elavl1", "Elavl2", "Elavl3", "Elavl4",
)


@dataclass
class ModuleAnnotationResult:
    events: pd.DataFrame
    genes: pd.DataFrame
    module_summary: pd.DataFrame
    rbp_enrichment: pd.DataFrame
    support_sources: pd.DataFrame


def _bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(False, index=frame.index, dtype=bool)
    values = frame[column]
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    numeric = pd.to_numeric(values, errors="coerce")
    text = values.astype(str).str.strip().str.casefold()
    return numeric.fillna(0).ne(0) | text.isin({"true", "yes", "y", "significant"})


def _gene_key(gene_id: object, gene_symbol: object) -> str:
    normalized = normalize_gene_id(gene_id)
    if normalized:
        return f"id:{normalized}"
    symbol = str(gene_symbol or "").strip()
    return f"sym:{symbol.casefold()}" if symbol else ""


def _read_candidate_support(support_root: Path | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "gene_key", "gene_id", "gene_symbol", "n_candidate_comparisons", "DEG_significant",
        "rMATS_significant", "DEXSeq_significant", "DTU_significant", "best_DE_FDR",
        "best_rMATS_FDR", "best_DEXSeq_qvalue", "best_DTU_qvalue",
    ]
    if support_root is None or not support_root.exists():
        return pd.DataFrame(columns=columns), pd.DataFrame(
            [{"source_type": "candidate_support", "path": str(support_root or ""), "status": "not_configured", "rows": 0}]
        )
    paths = sorted(support_root.glob("*/gene_level_integrated_candidates.tsv"))
    if not paths:
        return pd.DataFrame(columns=columns), pd.DataFrame(
            [{"source_type": "candidate_support", "path": str(support_root), "status": "missing_files", "rows": 0}]
        )
    frames: list[pd.DataFrame] = []
    source_rows: list[dict[str, object]] = []
    wanted_columns = {
        "comparison_id", "gene_id", "gene_symbol", "DEG_significant", "rMATS_significant",
        "DEXSeq_significant", "DTU_significant", "DE_FDR", "best_DE_FDR", "best_rMATS_FDR",
        "rMATS_FDR", "best_DEXSeq_qvalue", "DEXSeq_qvalue", "best_DTU_qvalue", "DTU_qvalue",
    }
    for path in paths:
        try:
            available = set(pd.read_csv(path, sep="\t", nrows=0).columns)
            frame = pd.read_csv(path, sep="\t", usecols=sorted(available & wanted_columns), low_memory=False)
        except Exception as exc:
            source_rows.append({"source_type": "candidate_support", "path": str(path), "status": f"read_failed: {exc}", "rows": 0})
            continue
        frame["source_comparison"] = frame.get("comparison_id", path.parent.name)
        frame["gene_id"] = frame.get("gene_id", pd.Series("", index=frame.index)).map(normalize_gene_id)
        frame["gene_symbol"] = frame.get("gene_symbol", pd.Series("", index=frame.index)).fillna("").astype(str).str.strip()
        frame["gene_key"] = [
            _gene_key(gene_id, symbol) for gene_id, symbol in zip(frame["gene_id"], frame["gene_symbol"], strict=False)
        ]
        frame = frame.loc[frame["gene_key"].ne("")].copy()
        for flag in ("DEG_significant", "rMATS_significant", "DEXSeq_significant", "DTU_significant"):
            frame[flag] = _bool_series(frame, flag)
        frames.append(frame)
        source_rows.append({"source_type": "candidate_support", "path": str(path), "status": "loaded", "rows": len(frame)})
    if not frames:
        return pd.DataFrame(columns=columns), pd.DataFrame(source_rows)
    combined = pd.concat(frames, ignore_index=True, sort=False)
    q_columns = {
        "best_DE_FDR": ("DE_FDR", "best_DE_FDR"),
        "best_rMATS_FDR": ("best_rMATS_FDR", "rMATS_FDR"),
        "best_DEXSeq_qvalue": ("best_DEXSeq_qvalue", "DEXSeq_qvalue"),
        "best_DTU_qvalue": ("best_DTU_qvalue", "DTU_qvalue"),
    }
    for target, candidates in q_columns.items():
        source = next((name for name in candidates if name in combined), None)
        combined[target] = pd.to_numeric(combined[source], errors="coerce") if source else np.nan
    grouped = combined.groupby("gene_key", sort=False, observed=True).agg(
        gene_id=("gene_id", "first"),
        gene_symbol=("gene_symbol", "first"),
        n_candidate_comparisons=("source_comparison", "nunique"),
        DEG_significant=("DEG_significant", "max"),
        rMATS_significant=("rMATS_significant", "max"),
        DEXSeq_significant=("DEXSeq_significant", "max"),
        DTU_significant=("DTU_significant", "max"),
        best_DE_FDR=("best_DE_FDR", "min"),
        best_rMATS_FDR=("best_rMATS_FDR", "min"),
        best_DEXSeq_qvalue=("best_DEXSeq_qvalue", "min"),
        best_DTU_qvalue=("best_DTU_qvalue", "min"),
    ).reset_index()
    return grouped.reindex(columns=columns), pd.DataFrame(source_rows)


def _load_rbp_sets(custom_path: Path | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = [
        {"gene_set": "katana_builtin_splicing_rbp_seed", "source": "Katana built-in generic seed", "gene_symbol": symbol, "gene_id": ""}
        for symbol in BUILTIN_RBP_SYMBOLS
    ]
    source_rows = [{"source_type": "rbp_gene_set", "path": "built-in", "status": "loaded", "rows": len(rows)}]
    if custom_path is not None:
        if not custom_path.exists():
            source_rows.append({"source_type": "rbp_gene_set", "path": str(custom_path), "status": "missing_file", "rows": 0})
        else:
            try:
                custom = pd.read_csv(custom_path, sep="," if custom_path.suffix.lower() == ".csv" else "\t", low_memory=False)
                symbol_col = next((name for name in ("gene_symbol", "gene", "symbol", "Gene", "GeneSymbol") if name in custom), None)
                id_col = next((name for name in ("gene_id", "GeneID", "ensembl_gene_id") if name in custom), None)
                set_col = next((name for name in ("gene_set", "set", "category", "RBP", "rbp_gene_set") if name in custom), None)
                source_col = next((name for name in ("source", "database") if name in custom), None)
                if symbol_col is None and id_col is None:
                    raise ValueError("custom RBP table needs gene_symbol or gene_id")
                for _, item in custom.iterrows():
                    rows.append(
                        {
                            "gene_set": str(item.get(set_col, "custom_rbp_gene_set") or "custom_rbp_gene_set"),
                            "source": str(item.get(source_col, custom_path) or custom_path),
                            "gene_symbol": str(item.get(symbol_col, "") or "") if symbol_col else "",
                            "gene_id": normalize_gene_id(item.get(id_col)) if id_col else "",
                        }
                    )
                source_rows.append({"source_type": "rbp_gene_set", "path": str(custom_path), "status": "loaded", "rows": len(custom)})
            except Exception as exc:
                source_rows.append({"source_type": "rbp_gene_set", "path": str(custom_path), "status": f"read_failed: {exc}", "rows": 0})
    frame = pd.DataFrame(rows)
    frame["gene_key"] = [_gene_key(gene_id, symbol) for gene_id, symbol in zip(frame["gene_id"], frame["gene_symbol"], strict=False)]
    frame = frame.loc[frame["gene_key"].ne("")].drop_duplicates(["gene_set", "source", "gene_key"])
    return frame.reset_index(drop=True), pd.DataFrame(source_rows)


def _enrichment_pvalue(universe_size: int, selected_size: int, module_size: int, overlap: int) -> float:
    if min(universe_size, selected_size, module_size, overlap) <= 0:
        return float("nan")
    return float(hypergeom.sf(overlap - 1, universe_size, selected_size, module_size))


def annotate_modules(
    assignments: pd.DataFrame,
    event_metadata: pd.DataFrame,
    *,
    support_root: Path | None = None,
    rbp_gene_set: Path | None = None,
) -> ModuleAnnotationResult:
    required_assignment = {"event_key", "module"}
    required_event = {"event_key", "gene_id", "gene_symbol", "event_type"}
    if not required_assignment.issubset(assignments.columns):
        raise ValueError(f"Module assignments missing columns: {sorted(required_assignment - set(assignments.columns))}")
    if not required_event.issubset(event_metadata.columns):
        raise ValueError(f"Event metadata missing columns: {sorted(required_event - set(event_metadata.columns))}")
    events = assignments.merge(event_metadata, on="event_key", how="left", validate="one_to_one")
    events = events.loc[events["module"].astype(str).ne("grey")].copy()
    events["gene_id"] = events["gene_id"].map(normalize_gene_id)
    events["gene_symbol"] = events["gene_symbol"].fillna("").astype(str).str.strip()
    events["gene_key"] = [_gene_key(gene_id, symbol) for gene_id, symbol in zip(events["gene_id"], events["gene_symbol"], strict=False)]
    events["gene_key"] = events["gene_key"].astype("string")
    events = events.loc[events["gene_key"].ne("")].copy()
    support, support_sources = _read_candidate_support(support_root)
    rbp_sets, rbp_sources = _load_rbp_sets(rbp_gene_set)

    genes = (
        events.groupby(["module", "gene_key"], as_index=False)
        .agg(
            gene_id=("gene_id", "first"),
            gene_symbol=("gene_symbol", "first"),
            n_module_events=("event_key", "nunique"),
            event_types=("event_type", lambda values: ";".join(sorted(set(map(str, values))))),
        )
    )
    genes["gene_key"] = genes["gene_key"].astype("string")
    support["gene_key"] = support["gene_key"].astype("string")
    genes = genes.merge(support, on=["gene_key"], how="left", suffixes=("", "_support"))
    for name in ("gene_id", "gene_symbol"):
        fallback = f"{name}_support"
        if fallback in genes:
            genes[name] = genes[name].where(genes[name].astype(str).ne(""), genes[fallback].fillna(""))
            genes = genes.drop(columns=fallback)
    for flag in ("DEG_significant", "rMATS_significant", "DEXSeq_significant", "DTU_significant"):
        genes[flag] = genes[flag].map(lambda value: bool(value) if pd.notna(value) else False).astype(bool)
    rbp_ids = {normalize_gene_id(value) for value in rbp_sets["gene_id"] if normalize_gene_id(value)}
    rbp_symbols = {str(value).strip().casefold() for value in rbp_sets["gene_symbol"] if str(value).strip()}
    genes["RBP_gene"] = genes["gene_id"].map(normalize_gene_id).isin(rbp_ids) | genes["gene_symbol"].astype(str).str.strip().str.casefold().isin(rbp_symbols)
    genes["support_count"] = genes[["DEG_significant", "rMATS_significant", "DEXSeq_significant", "DTU_significant", "RBP_gene"]].sum(axis=1)

    event_support_columns = [
        "gene_key", "DEG_significant", "rMATS_significant", "DEXSeq_significant",
        "DTU_significant", "RBP_gene", "support_count",
    ]
    events = events.merge(genes[event_support_columns].drop_duplicates("gene_key"), on="gene_key", how="left")

    universe = set(genes["gene_key"])
    evidence_sets = {
        "DEG": set(genes.loc[genes["DEG_significant"], "gene_key"]),
        "rMATS": set(genes.loc[genes["rMATS_significant"], "gene_key"]),
        "DEXSeq": set(genes.loc[genes["DEXSeq_significant"], "gene_key"]),
        "DTU": set(genes.loc[genes["DTU_significant"], "gene_key"]),
        "RBP": set(genes.loc[genes["RBP_gene"], "gene_key"]),
    }
    summary_rows: list[dict[str, object]] = []
    for module, module_events in events.groupby("module", sort=True):
        module_genes = set(module_events["gene_key"])
        per_gene_events = (
            module_events.groupby("gene_key", observed=True)["event_key"]
            .nunique()
            .sort_values(ascending=False)
        )
        event_count = int(per_gene_events.sum())
        gene_count = int(len(per_gene_events))
        gene_event_shares = per_gene_events.to_numpy(float) / max(event_count, 1)
        positive_shares = gene_event_shares[gene_event_shares > 0]
        effective_gene_count = float(
            np.exp(-(positive_shares * np.log(positive_shares)).sum())
        )
        row: dict[str, object] = {
            "module": module,
            "n_events": event_count,
            "n_genes": gene_count,
            "events_per_unique_gene": event_count / max(gene_count, 1),
            "top_gene_event_fraction": (
                float(gene_event_shares[0]) if len(gene_event_shares) else np.nan
            ),
            "effective_gene_count": effective_gene_count,
            "effective_gene_fraction": effective_gene_count / max(gene_count, 1),
            "multi_event_gene_fraction": (
                float(per_gene_events.gt(1).mean()) if gene_count else np.nan
            ),
            "event_types": ";".join(sorted(set(module_events["event_type"].astype(str)))),
        }
        for label, selected in evidence_sets.items():
            overlap = module_genes & selected
            key = label.lower()
            row[f"{key}_gene_count"] = len(overlap)
            row[f"{key}_gene_fraction"] = len(overlap) / len(module_genes) if module_genes else 0.0
            row[f"{key}_enrichment_pvalue"] = _enrichment_pvalue(len(universe), len(selected), len(module_genes), len(overlap))
        summary_rows.append(row)
    module_summary = pd.DataFrame(summary_rows)
    if "module" not in module_summary:
        module_summary["module"] = pd.Series(dtype="string")
    for label in evidence_sets:
        column = f"{label.lower()}_enrichment_pvalue"
        module_summary[f"{label.lower()}_enrichment_FDR"] = bh_adjust(module_summary[column]) if column in module_summary else np.nan

    rbp_rows: list[dict[str, object]] = []
    for (gene_set, source), subset in rbp_sets.groupby(["gene_set", "source"], sort=True):
        set_ids = {normalize_gene_id(value) for value in subset["gene_id"] if normalize_gene_id(value)}
        set_symbols = {str(value).strip().casefold() for value in subset["gene_symbol"] if str(value).strip()}
        selected = set(
            genes.loc[
                genes["gene_id"].map(normalize_gene_id).isin(set_ids)
                | genes["gene_symbol"].astype(str).str.strip().str.casefold().isin(set_symbols),
                "gene_key",
            ]
        ) & universe
        for module, module_genes_frame in genes.groupby("module", sort=True):
            module_genes = set(module_genes_frame["gene_key"])
            overlap = module_genes & selected
            symbols = sorted(set(module_genes_frame.loc[module_genes_frame["gene_key"].isin(overlap), "gene_symbol"].astype(str)) - {""})
            rbp_rows.append(
                {
                    "module": module,
                    "gene_set": gene_set,
                    "source": source,
                    "network_gene_universe": len(universe),
                    "gene_set_size_in_universe": len(selected),
                    "module_gene_count": len(module_genes),
                    "overlap_count": len(overlap),
                    "overlap_genes": ";".join(symbols),
                    "pvalue": _enrichment_pvalue(len(universe), len(selected), len(module_genes), len(overlap)),
                }
            )
    rbp_enrichment = pd.DataFrame(rbp_rows)
    if not rbp_enrichment.empty:
        rbp_enrichment["FDR"] = bh_adjust(rbp_enrichment["pvalue"])
        rbp_enrichment = rbp_enrichment.sort_values(["FDR", "pvalue", "overlap_count"], ascending=[True, True, False], na_position="last")

    sources = pd.concat([support_sources, rbp_sources], ignore_index=True)
    return ModuleAnnotationResult(
        events=events.reset_index(drop=True),
        genes=genes.sort_values(["module", "support_count", "n_module_events"], ascending=[True, False, False]).reset_index(drop=True),
        module_summary=module_summary.sort_values("module").reset_index(drop=True),
        rbp_enrichment=rbp_enrichment.reset_index(drop=True),
        support_sources=sources,
    )


def render_module_annotation_plot(summary: pd.DataFrame, output_dir: Path) -> list[Path]:
    if summary.empty:
        return []
    fraction_columns = [name for name in summary.columns if name.endswith("_gene_fraction")]
    if not fraction_columns:
        return []
    matrix = summary.set_index("module")[fraction_columns].copy()
    matrix.columns = [name.removesuffix("_gene_fraction").upper() for name in matrix.columns]
    fig, ax = plt.subplots(figsize=(max(8.5, matrix.shape[1] * 1.3), max(5.5, matrix.shape[0] * 0.55 + 2.2)))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="crest", vmin=0, vmax=1, linewidths=0.4, ax=ax, cbar_kws={"label": "Fraction of unique module genes"})
    ax.set_title("CoSplice module evidence support", fontweight="bold", pad=16)
    ax.set_xlabel("Support layer (annotation only; does not define modules)")
    ax.set_ylabel("Co-splicing module")
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "module_evidence_support_heatmap.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return [path]
