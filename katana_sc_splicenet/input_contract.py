"""Smart-seq3 / full-length single-cell splicing input contract and QC whitelist validation."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ValidatedSCUsageInputV2:
    """Immutable validated single-cell splicing usage bundle."""

    usage: pd.DataFrame          # cell_id x event_id, float64 in [0, 1] or NaN
    support: pd.DataFrame | None # cell_id x event_id, float64 >= 0 (0 -> usage is NaN)
    metadata: pd.DataFrame       # QC-passed cells aligned to usage.index
    event_genes: pd.Series       # event_id -> gene_id
    host_log_expr: pd.DataFrame | None = None  # optional cell_id x event_id host-gene log1p expression


def detect_jeffreys_or_pseudocount_artifact(
    usage_values: np.ndarray,
    support_values: np.ndarray | None = None,
    *,
    tol: float = 1e-6,
) -> bool:
    """Detect forbidden (a + 0.5)/(n + 1) Jeffreys shrinkage passed as observed usage."""
    u = np.asarray(usage_values, dtype=float)
    finite = np.isfinite(u)
    if not finite.any():
        return False
    if support_values is not None:
        s = np.asarray(support_values, dtype=float)
        zero_sup_finite = (s <= 0) & finite
        if np.any(zero_sup_finite):
            return True
        s1 = (np.abs(s - 1.0) <= tol) & finite
        if np.sum(s1) >= 20:
            u1 = u[s1]
            near_quarter = np.mean((np.abs(u1 - 0.25) <= 1e-3) | (np.abs(u1 - 0.75) <= 1e-3))
            near_binary = np.mean((np.abs(u1 - 0.0) <= 1e-3) | (np.abs(u1 - 1.0) <= 1e-3))
            if near_quarter > 0.80 and near_binary < 0.05:
                return True
    return False


def validate_sc_usage_input_v2(
    usage: pd.DataFrame,
    metadata: pd.DataFrame,
    event_genes: pd.Series | pd.DataFrame,
    *,
    support: pd.DataFrame | None = None,
    host_log_expr: pd.DataFrame | None = None,
    require_qc_whitelist: bool = True,
) -> ValidatedSCUsageInputV2:
    """Validate single-cell usage matrix, explicit QC whitelist, support denominators, and gene mapping.

    Extra rows in `usage`, `support`, or `host_log_expr` outside the approved `metadata`
    cell whitelist are never accessed downstream.
    """
    if not isinstance(usage, pd.DataFrame) or usage.empty:
        raise ValueError("Non-empty cell x event usage DataFrame required")
    if usage.index.has_duplicates or usage.columns.has_duplicates:
        raise ValueError("usage must have unique cell_id index and unique event_id columns")

    if isinstance(event_genes, pd.DataFrame):
        if not {"event_id", "gene_id"} <= set(event_genes.columns):
            raise ValueError("event_genes DataFrame requires event_id and gene_id columns")
        if event_genes["event_id"].duplicated().any():
            raise ValueError("Duplicate event_id in event_genes")
        eg = event_genes.set_index("event_id")["gene_id"]
    else:
        eg = event_genes.copy()
        if eg.index.has_duplicates:
            raise ValueError("Duplicate event_id in event_genes Series")

    eg = eg.reindex(usage.columns)
    if eg.isna().any() or eg.astype(str).str.strip().eq("").any():
        raise ValueError("Every event/isoform must map to a valid non-empty gene_id")
    eg = eg.astype(str)

    if not isinstance(metadata, pd.DataFrame) or metadata.empty:
        raise ValueError("Non-empty cell metadata DataFrame required")
    md = metadata.copy()
    if "cell_id" not in md.columns:
        if md.index.name == "cell_id":
            md = md.reset_index()
        else:
            raise ValueError("metadata must include explicit cell_id column")
    md["cell_id"] = md["cell_id"].astype(str)
    if md["cell_id"].duplicated().any():
        raise ValueError("Duplicate cell_id entries in metadata")

    required_cols = {"cell_id", "donor_id", "cell_type", "condition"}
    for alias_src, alias_dst in (("donor", "donor_id"), ("animal_ID", "donor_id"), ("treatment", "condition")):
        if alias_dst not in md.columns and alias_src in md.columns:
            md[alias_dst] = md[alias_src]
    if not required_cols <= set(md.columns):
        missing = sorted(required_cols - set(md.columns))
        raise ValueError(f"metadata missing required columns: {missing}")
    if md[list(required_cols)].isna().any().any():
        raise ValueError("metadata contains missing values in required cell_id/donor_id/cell_type/condition columns")

    if require_qc_whitelist:
        if "qc_pass" not in md.columns:
            raise ValueError("Explicit qc_pass column required when require_qc_whitelist=True")
        if md["qc_pass"].isna().any():
            raise ValueError("Missing qc_pass decisions cannot be coerced to True")
        md = md.loc[md["qc_pass"].eq(True)].copy()
        if md.empty:
            raise ValueError("Zero cells remain after applying qc_pass==True whitelist")

    if "subtype" not in md.columns:
        md["subtype"] = md["cell_type"].astype(str)
    for col in ("donor_id", "cell_type", "condition", "subtype"):
        md[col] = md[col].astype(str)

    usage_idx_str = usage.index.astype(str)
    if len(set(usage_idx_str)) != len(usage_idx_str):
        raise ValueError("String-coerced usage cell_id index has duplicates")
    u_aligned = usage.copy()
    u_aligned.index = usage_idx_str
    if not md["cell_id"].isin(u_aligned.index).all():
        missing_cells = sorted(set(md["cell_id"]) - set(u_aligned.index))[:5]
        raise ValueError(f"Whitelisted cells absent from usage matrix, e.g. {missing_cells}")
    u_aligned = u_aligned.loc[md["cell_id"]].astype(float)

    u_vals = u_aligned.to_numpy(dtype=float)
    finite = np.isfinite(u_vals)
    if np.isinf(u_vals).any() or np.any((u_vals[finite] < -1e-9) | (u_vals[finite] > 1.0 + 1e-9)):
        raise ValueError("Observed usage values must lie in [0, 1] or be NaN")
    if finite.any():
        u_vals = np.where(finite, np.clip(u_vals, 0.0, 1.0), np.nan)
        u_aligned = pd.DataFrame(u_vals, index=u_aligned.index, columns=u_aligned.columns)

    s_aligned: pd.DataFrame | None = None
    if support is not None:
        s_tmp = support.copy()
        s_tmp.index = s_tmp.index.astype(str)
        if not md["cell_id"].isin(s_tmp.index).all() or not set(u_aligned.columns) <= set(s_tmp.columns):
            raise ValueError("support matrix must cover all whitelisted cells and usage events")
        s_aligned = s_tmp.loc[md["cell_id"], u_aligned.columns].astype(float)
        s_vals = s_aligned.to_numpy(dtype=float)
        if np.isnan(s_vals).any() or np.isinf(s_vals).any() or np.any(s_vals < 0.0):
            raise ValueError("support denominators must be finite non-negative numbers (0 for unmeasured)")
        if detect_jeffreys_or_pseudocount_artifact(u_vals, s_vals):
            raise ValueError(
                "Forbidden pseudocount/Jeffreys usage detected (e.g., finite usage where support==0 "
                "or (a+0.5)/(n+1) values at support==1); pass unshrunk observed usage a/t with NaN at support==0"
            )
    else:
        if detect_jeffreys_or_pseudocount_artifact(u_vals, None):
            raise ValueError("Forbidden pseudocount usage artifact detected")

    h_aligned: pd.DataFrame | None = None
    if host_log_expr is not None:
        h_tmp = host_log_expr.copy()
        h_tmp.index = h_tmp.index.astype(str)
        if not md["cell_id"].isin(h_tmp.index).all() or not set(u_aligned.columns) <= set(h_tmp.columns):
            raise ValueError("host_log_expr must cover all whitelisted cells and usage events")
        h_aligned = h_tmp.loc[md["cell_id"], u_aligned.columns].astype(float)

    md = md.reset_index(drop=True)
    return ValidatedSCUsageInputV2(
        usage=u_aligned,
        support=s_aligned,
        metadata=md,
        event_genes=eg,
        host_log_expr=h_aligned,
    )
