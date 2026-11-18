"""Precision-weighted single-cell co-splicing network estimation and Signed TOM module discovery."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
import warnings
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import average, fcluster
from scipy.spatial.distance import squareform

from katana_sc_splicenet.input_contract import ValidatedSCUsageInputV2, validate_sc_usage_input_v2
from katana_sc_splicenet.qc_and_selection import audit_sc_events_v2, estimate_sc_observation_weights_v2


@dataclass(frozen=True)
class SCCoSpliceNetworkV2:
    """Fitted single-cell co-splicing network and module decomposition."""

    edges: pd.DataFrame
    modules: pd.DataFrame          # event_id, gene_id, module, kME, excess_bio_var, precision_mean
    event_audit: pd.DataFrame
    shrunk_corr: np.ndarray        # E x E symmetric cross-gene shrunk correlation (0 on diagonal & intra-gene)
    tom_matrix: np.ndarray         # E x E Signed Topological Overlap Matrix
    diagnostics: dict


def precision_pairwise_correlation_v2(
    values: np.ndarray,
    weights: np.ndarray | None = None,
    *,
    min_pairs: int = 5,
    min_eff_n: float = 4.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Shift-invariant, BLAS-vectorized weighted pairwise Pearson correlation with Kish effective N.

    Pair weight on unit c for pair (i, j) is p_{cij} = sqrt(w_{ci} * w_{cj}), allowing exact
    shared-subset weighted moments and Kish effective sample sizes via five dense GEMMs.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or min_pairs < 3 or np.isinf(x).any():
        raise ValueError("2D finite-or-NaN array and min_pairs >= 3 required")
    finite = np.isfinite(x)
    if weights is None:
        w = finite.astype(float)
    else:
        w = np.where(finite, np.maximum(np.asarray(weights, dtype=float), 0.0), 0.0)
    active = finite & (w > 0.0)
    m_bool = active.astype(float)
    n_obs = (m_bool.T @ m_bool).astype(int)

    # Pre-center each column by its own weighted mean to prevent catastrophic cancellation on near-constant columns
    w_col = w.sum(axis=0)
    x_safe = np.where(active, x, 0.0)
    col_mean = np.divide((w * x_safe).sum(axis=0), w_col, out=np.zeros(x.shape[1], dtype=float), where=w_col > 0.0)
    x0 = np.where(active, x - col_mean[None, :], 0.0)

    u = np.sqrt(w)
    z = u * x0
    zx = z * x0  # sqrt(w) * x0^2

    p_sum = u.T @ u
    p_sq_sum = w.T @ w
    n_eff = np.divide(p_sum * p_sum, p_sq_sum, out=np.zeros_like(p_sum), where=p_sq_sum > 0.0)

    denom = np.maximum(p_sum, 1e-15)
    mu_ij = (z.T @ u) / denom          # mu_{i | pair(i,j)}
    ex2_ij = (zx.T @ u) / denom        # E[x0_i^2 | pair(i,j)]
    var_ij = np.maximum(ex2_ij - mu_ij * mu_ij, 0.0)
    cov_ij = (z.T @ z) / denom - mu_ij * mu_ij.T
    scale = np.sqrt(var_ij * var_ij.T)

    valid_pair = (n_obs >= int(min_pairs)) & (n_eff >= float(min_eff_n)) & (scale > 1e-11) & (var_ij > 1e-12) & (var_ij.T > 1e-12)
    r = np.divide(cov_ij, scale, out=np.full_like(cov_ij, np.nan), where=valid_pair)
    r = np.clip(r, -1.0, 1.0)
    np.fill_diagonal(r, 1.0)
    return r, n_obs, n_eff


def precision_group_profiles_v2(
    values: np.ndarray,
    weights: np.ndarray,
    groups: np.ndarray,
    ordered_groups: list[str] | np.ndarray,
    *,
    min_cells: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate cell-level usage into precision-weighted group/donor/micro-state profiles."""
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    g = np.asarray(groups, dtype=str)
    finite = np.isfinite(v) & (w > 0.0)

    prof_rows: list[np.ndarray] = []
    wt_rows: list[np.ndarray] = []
    for grp in ordered_groups:
        sel = g == str(grp)
        if not np.any(sel):
            prof_rows.append(np.full(v.shape[1], np.nan, dtype=float))
            wt_rows.append(np.zeros(v.shape[1], dtype=float))
            continue
        f_sub = finite[sel]
        n_c = f_sub.sum(axis=0)
        w_sub = np.where(f_sub, w[sel], 0.0)
        v_sub = np.where(f_sub, v[sel], 0.0)
        w_tot = w_sub.sum(axis=0)
        ok = (n_c >= int(min_cells)) & (w_tot > 0.0)
        p_row = np.divide((w_sub * v_sub).sum(axis=0), w_tot, out=np.full(v.shape[1], np.nan), where=ok)
        prof_rows.append(p_row)
        wt_rows.append(np.where(ok, w_tot, 0.0))
    return np.stack(prof_rows, axis=0), np.stack(wt_rows, axis=0)


def stratified_split_halves_v2(strata: np.ndarray, seed: int) -> np.ndarray:
    """Split cells into balanced halves within each stratum."""
    s = np.asarray(strata, dtype=str)
    rng = np.random.default_rng(int(seed))
    half = np.zeros(len(s), dtype=bool)
    for label in np.unique(s):
        ix = np.flatnonzero(s == label)
        ix = rng.permutation(ix)
        half[ix[::2]] = True
    return half


def cross_half_correlation_v2(
    left_prof: np.ndarray,
    left_wt: np.ndarray,
    right_prof: np.ndarray,
    right_wt: np.ndarray,
    *,
    min_pairs: int = 4,
    min_eff_n: float = 3.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric precision-weighted cross-half correlation across aligned groups/donors."""
    n_events = left_prof.shape[1]
    joint_v = np.column_stack([left_prof, right_prof])
    joint_w = np.column_stack([left_wt, right_wt])
    r_all, _, neff_all = precision_pairwise_correlation_v2(
        joint_v, joint_w, min_pairs=min_pairs, min_eff_n=min_eff_n
    )
    c_ab = r_all[:n_events, n_events:]
    n_ab = neff_all[:n_events, n_events:]
    both = np.isfinite(c_ab) & np.isfinite(c_ab.T)
    c_sym = np.where(both, 0.5 * (c_ab + c_ab.T), np.nan)
    n_sym = np.where(both, np.minimum(n_ab, n_ab.T), np.nan)
    return c_sym, n_sym


def empirical_bayes_fisher_shrink_v2(
    r_vec: np.ndarray,
    n_eff_vec: np.ndarray,
    eligible_mask: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Shrink noisy pairwise correlations via Empirical-Bayes Fisher-z with Kish effective support."""
    r = np.asarray(r_vec, dtype=float)
    n = np.asarray(n_eff_vec, dtype=float)
    ok = np.asarray(eligible_mask, dtype=bool) & np.isfinite(r) & np.isfinite(n) & (n > 3.2)
    if not np.any(ok):
        return np.full_like(r, np.nan), {"prior_tau2": 1e-4, "eligible_pairs": 0}
    z = np.arctanh(np.clip(r, -0.999, 0.999))
    var_z = 1.0 / np.maximum(n - 3.0, 0.5)
    tau2 = max(float(np.mean((z[ok] * z[ok] - var_z[ok]))), 1e-4)
    shrink_factor = tau2 / (tau2 + var_z)
    z_post = np.where(ok, z * shrink_factor, np.nan)
    r_post = np.where(ok, np.tanh(z_post), np.nan)
    return r_post, {"prior_tau2": tau2, "eligible_pairs": int(ok.sum())}


def _build_training_microstates_v2(
    md: pd.DataFrame,
    support_totals_per_cell: np.ndarray,
    *,
    target_cells_per_bin: int = 6,
) -> np.ndarray:
    """Partition each (donor_id, condition, cell_type, subtype) stratum into deterministic depth-stratified micro-clusters.

    Used ONLY on training cells in `cell_state_distribution` mode to capture multi-state and
    fine-grained substate isoform shifts while averaging out single-cell binomial noise.
    """
    base_key = (
        md["donor_id"].astype(str)
        + "|" + md["condition"].astype(str)
        + "|" + md["cell_type"].astype(str)
        + "|" + md["subtype"].astype(str)
    ).to_numpy(dtype=str)
    micro = np.empty(len(md), dtype=object)
    for key in np.unique(base_key):
        ix = np.flatnonzero(base_key == key)
        n_bins = max(1, len(ix) // max(int(target_cells_per_bin), 3))
        if n_bins == 1:
            micro[ix] = f"{key}#0"
        else:
            # Order by cell total support so cells of comparable coverage group cleanly without touching usage values
            order = ix[np.argsort(support_totals_per_cell[ix], kind="mergesort")]
            bins = np.array_split(order, n_bins)
            for b_id, b_ix in enumerate(bins):
                micro[b_ix] = f"{key}#{b_id}"
    return micro.astype(str)


def compute_signed_tom_v2(
    shrunk_corr: np.ndarray,
    genes: np.ndarray,
    *,
    soft_power: int | None = None,
) -> tuple[np.ndarray, int]:
    """Compute cross-gene Signed Topological Overlap Matrix (TOM) zeroing intra-gene closure."""
    n = shrunk_corr.shape[0]
    same_gene = genes[:, None] == genes[None, :]
    r_clean = np.where(np.isfinite(shrunk_corr) & (~same_gene), shrunk_corr, 0.0)
    # Positive co-splicing adjacency (suppressing noisy negative compositional artifacts)
    pos_sim = np.clip(r_clean, 0.0, 1.0)
    np.fill_diagonal(pos_sim, 0.0)

    if soft_power is None:
        chosen_p = 4
        for p in (2, 3, 4, 5, 6):
            adj_try = np.power(pos_sim, p)
            k = adj_try.sum(axis=1)
            pos_k = k[k > 1e-6]
            if len(pos_k) >= 15:
                # Pick soft power balancing sparsity and connected hub structure
                mean_k = float(np.mean(pos_k))
                if mean_k <= max(12.0, 0.08 * n):
                    chosen_p = p
                    break
        soft_power = chosen_p

    adj = np.power(pos_sim, int(soft_power))
    np.fill_diagonal(adj, 0.0)
    shared = adj @ adj
    k_vec = adj.sum(axis=1)
    denom = np.minimum(k_vec[:, None], k_vec[None, :]) + 1.0 - adj
    tom = np.divide(shared + adj, denom, out=np.zeros_like(adj), where=denom > 1e-12)
    tom = np.clip(tom, 0.0, 1.0)
    tom[same_gene] = 0.0
    np.fill_diagonal(tom, 1.0)
    return tom, int(soft_power)


def _dynamic_cut_via_r_or_python(
    dissim_tom: np.ndarray,
    min_module_size: int = 10,
    deep_split: int = 4,
    use_r_dynamic_cut: bool = False,
    rscript_exe: str | None = None,
) -> np.ndarray:
    """Cluster 1 - TOM into co-splicing modules via pure-Python dynamicTreeCut (default, zero R needed!) or optional R."""
    import shutil
    n = dissim_tom.shape[0]
    if n < max(2 * min_module_size, 6):
        return np.zeros(n, dtype=int)

    d_clean = np.clip(0.5 * (dissim_tom + dissim_tom.T), 0.0, 1.0)
    np.fill_diagonal(d_clean, 0.0)

    r_bin = rscript_exe or (shutil.which("Rscript") if use_r_dynamic_cut else None)
    if use_r_dynamic_cut and r_bin and Path(r_bin).is_file():
        with tempfile.TemporaryDirectory(prefix="sc_cosplice_") as tmp_dir:
            tmp = Path(tmp_dir)
            bin_in, tsv_out, r_file = tmp / "dissim.bin", tmp / "labels.tsv", tmp / "cut.R"
            np.ascontiguousarray(d_clean, dtype="<f8").tofile(bin_in)
            r_code = (
                f'suppressPackageStartupMessages(library(dynamicTreeCut))\n'
                f'n <- {n}\n'
                f'raw <- readBin("{bin_in.as_posix()}", what="numeric", n=n*n, size=8, endian="little")\n'
                f'dmat <- matrix(raw, nrow=n, ncol=n)\n'
                f'hc <- hclust(as.dist(dmat), method="average")\n'
                f'mods <- cutreeDynamic(dendro=hc, distM=dmat, deepSplit={int(deep_split)}, pamRespectsDendro=FALSE, minClusterSize={int(min_module_size)})\n'
                f'write.table(data.frame(label=as.integer(mods)), file="{tsv_out.as_posix()}", sep="\t", row.names=FALSE, quote=FALSE)\n'
            )
            r_file.write_text(r_code, encoding="utf-8")
            proc = subprocess.run([r_bin, "--vanilla", str(r_file)], capture_output=True, text=True, timeout=120)
            if proc.returncode == 0 and tsv_out.is_file():
                lbl = pd.read_csv(tsv_out, sep="\t")["label"].to_numpy(dtype=int)
                if len(lbl) == n:
                    return lbl

    Z = average(squareform(d_clean, checks=False))
    try:
        from dynamicTreeCut import cutreeHybrid
        res = cutreeHybrid(
            Z, d_clean, minClusterSize=int(min_module_size), deepSplit=int(min(deep_split, 4)), pamRespectsDendro=False, verbose=0
        )
        lbl = np.asarray(res["labels"] if isinstance(res, dict) else res, dtype=int)
        if len(lbl) == n:
            return lbl
    except Exception:
        pass

    best_labels = np.zeros(n, dtype=int)
    best_score = -1.0
    for q in (0.70, 0.80, 0.85, 0.90, 0.94, 0.97):
        h = float(np.quantile(Z[:, 2], q))
        raw_cl = fcluster(Z, t=h, criterion="distance")
        vc = pd.Series(raw_cl).value_counts()
        valid_c = vc[(vc >= min_module_size) & (vc <= int(0.75 * n))].index.tolist()
        if not valid_c:
            continue
        mapped = np.zeros(n, dtype=int)
        for new_id, cid in enumerate(valid_c, start=1):
            mapped[raw_cl == cid] = new_id
        score = float(np.mean(mapped > 0)) * np.log1p(len(valid_c))
        if score > best_score:
            best_score = score
            best_labels = mapped
    return best_labels


def discover_sc_modules_v2(
    shrunk_corr: np.ndarray,
    usage_values: np.ndarray,
    obs_weights: np.ndarray,
    event_ids: np.ndarray,
    genes: np.ndarray,
    audit_df: pd.DataFrame,
    *,
    min_module_size: int = 10,
    soft_power: int | None = None,
    use_r_dynamic_cut: bool = False,
) -> tuple[pd.DataFrame, np.ndarray, int]:
    """Discover cross-gene co-splicing modules on Signed TOM and compute precision-weighted kME."""
    same_gene = genes[:, None] == genes[None, :]
    pos_vals = shrunk_corr[(~same_gene) & np.isfinite(shrunk_corr) & (shrunk_corr > 0.0)]
    if len(pos_vals) >= 10:
        q995 = max(float(np.quantile(pos_vals, 0.995)), 1e-6)
        scaled_corr = np.clip(np.where(np.isfinite(shrunk_corr), shrunk_corr / q995, 0.0), -1.0, 1.0)
    else:
        scaled_corr = shrunk_corr
    eff_power = 2 if soft_power is None else int(soft_power)
    tom, chosen_power = compute_signed_tom_v2(scaled_corr, genes, soft_power=eff_power)
    pos_scaled = np.where((~same_gene) & np.isfinite(scaled_corr), np.clip(scaled_corr, 0.0, 1.0), 0.0)
    hybrid_sim = np.maximum(tom, 0.90 * pos_scaled)
    np.fill_diagonal(hybrid_sim, 1.0)
    dissim = 1.0 - hybrid_sim
    np.fill_diagonal(dissim, 0.0)

    # Exclude completely disconnected nodes from clustering so noise singletons stay M0 (unassigned)
    max_tom = np.where(np.eye(len(event_ids), dtype=bool), 0.0, tom).max(axis=1)
    connected = max_tom > 1e-4
    raw_labels = np.zeros(len(event_ids), dtype=int)
    if int(connected.sum()) >= max(2 * min_module_size, 6):
        sub_idx = np.flatnonzero(connected)
        sub_lbl = _dynamic_cut_via_r_or_python(
            dissim[np.ix_(sub_idx, sub_idx)],
            min_module_size=min_module_size,
            use_r_dynamic_cut=use_r_dynamic_cut,
        )
        raw_labels[sub_idx] = sub_lbl

    # Enforce >=2 distinct genes per module and relabel M1, M2, ... by descending size
    vc = pd.Series(raw_labels[raw_labels > 0]).value_counts()
    clean_labels = np.zeros(len(event_ids), dtype=object)
    clean_labels[:] = "M0"
    mod_counter = 0
    for old_id in vc.index:
        members = raw_labels == int(old_id)
        if int(members.sum()) >= min_module_size and len(set(genes[members])) >= 2:
            mod_counter += 1
            clean_labels[members] = f"M{mod_counter}"

    # Compute leave-same-gene-out correlation to module eigengene/centroid (kME) on training cells
    kme = np.zeros(len(event_ids), dtype=float)
    same_gene = genes[:, None] == genes[None, :]
    r_clean = np.where(np.isfinite(shrunk_corr) & (~same_gene), shrunk_corr, 0.0)
    for m_idx in range(1, mod_counter + 1):
        m_name = f"M{m_idx}"
        in_mod = clean_labels == m_name
        for j in np.flatnonzero(in_mod):
            peers = in_mod & (~same_gene[j])
            if np.any(peers):
                kme[j] = float(np.mean(r_clean[j, peers]))

    audit_idx = audit_df.set_index("event_id").reindex(event_ids)
    mod_df = pd.DataFrame({
        "event_id": event_ids.astype(str),
        "gene_id": genes.astype(str),
        "module": clean_labels.astype(str),
        "kME": kme,
        "excess_bio_var": audit_idx["excess_bio_var"].to_numpy(dtype=float),
        "precision_mean": audit_idx["precision_mean"].to_numpy(dtype=float),
    })
    return mod_df, tom, chosen_power


def fit_sc_cosplice_v2(
    usage: pd.DataFrame,
    metadata: pd.DataFrame,
    event_genes: pd.Series | pd.DataFrame,
    *,
    support: pd.DataFrame | None = None,
    host_log_expr: pd.DataFrame | None = None,
    level: str = "cell_state_distribution",
    top_fraction: float = 0.05,
    edge_budget: int | None = None,
    eligible_pairs: np.ndarray | None = None,
    repeats: int = 8,
    random_seed: int = 20261117,
    min_module_size: int = 10,
    compute_modules: bool = True,
    use_r_dynamic_cut: bool = False,
    require_qc_whitelist: bool = True,
) -> SCCoSpliceNetworkV2:
    """Fit Standalone SC CoSplice v2 precision-weighted co-splicing network and modules.

    Supported `level` estimands:
      - `"between_donor"`: Multi-donor (>=5 training donors within a single cell_type/condition)
        precision-weighted donor profile network with split-half cross-correlation q25 and
        Empirical-Bayes Kish-support shrinkage.
      - `"cell_state_distribution"`: Multi-state / Smart-seq3 precision-weighted micro-state +
        cross-cell shrinkage usage network (works with >=1 training donor when >=2 states/subtypes
        or >=30 cells exist).
      - `"within_state_cell"`: Precision-weighted state-residual cell-level co-splicing network.
    """
    if level not in {"between_donor", "cell_state_distribution", "within_state_cell"}:
        raise ValueError(f"Unsupported network level: {level!r}")
    if repeats < 4 or not (0.0 < top_fraction < 1.0):
        raise ValueError("repeats >= 4 and 0 < top_fraction < 1 required")

    bundle: ValidatedSCUsageInputV2 = validate_sc_usage_input_v2(
        usage,
        metadata,
        event_genes,
        support=support,
        host_log_expr=host_log_expr,
        require_qc_whitelist=require_qc_whitelist,
    )
    md = bundle.metadata
    u_df = bundle.usage
    s_df = bundle.support
    genes = bundle.event_genes.to_numpy(dtype=str)
    event_ids = u_df.columns.to_numpy(dtype=str)
    n_events = len(event_ids)

    audit_df = audit_sc_events_v2(u_df, bundle.event_genes, support=s_df, min_cells=3, min_total_support=3.0)
    u_vals = u_df.to_numpy(dtype=float)
    s_vals = s_df.to_numpy(dtype=float) if s_df is not None else None
    obs_w = estimate_sc_observation_weights_v2(u_vals, s_vals, audit_df["excess_bio_var"].to_numpy(dtype=float))

    ia, ib = np.triu_indices(n_events, k=1)
    cross_gene_pair = genes[ia] != genes[ib]
    if eligible_pairs is not None:
        ep = np.asarray(eligible_pairs, dtype=bool)
        if ep.shape != (n_events, n_events) or not np.array_equal(ep, ep.T):
            raise ValueError("eligible_pairs must be a symmetric E x E boolean matrix")
        cross_gene_pair = cross_gene_pair & ep[ia, ib]

    donors = sorted(md["donor_id"].unique())
    split_strata = (
        md["donor_id"] + "|" + md["condition"] + "|" + md["cell_type"] + "|" + md["subtype"]
    ).to_numpy(dtype=str)

    if level == "between_donor":
        if md["cell_type"].nunique() != 1 or md["condition"].nunique() != 1:
            raise ValueError("between_donor level requires single cell_type and condition; stratify before fitting")
        if len(donors) < 5:
            raise ValueError("between_donor level requires >= 5 distinct training donors")

        full_prof, full_wt = precision_group_profiles_v2(u_vals, obs_w, md["donor_id"].to_numpy(), donors, min_cells=3)
        full_r, full_n, full_neff = precision_pairwise_correlation_v2(full_prof, full_wt, min_pairs=5, min_eff_n=3.5)

        cross_draws: list[np.ndarray] = []
        neff_draws: list[np.ndarray] = []
        left_draws: list[np.ndarray] = []
        right_draws: list[np.ndarray] = []
        for rep in range(repeats):
            h = stratified_split_halves_v2(split_strata, random_seed + rep)
            lp, lw = precision_group_profiles_v2(u_vals[h], obs_w[h], md["donor_id"].to_numpy()[h], donors, min_cells=2)
            rp, rw = precision_group_profiles_v2(u_vals[~h], obs_w[~h], md["donor_id"].to_numpy()[~h], donors, min_cells=2)
            lr, _, _ = precision_pairwise_correlation_v2(lp, lw, min_pairs=4, min_eff_n=3.0)
            rr, _, _ = precision_pairwise_correlation_v2(rp, rw, min_pairs=4, min_eff_n=3.0)
            c_sym, n_sym = cross_half_correlation_v2(lp, lw, rp, rw, min_pairs=4, min_eff_n=3.2)
            left_draws.append(lr[ia, ib])
            right_draws.append(rr[ia, ib])
            cross_draws.append(c_sym[ia, ib])
            neff_draws.append(n_sym[ia, ib])

        valid_splits = np.isfinite(cross_draws).sum(axis=0)
        min_valid = int(np.ceil(0.75 * repeats))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            q25 = np.nanquantile(cross_draws, 0.25, axis=0)
            eff_sup = np.nanmedian(neff_draws, axis=0)

        elig = cross_gene_pair & np.isfinite(full_r[ia, ib]) & (valid_splits >= min_valid) & np.isfinite(q25) & (eff_sup > 3.2)
        shrunk_vec, shrink_info = empirical_bayes_fisher_shrink_v2(q25, eff_sup, elig)
        raw_vec = full_r[ia, ib]
        sup_vec = full_n[ia, ib].astype(float)

        # Split-half stability gate for default budget cap when edge_budget is None
        l_arr = np.asarray(left_draws, dtype=float)
        r_arr = np.asarray(right_draws, dtype=float)
        v_lr = np.isfinite(l_arr) & np.isfinite(r_arr)
        stab_ok = v_lr & (l_arr > 0.0) & (r_arr > 0.0) & (l_arr >= 0.5 * raw_vec[None, :]) & (r_arr >= 0.5 * raw_vec[None, :])
        n_lr = v_lr.sum(axis=0)
        stability = np.divide(stab_ok.sum(axis=0), n_lr, out=np.zeros_like(raw_vec), where=n_lr > 0)
        gate_pass = elig & (raw_vec > 0.0) & (q25 > 0.0) & (stability >= 0.75)

    elif level == "cell_state_distribution":
        cell_tot_sup = obs_w.sum(axis=1)
        micro_groups = _build_training_microstates_v2(md, cell_tot_sup, target_cells_per_bin=6)
        ordered_micro = sorted(np.unique(micro_groups))

        full_cell_r, full_cell_n, full_cell_neff = precision_pairwise_correlation_v2(
            u_vals, obs_w, min_pairs=8, min_eff_n=5.0
        )
        if len(ordered_micro) >= 6:
            m_prof, m_wt = precision_group_profiles_v2(u_vals, obs_w, micro_groups, ordered_micro, min_cells=2)
            full_micro_r, _, full_micro_neff = precision_pairwise_correlation_v2(
                m_prof, m_wt, min_pairs=5, min_eff_n=4.0
            )
        else:
            full_micro_r = full_cell_r
            full_micro_neff = full_cell_neff

        # Split-half stability across stratified cell halves
        split_r_draws: list[np.ndarray] = []
        for rep in range(repeats):
            h = stratified_split_halves_v2(split_strata, random_seed + rep)
            r_h1, _, _ = precision_pairwise_correlation_v2(u_vals[h], obs_w[h], min_pairs=5, min_eff_n=3.5)
            r_h2, _, _ = precision_pairwise_correlation_v2(u_vals[~h], obs_w[~h], min_pairs=5, min_eff_n=3.5)
            if len(ordered_micro) >= 6:
                mp1, mw1 = precision_group_profiles_v2(u_vals[h], obs_w[h], micro_groups[h], ordered_micro, min_cells=1)
                mp2, mw2 = precision_group_profiles_v2(u_vals[~h], obs_w[~h], micro_groups[~h], ordered_micro, min_cells=1)
                c_m, _ = cross_half_correlation_v2(mp1, mw1, mp2, mw2, min_pairs=4, min_eff_n=3.2)
                comb_rep = np.where(
                    np.isfinite(c_m) & np.isfinite(r_h1) & np.isfinite(r_h2),
                    0.65 * c_m + 0.35 * 0.5 * (r_h1 + r_h2),
                    np.where(np.isfinite(c_m), c_m, 0.5 * (r_h1 + r_h2)),
                )
            else:
                comb_rep = 0.5 * (r_h1 + r_h2)
            split_r_draws.append(comb_rep[ia, ib])

        valid_splits = np.isfinite(split_r_draws).sum(axis=0)
        min_valid = int(np.ceil(0.75 * repeats))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            q25 = np.nanquantile(split_r_draws, 0.25, axis=0)

        raw_blend = np.where(
            np.isfinite(full_micro_r) & np.isfinite(full_cell_r),
            0.65 * full_micro_r + 0.35 * full_cell_r,
            np.where(np.isfinite(full_micro_r), full_micro_r, full_cell_r),
        )
        raw_vec = raw_blend[ia, ib]
        # Effective support bounded by micro-state count when state correlation dominates, preventing cell pseudoreplication
        eff_sup = np.where(
            np.isfinite(full_micro_neff[ia, ib]),
            np.minimum(full_cell_neff[ia, ib], np.maximum(full_micro_neff[ia, ib] * 1.5, 5.0)),
            full_cell_neff[ia, ib],
        )
        elig = cross_gene_pair & np.isfinite(raw_vec) & (valid_splits >= min_valid) & np.isfinite(q25) & (eff_sup > 3.5)
        shrunk_vec, shrink_info = empirical_bayes_fisher_shrink_v2(q25, eff_sup, elig)
        sup_vec = full_cell_n[ia, ib].astype(float)
        s_arr = np.asarray(split_r_draws, dtype=float)
        stability = np.divide((np.isfinite(s_arr) & (s_arr > 0.0)).sum(axis=0), np.maximum(valid_splits, 1))
        gate_pass = elig & (raw_vec > 0.0) & (q25 > 0.0) & (stability >= 0.75)

    else:  # within_state_cell
        ordered_states = sorted(np.unique(split_strata))
        st_prof, _ = precision_group_profiles_v2(u_vals, obs_w, split_strata, ordered_states, min_cells=3)
        st_map = {s: idx for idx, s in enumerate(ordered_states)}
        row_idx = np.array([st_map[s] for s in split_strata], dtype=int)
        u_res = u_vals - st_prof[row_idx]

        full_r, full_n, full_neff = precision_pairwise_correlation_v2(u_res, obs_w, min_pairs=8, min_eff_n=5.0)
        split_draws: list[np.ndarray] = []
        for rep in range(repeats):
            h = stratified_split_halves_v2(split_strata, random_seed + rep)
            r1, _, _ = precision_pairwise_correlation_v2(u_res[h], obs_w[h], min_pairs=5, min_eff_n=3.5)
            r2, _, _ = precision_pairwise_correlation_v2(u_res[~h], obs_w[~h], min_pairs=5, min_eff_n=3.5)
            split_draws.append(0.5 * (r1[ia, ib] + r2[ia, ib]))

        valid_splits = np.isfinite(split_draws).sum(axis=0)
        min_valid = int(np.ceil(0.75 * repeats))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            q25 = np.nanquantile(split_draws, 0.25, axis=0)
        raw_vec = full_r[ia, ib]
        eff_sup = full_neff[ia, ib]
        elig = cross_gene_pair & np.isfinite(raw_vec) & (valid_splits >= min_valid) & np.isfinite(q25) & (eff_sup > 4.0)
        shrunk_vec, shrink_info = empirical_bayes_fisher_shrink_v2(q25, eff_sup, elig)
        sup_vec = full_n[ia, ib].astype(float)
        s_arr = np.asarray(split_draws, dtype=float)
        stability = np.divide((np.isfinite(s_arr) & (s_arr > 0.0)).sum(axis=0), np.maximum(valid_splits, 1))
        gate_pass = elig & (raw_vec > 0.0) & (q25 > 0.0) & (stability >= 0.75)

    # Assemble full symmetric cross-gene shrunk correlation matrix
    shrunk_mat = np.zeros((n_events, n_events), dtype=float)
    fin_shrunk = np.isfinite(shrunk_vec) & cross_gene_pair
    shrunk_mat[ia[fin_shrunk], ib[fin_shrunk]] = shrunk_vec[fin_shrunk]
    shrunk_mat[ib[fin_shrunk], ia[fin_shrunk]] = shrunk_vec[fin_shrunk]

    candidate_idx = np.flatnonzero(elig & np.isfinite(shrunk_vec) & (shrunk_vec > 0.0))
    if edge_budget is None:
        max_frac_budget = int(np.ceil(int(cross_gene_pair.sum()) * top_fraction))
        budget = min(max_frac_budget, int(gate_pass.sum()))
    else:
        budget = int(edge_budget)
    if budget > 0 and len(candidate_idx) > 0:
        take = min(budget, len(candidate_idx))
        chosen = candidate_idx[np.argsort(shrunk_vec[candidate_idx], kind="mergesort")[-take:]][::-1]
    else:
        chosen = np.array([], dtype=int)

    edges_df = pd.DataFrame({
        "network_level": level,
        "edge_index": chosen,
        "source_event": event_ids[ia[chosen]],
        "target_event": event_ids[ib[chosen]],
        "source_gene": genes[ia[chosen]],
        "target_gene": genes[ib[chosen]],
        "raw_correlation": raw_vec[chosen],
        "split_q25": q25[chosen],
        "shrunk_correlation": shrunk_vec[chosen],
        "effective_support": eff_sup[chosen],
        "observed_pairs": sup_vec[chosen],
        "split_stability": stability[chosen],
    })

    if compute_modules:
        mod_df, tom_mat, chosen_power = discover_sc_modules_v2(
            shrunk_mat,
            u_vals,
            obs_w,
            event_ids,
            genes,
            audit_df,
            min_module_size=min_module_size,
            use_r_dynamic_cut=use_r_dynamic_cut,
        )
    else:
        mod_df = pd.DataFrame({
            "event_id": event_ids,
            "gene_id": genes,
            "module": "M0",
            "kME": 0.0,
            "excess_bio_var": audit_df["excess_bio_var"].to_numpy(dtype=float),
            "precision_mean": audit_df["precision_mean"].to_numpy(dtype=float),
        })
        tom_mat = np.eye(n_events, dtype=float)
        chosen_power = 0

    covered_nodes = set(edges_df["source_event"]).union(set(edges_df["target_event"]))
    n_mods = int(mod_df.loc[mod_df["module"].ne("M0"), "module"].nunique())
    diagnostics = {
        "status": "complete" if (len(edges_df) > 0 or n_mods > 0) else "no_reliable_edges",
        "method": "sc_cosplice_v2_precision_tom",
        "network_level": level,
        "training_cells": int(len(md)),
        "training_donors": donors,
        "events": int(n_events),
        "genes": int(len(set(genes))),
        "candidate_cross_gene_pairs": int(cross_gene_pair.sum()),
        "selected_edges": int(len(edges_df)),
        "edge_node_coverage": float(len(covered_nodes) / n_events) if n_events else 0.0,
        "modules_computed": bool(compute_modules),
        "discovered_modules": n_mods,
        "module_assigned_events": int(mod_df["module"].ne("M0").sum()),
        "tom_soft_power": int(chosen_power),
        "split_repeats": int(repeats),
        **shrink_info,
    }
    return SCCoSpliceNetworkV2(
        edges=edges_df,
        modules=mod_df,
        event_audit=audit_df,
        shrunk_corr=shrunk_mat,
        tom_matrix=tom_mat,
        diagnostics=diagnostics,
    )
