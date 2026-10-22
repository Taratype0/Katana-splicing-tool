"""Matched-gene, whole-pair swaps with maximum cluster-mass correction.

This is a new conditional sequence association analysis, not the saved native
rMAPS3 rank-sum test. All models, methods, regions, widths and AS types in one
comparison share the maximum-cluster null. Windows have a separate BH family.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from katana_next.science.motiftrace_events import score_positions, validate_models
from katana_next.science.motiftrace_custom import fixed_pattern
from katana_next.science.statistics import bh_adjust

CLUSTER_Z = float(norm.isf(.01))
POLICY = {
    "schema": "motiftrace.event_consequences.v1",
    "foreground": "exact selected read-supported DAS; no refit",
    "background": "same AS type; measured non-DAS FDR>0.1, abs(dPSI)<0.05; exclude every significant/selected gene",
    "matching": "unique foreground/control genes, affected length +/-20%, GC +/-0.05, denominator PSI +/-0.15, minimum group mean reads within two-fold",
    "RNA": "reference local event forms, not phased full transcripts; RNA-oriented half-open genomic blocks",
    "short_sequences": "complete-sequence scalar at every eligible width; unavailable windows remain NA",
    "methods": "finite nonoverlapping counts and overlapping union coverage; PWM log2(1+mean odds), no hit cutoff",
    "PWM_candidates": "relative likelihood >=0.8 for position display ONLY; no site P/affinity claim; not used for enrichment",
    "PWM_background": "fixed A/C/G/T=0.25; same in foreground/control; pseudocount (20*PWM+0.1*background)/20.1",
    "statistic": "sum signed ranks / sqrt(sum squared ranks), zero differences excluded; paired whole-event vectors swapped",
    "exchangeability": "conditional on motif-blind matching; all windows/models of a gene pair move together; observational, not randomized biology",
    "cluster_forming": "one-sided positive Z > norm.isf(0.01), adjacent starts at step=1 within the SAME region/width/model/method/type/direction",
    "cluster_mass": "sum(Z - cluster_forming_Z); repeat cluster search on each null swap",
    "cluster_family": "maximum over all supplied models/methods/types/directions/regions/window widths in ONE comparison; corrected P is FWER, NOT FDR or individual-base significance",
    "window_family": "separate BH over all estimable new paired-swap scalar/window P values; historical P/FDR unchanged",
    "P": "Monte Carlo (1 + exceedances)/(B+1); never zero",
    "claim_boundary": "sequence opportunity/association, not measured occupancy, affinity, causality or superiority",
    "quantification_rerun": False,
}


def cluster_runs(z, threshold=CLUSTER_Z):
    eligible = np.flatnonzero(np.isfinite(z) & (z > threshold))
    if not len(eligible):
        return []
    parts = np.split(eligible, np.flatnonzero(np.diff(eligible) != 1)+1)
    return [(int(part[0]), int(part[-1]+1), float(np.sum(z[part]-threshold))) for part in parts]


def _signed_ranks(differences, minimum):
    valid = np.isfinite(differences)
    n = valid.sum(axis=0)
    ranks = np.zeros(differences.shape, dtype=float)
    for col in range(differences.shape[1]):
        chosen = valid[:, col] & (differences[:, col] != 0)
        ranks[chosen, col] = rankdata(np.abs(differences[chosen, col]), method="average") * np.sign(differences[chosen, col])
    denominator = np.sqrt((ranks*ranks).sum(axis=0))
    denominator[denominator == 0] = 1
    z = ranks.sum(axis=0)/denominator
    z[n < minimum] = np.nan
    return ranks/denominator, z, n


def paired_cluster_family(datasets, pair_ids, *, permutations=999, seed=20261015,
                           minimum_pairs=20, check_cancelled=lambda: None):
    """Input datasets carry complete pair x ordered-position matrices and metadata."""
    if not 99 <= permutations <= 9999 or minimum_pairs < 2:
        raise ValueError("Require 99-9999 fixed permutations and at least two independent pairs")
    if len(pair_ids) != len(set(pair_ids)):
        raise ValueError("Gene-pair identities must be unique")
    if permutations*len(pair_ids) > 8_000_000:
        raise ValueError("Explicit whole-pair permutation memory budget exceeded")
    if any(len(d["differences"]) != len(d["pair_ids"]) for d in datasets):
        raise ValueError("Pair IDs/matrices disagree")
    if any(len(d["pair_ids"]) != len(set(d["pair_ids"])) for d in datasets):
        raise ValueError("Dependent duplicate pairs within a dataset")
    if not set(p for d in datasets for p in d["pair_ids"]) <= set(pair_ids):
        raise ValueError("Unknown gene pair in cluster family")
    rng = np.random.default_rng(seed)
    signs = rng.choice(np.array([-1., 1.]), size=(permutations, len(pair_ids)))
    pair_lookup = {key: i for i, key in enumerate(pair_ids)}
    null_max = np.zeros(permutations)
    windows, clusters = [], []
    for data in datasets:
        check_cancelled()
        difference = np.asarray(data["differences"], float)
        if difference.ndim != 2 or len(data["positions"]) != difference.shape[1]:
            raise ValueError("Expected pair x position matrix")
        if data["positions"] != list(range(len(data["positions"]))):
            raise ValueError("Cluster search requires the complete step-one position grid including NA")
        ranks, observed, n = _signed_ranks(difference, minimum_pairs)
        indices = [pair_lookup[p] for p in data["pair_ids"]]
        exceed = np.zeros(len(observed), dtype=int)
        is_scalar = data["meta"]["window_nt"] == 0
        for start in range(0, permutations, 32):
            check_cancelled()
            # Do not invoke a BLAS thread pool or change global thread settings
            # in a server that may be running the user's other analyses.
            values = np.einsum("bi,ij->bj", signs[start:start+32, indices], ranks, optimize=False)
            values[:, n < minimum_pairs] = np.nan
            exceed += np.sum(values >= observed-1e-12, axis=0)
            if not is_scalar:
                for j, value in enumerate(values):
                    runs = cluster_runs(value)
                    if runs:
                        null_max[start+j] = max(null_max[start+j], max(run[2] for run in runs))
        for position in data["positions"]:
            finite = np.isfinite(difference[:, position])
            row = {**data["meta"], "window_start0": position if not is_scalar else -1,
                   "window_end0": position+data["meta"]["window_nt"] if not is_scalar else -1,
                   "n_matched_gene_pairs": int(n[position]), "paired_rank_Z": observed[position],
                   "mean_paired_score_difference": float(np.mean(difference[finite, position])) if finite.any() else np.nan,
                   "raw_P_new_paired_swaps": (1+exceed[position])/(permutations+1) if np.isfinite(observed[position]) else np.nan,
                   "status": "estimable" if np.isfinite(observed[position]) else "insufficient_independent_pairs"}
            windows.append(row)
        if not is_scalar:
            for a, b, mass in cluster_runs(observed):
                contributors = set()
                for i, name in enumerate(data["pair_ids"]):
                    if np.any(np.isfinite(difference[i, a:b]) & (difference[i, a:b] > 0)):
                        contributors.add(name)
                clusters.append({**data["meta"], "first_window_start0": a, "last_window_start0": b-1,
                    "covered_region_start0": a, "covered_region_end0": b-1+data["meta"]["window_nt"],
                    "number_correlated_windows": b-a, "cluster_mass": mass,
                    "minimum_matched_pairs_in_cluster": int(n[a:b].min()),
                    "supporting_independent_gene_pairs": len(contributors), "supporting_pair_ids_json": json.dumps(sorted(contributors))})
    window_table = pd.DataFrame(windows)
    if len(window_table):
        p = window_table.raw_P_new_paired_swaps.to_numpy()
        q = np.full(len(p), np.nan)
        valid = np.isfinite(p)
        if valid.any():
            q[valid] = bh_adjust(p[valid])
        window_table["q_new_window_family"] = q
    cluster_table = pd.DataFrame(clusters)
    if len(cluster_table):
        cluster_table["cluster_P_FWER_new_family"] = [
            (1+np.sum(null_max >= mass-1e-12))/(permutations+1) for mass in cluster_table.cluster_mass]
    return window_table, cluster_table, null_max


def _profile(sequence, model, method, width):
    position, _, motif_width = score_positions(sequence, model)
    span = max(0, len(sequence)-width+1) if width else 1
    output = np.full(span, np.nan)
    pattern = fixed_pattern(model["pattern"])[0] if model["kind"] == "consensus" else ""
    for start in range(span):
        end = start+width if width else len(sequence)
        word = sequence[start:end]
        if len(word) < motif_width:
            continue
        valid = np.array([set(word[j:j+motif_width]) <= set("ACGT") for j in range(len(word)-motif_width+1)])
        if not valid.any():
            continue
        if method == "PWM_mean_odds":
            output[start] = np.log2(1+np.nanmean(position[start:end-motif_width+1]))
        elif method == "finite_nonoverlap_counts":
            output[start] = len(re.findall(pattern, word))
        else:
            covered, available = np.zeros(len(word), bool), np.zeros(len(word), bool)
            for j in np.flatnonzero(valid):
                available[j:j+motif_width] = True
            for j in np.flatnonzero(position[start:end-motif_width+1] == 1):
                covered[j:j+motif_width] = True
            output[start] = covered.sum()/available.sum()
    return output


def build_paired_datasets(rna, matching, models, *, widths=(10, 20, 50),
                          max_tests=100_000, max_matrix_cells=8_000_000, check_cancelled=lambda: None):
    validate_models(models)
    if not widths or any(not isinstance(w, int) or w < 1 for w in widths):
        raise ValueError("Explicit positive integer window widths required")
    pairs = matching[matching.status.eq("matched")].copy()
    if pairs.foreground_gene.duplicated().any() or pairs.background_gene.duplicated().any() or set(pairs.foreground_gene) & set(pairs.background_gene):
        raise ValueError("Dependent/reused foreground/background genes are not exchangeable independent pairs")
    pairs["pair_id"] = pairs.foreground_event + "::" + pairs.background_event
    # Junction consequences are traced per actual instance, not pooled with
    # body/regulatory profiles whose RNA origins differ for short flanks.
    regions = rna[rna.layer.isin(["affected_sequence", "regulatory_pre_mRNA"])]
    lookup = {(r["event_id"], r["region"]): r["sequence"] for r in regions.to_dict("records")}
    datasets, tests, cells = [], 0, 0
    for (kind, direction), group in pairs.groupby(["event_type", "direction"], sort=True):
        for region in sorted(set(regions.loc[regions.event_type.eq(kind), "region"])):
            seq_pairs = [(lookup.get((r.foreground_event, region), ""), lookup.get((r.background_event, region), "")) for r in group.itertuples()]
            for token, model in models.items():
                methods = ("PWM_mean_odds",) if model["kind"] == "PWM" else ("finite_nonoverlap_counts", "overlapping_union_coverage")
                for width in (0, *sorted(set(widths))):
                    span = 1 if not width else max([max(0, min(len(a), len(b))-width+1) for a, b in seq_pairs], default=0)
                    # Keep an explicit NA row when no sequence can hold this window.
                    span = max(1, span)
                    for method in methods:
                        tests += span
                        cells += len(group)*span
                        if tests > max_tests:
                            raise ValueError("Explicit new-test budget exceeded; no models/windows are silently dropped")
                        if cells > max_matrix_cells:
                            raise ValueError("Explicit paired-matrix memory budget exceeded; no gene pairs silently dropped")
                        check_cancelled()
                        diff = np.full((len(group), span), np.nan)
                        for i, (a, b) in enumerate(seq_pairs):
                            if not a or not b:
                                continue
                            x, y = _profile(a, model, method, width), _profile(b, model, method, width)
                            if region.endswith("acceptor_intron") and width:
                                x, y = x[::-1], y[::-1]
                            n = min(len(x), len(y), span)
                            diff[i, :n] = x[:n]-y[:n]
                        datasets.append({"meta": {"event_type": kind, "direction": direction, "model_id": token,
                            "RBP": model.get("RBPs", ""), "region": region, "method": method, "window_nt": width,
                            "coordinate_system": "distance upstream from masked acceptor edge (-20nt), 3prime-aligned" if region.endswith("acceptor_intron") else "distance from masked donor edge (+6nt)" if region.endswith("donor_intron") else "affected RNA 5prime-relative"},
                            "pair_ids": group.pair_id.tolist(), "positions": list(range(span)), "differences": diff})
    return datasets, pairs.pair_id.tolist()
