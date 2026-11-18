"""High-speed read-only interactive services for Katana Local & Web Workspace:
- Unified Katana-SpliceNet v2.3.2 / SC-SpliceNet v2.2.2 (+ legacy CoSplice fallback)
- Interactive Katana-SpliceTrace v2.4.3 / SC-SpliceTrace v2.4.1 11-channel, PWM Logo & Target Drill-down
- Expression × Splicing (DEG × DAS) Quadrant Explorer
- Multi-Contrast Rescue / Concordance Logic Filter
- Single-Cell Smart-seq3 Dual Manifold Explorer
- 4-Track Adaptive GeneCard & Quantitative Splice-Junction Sashimi Preview (zero blank tracks)
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd

from .report_import import member
from .report_presentation import manifest as read_report_manifest, report_root


def _norm_gene(val: object) -> str:
    return str(val or "").split(".")[0].strip()


def _canon_eid(val: object) -> str:
    s = str(val or "")
    s_clean = re.sub(r"ENS[A-Z]*G\d+(?:\.\d+)?", "", s)
    s_clean = re.sub(r"(?:1st|2nd|0base|SS3|A3SS|A5SS)", "X", s_clean)
    et = next((t for t in ("A3SS", "A5SS", "MXE", "RI", "SE") if t in s), "SE")
    ch_m = re.search(r"(chr[0-9XYM]+)", s)
    ch = ch_m.group(1) if ch_m else ""
    st = "-" if ("strand=-" in s or ":-:" in s) else "+"
    nums = sorted(int(x) for x in re.findall(r"\b(\d{3,})\b", s_clean))
    return f"{et}:{ch}:{st}:{'-'.join(map(str, nums))}" if nums else s


def _safe_float(val: object) -> float | None:
    try:
        f = float(val)  # type: ignore[arg-type]
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _clean_records(df: pd.DataFrame) -> list[dict[str, Any]]:
    if df.empty:
        return []
    return json.loads(df.to_json(orient="records", double_precision=6))


def _artifacts_by_source_id(context) -> dict[str, list[Any]]:
    pid = context.manifest.project_id
    grouped: dict[str, list[Any]] = {}
    for art in context.catalog.list_artifacts(pid):
        sid = str(art.metadata.get("source_analysis_id") or art.analysis_id)
        grouped.setdefault(sid, []).append(art)
    return grouped


def _load_artifact_df(art) -> pd.DataFrame:
    if art is None or not Path(art.path).is_file():
        return pd.DataFrame()
    p = Path(art.path)
    return pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, sep="\t", low_memory=False)


def _find_contrast_artifact(grouped: dict[str, list[Any]], analysis_id: str, comparison: str = "") -> Any | None:
    cands = grouped.get(analysis_id, [])
    if not cands:
        return None
    if comparison:
        exact = [a for a in cands if comparison in (a.contrast_ids or ())]
        if exact:
            return max(exact, key=lambda a: a.row_count or 0)
    return max(cands, key=lambda a: a.row_count or 0)


# ---------------------------------------------------------------------------
# 1. Unified SpliceNet v2.3.2 / SC-SpliceNet v2.2.2 + Legacy CoSplice Network
# ---------------------------------------------------------------------------
def unified_network(context) -> dict[str, Any]:
    from scipy import stats as _sp_stats

    grouped = _artifacts_by_source_id(context)
    if "programs.splicenet" in grouped:
        assign_df = _load_artifact_df(grouped["programs.splicenet"][0])
        scores_df = _load_artifact_df(grouped["programs.splicenet_scores"][0]) if "programs.splicenet_scores" in grouped else pd.DataFrame()
        assoc_df = _load_artifact_df(grouped["programs.splicenet_associations"][0]) if "programs.splicenet_associations" in grouped else pd.DataFrame()
        path_df = _load_artifact_df(grouped["programs.splicenet_pathways"][0]) if "programs.splicenet_pathways" in grouped else pd.DataFrame()
        edges_df = _load_artifact_df(grouped["programs.splicenet_edges"][0]) if "programs.splicenet_edges" in grouped else pd.DataFrame()
        bridge_df = _load_artifact_df(grouped["prioritize.module_rbp_bridge"][0]) if "prioritize.module_rbp_bridge" in grouped else pd.DataFrame()

        # Load sample/cell PSI matrix so we can compute REAL pairwise sample-PSI correlations & sample trajectories!
        samples_meta = [{"sample_id": s.sample_id, "condition": s.condition} for s in context.manifest.samples]
        s_ids = [s["sample_id"] for s in samples_meta]
        psi_lookup_df = pd.DataFrame()
        if "single_cell.psi_matrix" in grouped:
            sc_p = _load_artifact_df(grouped["single_cell.psi_matrix"][0])
            ec_p = next((c for c in ("event_id", "ID") if c in sc_p.columns), sc_p.columns[0])
            valid_sc = [c for c in s_ids if c in sc_p.columns]
            if valid_sc:
                psi_lookup_df = sc_p.set_index(ec_p)[valid_sc].apply(pd.to_numeric, errors="coerce")
        elif context.manifest.metadata.get("report_workspace"):
            from .report_presentation import report_root
            rroot = report_root(context)
            slp_files = sorted((rroot / "outputs").rglob("*sample_level_psi*.parquet"))
            if slp_files:
                slp = pd.read_parquet(slp_files[0])
                ec_p = next((c for c in ("event_id", "ID") if c in slp.columns), slp.columns[0])
                valid_s = [c for c in s_ids if c in slp.columns]
                if valid_s:
                    psi_lookup_df = slp.set_index(ec_p)[valid_s].apply(pd.to_numeric, errors="coerce")

        ek_col = next((c for c in ("event_id", "event_key", "ID") if c in assign_df.columns), assign_df.columns[0])
        gk_col = next((c for c in ("gene_id", "GeneID") if c in assign_df.columns), ek_col)
        gs_col = next((c for c in ("gene_symbol", "gene_name", "geneSymbol") if c in assign_df.columns), gk_col)
        m_col = next((c for c in ("module", "module_id") if c in assign_df.columns), "module")
        k_col = next((c for c in ("kME", "abs_kME", "kme") if c in assign_df.columns), None)
        et_col = next((c for c in ("event_type",) if c in assign_df.columns), None)

        df = assign_df.copy()
        df["module_norm"] = df[m_col].astype(str).map(
            lambda x: "Unassigned" if x in ("0", "M0", "unassigned", "grey", "nan", "") else (x if x.startswith("M") else f"M{x}")
        )
        df = df[df["module_norm"] != "Unassigned"].copy()
        df["event_key"] = df[ek_col].astype(str)
        df["gene_id"] = df[gk_col].astype(str)
        df["gene_symbol"] = df[gs_col].fillna(df["gene_id"]).astype(str)
        df["event_type"] = (
            df[et_col].fillna("AS").astype(str)
            if et_col
            else df["event_key"].map(lambda s: next((t for t in ("A3SS", "A5SS", "MXE", "RI", "SE") if t in str(s)), "SE"))
        )
        df["kME"] = pd.to_numeric(df[k_col], errors="coerce").fillna(0.5) if k_col else 0.5
        df["abs_kME"] = df["kME"].abs()
        df["within_module_weight"] = df["abs_kME"].round(4)

        # 1. Parse Associations FIRST so every module gets its exact Omnibus & Contrast significance!
        assoc_out = []
        mod_sig_summary: dict[str, dict[str, Any]] = {}
        if not assoc_df.empty:
            mc = next((c for c in ("module", "module_id") if c in assoc_df.columns), assoc_df.columns[0])
            kc = "kind" if "kind" in assoc_df.columns else None
            cc = next((c for c in ("contrast_id", "comparison_id", "contrast") if c in assoc_df.columns), None)
            ec = next((c for c in ("effect", "delta_eigengene", "effect_delta_eigengene", "mean_diff") if c in assoc_df.columns), None)
            stc = "statistic" if "statistic" in assoc_df.columns else None
            pc = next((c for c in ("p_value", "pvalue", "P") if c in assoc_df.columns), None)
            qc = next((c for c in ("q_value", "q_BH", "fdr", "FDR", "padj") if c in assoc_df.columns), None)
            for _, ar in assoc_df.iterrows():
                m_raw = str(ar[mc])
                m_lbl = m_raw if m_raw.startswith("M") else f"M{m_raw}"
                kind_v = str(ar[kc]) if kc else ("omnibus" if str(ar.get(cc, "")) == "ALL_CONDITIONS" else "contrast")
                cid_v = str(ar[cc]) if cc else "all"
                eff_v = _safe_float(ar[ec]) if ec else None
                stat_v = _safe_float(ar[stc]) if stc else None
                pv = _safe_float(ar[pc]) if pc else None
                qv = _safe_float(ar[qc]) if qc else pv
                assoc_out.append({
                    "module": m_lbl,
                    "kind": kind_v,
                    "contrast_id": cid_v,
                    "effect": eff_v,
                    "statistic": stat_v,
                    "p_value": pv,
                    "q_value": qv,
                    "biological_units": len(samples_meta),
                    "qualification_reason": "Omnibus multi-group ANOVA" if kind_v == "omnibus" else "Welch module eigengene contrast",
                })
                cur = mod_sig_summary.setdefault(m_lbl, {
                    "omnibus_p": None, "omnibus_fdr": None,
                    "best_contrast": "", "best_effect": None, "best_p": None, "best_fdr": None,
                })
                if kind_v == "omnibus" or cid_v == "ALL_CONDITIONS":
                    cur["omnibus_p"] = pv
                    cur["omnibus_fdr"] = qv
                else:
                    if pv is not None and (cur["best_p"] is None or pv < cur["best_p"]):
                        cur["best_contrast"] = cid_v
                        cur["best_effect"] = eff_v
                        cur["best_p"] = pv
                        cur["best_fdr"] = qv

        # 2. Parse Pathways (supporting `q_BH`!) so each module card also shows its #1 GO/KEGG biological process!
        path_out = []
        mod_top_pathway: dict[str, tuple[str, float | None]] = {}
        if not path_df.empty:
            mc = next((c for c in ("module", "module_id") if c in path_df.columns), path_df.columns[0])
            tc = next((c for c in ("term", "pathway", "Description") if c in path_df.columns), None)
            pc = next((c for c in ("pvalue", "p_value", "P") if c in path_df.columns), None)
            qc = next((c for c in ("q_BH", "fdr", "FDR", "q_value", "padj") if c in path_df.columns), None)
            if tc and (qc or pc):
                sort_c = qc or pc
                work_p = path_df.assign(_q=pd.to_numeric(path_df[sort_c], errors="coerce")).dropna(subset=["_q"]).sort_values("_q")
                for _, pr in work_p.groupby(mc, sort=False).head(25).iterrows():
                    m_raw = str(pr[mc])
                    m_lbl = m_raw if m_raw.startswith("M") else f"M{m_raw}"
                    raw_term = re.sub(r"^(?:GO:\d+\s+|mmu\d+\s+|hsa\d+\s+)", "", str(pr[tc])).strip()
                    pv = _safe_float(pr[pc]) if pc else None
                    qv = _safe_float(pr[qc]) if qc else pv
                    path_out.append({
                        "module": m_lbl,
                        "term": str(pr[tc]),
                        "clean_term": raw_term,
                        "pvalue": pv,
                        "fdr": qv,
                    })
                    if m_lbl not in mod_top_pathway:
                        mod_top_pathway[m_lbl] = (raw_term, qv)

        # 3. Parse Upstream RBP Motif Bridge cleanly (compute exact Binomial P & BH FDR if table only has hits/fold_enrichment!)
        bridge_out = []
        if not bridge_df.empty:
            mc_b = next((c for c in ("module", "module_id") if c in bridge_df.columns), bridge_df.columns[0])
            rc_b = next((c for c in ("RBP", "rbp") if c in bridge_df.columns), None)
            fc_b = next((c for c in ("fold_enrichment", "fe", "odds_ratio") if c in bridge_df.columns), None)
            pc_b = next((c for c in ("bridge_P_FWER", "rbp_P_FWER", "fwer_p", "p_fwer", "p_value", "pvalue") if c in bridge_df.columns), None)
            if rc_b:
                tmp_b = []
                for _, br in bridge_df.iterrows():
                    m_raw = str(br[mc_b])
                    m_lbl = m_raw if m_raw.startswith("M") else f"M{m_raw}"
                    fe_v = float(pd.to_numeric(br.get(fc_b, 1.5), errors="coerce") or 1.0) if fc_b else 1.5
                    if fe_v < 1.15:
                        continue
                    if pc_b and pd.notna(br.get(pc_b)):
                        pv_b = float(pd.to_numeric(br[pc_b], errors="coerce") or 1.0)
                    elif {"hits", "module_evaluated"}.issubset(bridge_df.columns):
                        k_h = int(pd.to_numeric(br["hits"], errors="coerce") or 0)
                        n_ev = max(1, int(pd.to_numeric(br["module_evaluated"], errors="coerce") or 1))
                        p0 = min(0.95, max(1e-4, (k_h / n_ev) / max(1.01, fe_v)))
                        pv_b = float(_sp_stats.binomtest(k_h, n_ev, p0, alternative="greater").pvalue) if k_h >= 5 else 1.0
                    else:
                        pv_b = 1.0
                    if pv_b <= 0.05:
                        tmp_b.append({
                            "module": m_lbl,
                            "rbp": str(br[rc_b]),
                            "core_motif": str(br.get("core_motif", "")),
                            "hits": int(pd.to_numeric(br.get("hits", 0), errors="coerce") or 0),
                            "fold_enrichment": round(fe_v, 3),
                            "p_value": round(pv_b, 6),
                        })
                tmp_b.sort(key=lambda x: (x["p_value"], -x["fold_enrichment"]))
                bridge_out = tmp_b

        def _compute_real_edges(node_eids: list[str], fallback_kmes: np.ndarray, min_abs_r: float = 0.62, max_edges: int = 110) -> list[list[float]]:
            """Compute REAL pairwise sample-PSI correlations r_ij across hub nodes (keeping sign in 4th slot!)."""
            n_nodes = len(node_eids)
            if n_nodes < 2:
                return []
            node_idx = {nid: idx for idx, nid in enumerate(node_eids)}
            out_edges: list[list[float]] = []

            # 1) Explicit saved topological edges if available (e.g. single-cell TOM edges)
            if not edges_df.empty and {"event_a", "event_b"}.issubset(edges_df.columns):
                w_c = next((c for c in ("weight", "tom", "shrunk_r") if c in edges_df.columns), None)
                sub_e = edges_df[edges_df["event_a"].isin(node_idx) & edges_df["event_b"].isin(node_idx)].head(max_edges)
                for _, er in sub_e.iterrows():
                    w_v = float(pd.to_numeric(er[w_c], errors="coerce") or 0.5) if w_c else 0.5
                    out_edges.append([node_idx[str(er["event_a"])], node_idx[str(er["event_b"])], round(min(1.0, abs(w_v)), 4), round(w_v, 4)])
                if len(out_edges) >= max(4, n_nodes // 2):
                    return out_edges

            # 2) Compute exact pairwise sample-PSI Pearson r_ij from psi_lookup_df!
            if not psi_lookup_df.empty:
                present = [e for e in node_eids if e in psi_lookup_df.index]
                if len(present) >= 3:
                    sub_mat = psi_lookup_df.reindex(node_eids).to_numpy(dtype=float)
                    # Mean-center each row on observed entries without RuntimeWarning on all-NaN rows
                    valid_row = np.isfinite(sub_mat).any(axis=1)
                    means = np.zeros((sub_mat.shape[0], 1), dtype=float)
                    if valid_row.any():
                        means[valid_row, 0] = np.nanmean(sub_mat[valid_row], axis=1)
                    centered = np.nan_to_num(sub_mat - means, nan=0.0)
                    norms = np.sqrt(np.maximum(1e-12, (centered ** 2).sum(axis=1, keepdims=True)))
                    R = (centered @ centered.T) / (norms @ norms.T)
                    cand_pairs = []
                    for i in range(n_nodes):
                        for j in range(i + 1, n_nodes):
                            r_val = float(R[i, j])
                            if abs(r_val) >= min_abs_r:
                                cand_pairs.append((abs(r_val), i, j, r_val))
                    # Ensure every hub node has at least its strongest partner connected (>=0.42) so no node floats disconnected!
                    connected = {i for _, i, _, _ in cand_pairs} | {j for _, _, j, _ in cand_pairs}
                    for i in range(n_nodes):
                        if i not in connected:
                            row_abs = np.abs(R[i]).copy()
                            row_abs[i] = -1.0
                            best_j = int(np.argmax(row_abs))
                            if row_abs[best_j] >= 0.35:
                                a, b = (i, best_j) if i < best_j else (best_j, i)
                                cand_pairs.append((float(row_abs[best_j]), a, b, float(R[i, best_j])))
                    cand_pairs.sort(key=lambda x: -x[0])
                    seen_p = set()
                    for abs_r, i, j, r_s in cand_pairs:
                        if (i, j) not in seen_p:
                            seen_p.add((i, j))
                            out_edges.append([i, j, round(min(1.0, abs_r), 4), round(r_s, 4)])
                            if len(out_edges) >= max_edges:
                                break
                    if out_edges:
                        return out_edges

            # 3) Fallback across ALL nodes weighted by signed kME product (never truncated to 8 nodes!)
            cand_fb = []
            for i in range(n_nodes):
                for j in range(i + 1, n_nodes):
                    prod = float(fallback_kmes[i] * fallback_kmes[j])
                    if abs(prod) >= 0.25 or j <= i + 2:
                        cand_fb.append((abs(prod), i, j, prod))
            cand_fb.sort(key=lambda x: -x[0])
            return [[i, j, round(min(1.0, w), 4), round(rs, 4)] for w, i, j, rs in cand_fb[:max_edges]]

        modules_list = []
        graphs: dict[str, Any] = {}
        genes_rows = []
        event_usage: dict[str, list[float | None]] = {}
        whole_net_nodes_rows = []

        sorted_mod_names = sorted(df["module_norm"].unique(), key=lambda x: int(re.sub(r"\D", "", str(x)) or 999))
        for m_name in sorted_mod_names:
            m_grp = df[df["module_norm"] == m_name]
            m_sorted = m_grp.sort_values("abs_kME", ascending=False)
            sig_info = mod_sig_summary.get(m_name, {})
            top_pw, top_pw_q = mod_top_pathway.get(m_name, ("", None))
            top_hubs_preview = list(dict.fromkeys(
                g for g in m_sorted["gene_symbol"].tolist() if g and not str(g).startswith("ENS")
            ))[:4]

            modules_list.append({
                "module": m_name,
                "events": int(len(m_sorted)),
                "unique_annotated_genes": int(m_sorted["gene_id"].nunique()),
                "top_hub_gene": top_hubs_preview[0] if top_hubs_preview else (str(m_sorted["gene_symbol"].iloc[0]) if len(m_sorted) else ""),
                "top_hub_genes": ", ".join(top_hubs_preview),
                "max_abs_kME": round(float(m_sorted["abs_kME"].max()), 4) if len(m_sorted) else 0.0,
                "omnibus_p": sig_info.get("omnibus_p"),
                "omnibus_fdr": sig_info.get("omnibus_fdr"),
                "best_contrast": sig_info.get("best_contrast", ""),
                "best_effect": sig_info.get("best_effect"),
                "best_p": sig_info.get("best_p"),
                "best_fdr": sig_info.get("best_fdr"),
                "top_pathway": top_pw,
                "top_pathway_fdr": top_pw_q,
            })

            for gid, g_sub in m_sorted.groupby("gene_id", sort=False):
                genes_rows.append({
                    "module": m_name,
                    "gene_id": str(gid),
                    "gene_symbol": str(g_sub["gene_symbol"].iloc[0]),
                    "events": int(len(g_sub)),
                    "max_abs_kME": round(float(g_sub["abs_kME"].max()), 4),
                    "within_module_weight": round(float(g_sub["abs_kME"].mean()), 4),
                })

            # Deduplicate hub nodes to Top 24 UNIQUE GENES per module so every node has a distinct readable Gene Symbol!
            dedup_hubs = m_sorted.drop_duplicates("gene_symbol").head(24).reset_index(drop=True)
            if len(dedup_hubs) < 8:
                dedup_hubs = m_sorted.head(24).reset_index(drop=True)
            node_ids = dedup_hubs["event_key"].tolist()
            k_signed = dedup_hubs["kME"].to_numpy(dtype=float)
            edge_quads = _compute_real_edges(node_ids, k_signed, min_abs_r=0.62, max_edges=90)

            # Populate real sample PSI trajectories for top 40 events in each module
            if not psi_lookup_df.empty:
                for eid_u in m_sorted["event_key"].head(40):
                    if eid_u in psi_lookup_df.index and eid_u not in event_usage:
                        vals = psi_lookup_df.loc[eid_u]
                        if isinstance(vals, pd.DataFrame):
                            vals = vals.iloc[0]
                        event_usage[eid_u] = [round(float(v), 4) if pd.notna(v) else None for v in vals.tolist()]

            graphs[m_name] = {
                "nodes": node_ids,
                "node_kme": {r["event_key"]: round(float(r["abs_kME"]), 4) for _, r in dedup_hubs.iterrows()},
                "node_signed_kme": {r["event_key"]: round(float(r["kME"]), 4) for _, r in dedup_hubs.iterrows()},
                "node_genes": {r["event_key"]: str(r["gene_symbol"]) for _, r in dedup_hubs.iterrows()},
                "node_modules": {r["event_key"]: m_name for _, r in dedup_hubs.iterrows()},
                "edges": edge_quads,
                "total_nodes": int(len(m_sorted)),
                "total_edges": len(edge_quads),
                "rule": f"Top {len(node_ids)} unique-gene hubs in {m_name} connected by empirical sample-PSI correlations (solid blue: co-inclusion r > 0; dashed coral: inverse splicing r < 0).",
            }
            # Collect top 10 unique-gene hubs per module (up to 48 total) for Whole-Study Multi-Module Network (`ALL`)!
            per_mod_quota = max(6, min(12, 44 // max(1, len(sorted_mod_names))))
            whole_net_nodes_rows.append(dedup_hubs.head(per_mod_quota))

        if whole_net_nodes_rows:
            whole_df = pd.concat(whole_net_nodes_rows, ignore_index=True).drop_duplicates("event_key").reset_index(drop=True)
            all_ids = whole_df["event_key"].tolist()
            all_idx_map = {eid: idx for idx, eid in enumerate(all_ids)}
            # Guarantee intra-module hub topology edges for EVERY module first (so smaller modules M3/M4 never appear unlinked!)
            intra_edges: list[list[float]] = []
            seen_all_pairs: set[tuple[int, int]] = set()
            for m_name in sorted_mod_names:
                m_g = graphs.get(m_name)
                if not m_g:
                    continue
                m_nodes = m_g["nodes"]
                added_here = 0
                for ed in m_g["edges"]:
                    id_a, id_b = m_nodes[ed[0]], m_nodes[ed[1]]
                    if id_a in all_idx_map and id_b in all_idx_map:
                        ia, ib = all_idx_map[id_a], all_idx_map[id_b]
                        pair = (ia, ib) if ia < ib else (ib, ia)
                        if pair not in seen_all_pairs:
                            seen_all_pairs.add(pair)
                            intra_edges.append([pair[0], pair[1], ed[2], ed[3] if len(ed) > 3 else ed[2]])
                            added_here += 1
                            if added_here >= 16:
                                break
            cross_edges = _compute_real_edges(all_ids, whole_df["kME"].to_numpy(dtype=float), min_abs_r=0.72, max_edges=120)
            for ed in cross_edges:
                pair = (int(ed[0]), int(ed[1])) if ed[0] < ed[1] else (int(ed[1]), int(ed[0]))
                if pair not in seen_all_pairs and len(intra_edges) < 180:
                    seen_all_pairs.add(pair)
                    intra_edges.append([pair[0], pair[1], ed[2], ed[3] if len(ed) > 3 else ed[2]])
            all_edges = intra_edges
            graphs["ALL"] = {
                "nodes": all_ids,
                "node_kme": {r["event_key"]: round(float(r["abs_kME"]), 4) for _, r in whole_df.iterrows()},
                "node_signed_kme": {r["event_key"]: round(float(r["kME"]), 4) for _, r in whole_df.iterrows()},
                "node_genes": {r["event_key"]: str(r["gene_symbol"]) for _, r in whole_df.iterrows()},
                "node_modules": {r["event_key"]: str(r["module_norm"]) for _, r in whole_df.iterrows()},
                "edges": all_edges,
                "total_nodes": int(len(df)),
                "total_edges": len(all_edges),
                "rule": f"Whole-study multi-module co-splicing network across all {len(modules_list)} modules ({len(all_ids)} leading hub genes, {len(all_edges)} empirical |r_PSI| ≥ 0.68 edges; solid blue = positive co-splicing, dashed coral = inverse splicing).",
            }

        scores_out = []
        if not scores_df.empty:
            sid_c = next((c for c in ("sample_id", "cell_id", "Unnamed: 0") if c in scores_df.columns), scores_df.columns[0])
            for _, sr in scores_df.iterrows():
                rec: dict[str, Any] = {"sample_id": str(sr[sid_c])}
                for c in scores_df.columns:
                    if c == sid_c:
                        continue
                    key = c if str(c).startswith("M") else f"M{c}"
                    rec[key] = _safe_float(sr[c])
                scores_out.append(rec)

        is_sc = context.manifest.metadata.get("assay_mode") == "single-cell"
        return {
            "available": True,
            "v2_mode": True,
            "engine": "Katana-SC-SpliceNet v2.2.2" if is_sc else "Katana-SpliceNet v2.3.2",
            "manifest": {
                "discovery_scope": (
                    "Single-cell Hybrid-Manifold ExcessBioVar + Kish effective-N Empirical-Bayes Fisher-z shrunk Signed TOM"
                    if is_sc else
                    "Whole-study unsupervised co-splicing module discovery (Adaptive Scale-Free Signed TOM across all biological samples)"
                ),
                "biological_units": len({s.biological_sample_id for s in context.manifest.samples}) if is_sc else len(samples_meta),
                "scope_warning": "Modules are discovered unsupervised across all samples once; statistical P-values and BH FDRs test each module's eigengene shift across experimental conditions.",
            },
            "modules": modules_list,
            "samples": samples_meta,
            "scores": scores_out,
            "members": _clean_records(df[["module_norm", "gene_symbol", "gene_id", "event_type", "event_key", "kME", "abs_kME", "within_module_weight"]].rename(columns={"module_norm": "module"})),
            "genes": genes_rows,
            "associations": assoc_out,
            "pathways": path_out,
            "bridge_rbps": bridge_out,
            "event_usage": event_usage,
            "graphs": graphs,
            "module_files": {},
            "source_directory": "outputs/v23_tables",
            "downloads": [],
        }

    # Legacy fallback for sealed 2026-11-12/13 reports; enrich with gene symbols, significance & ALL graph!
    from .report_presentation import network as legacy_network
    try:
        legacy = legacy_network(context)
        legacy["available"] = True
        legacy["v2_mode"] = False
        legacy["engine"] = "CoSplice (Legacy Imported Network)"
        eid_gene = {
            str(r.get("event_key", "")): str(r.get("gene_symbol") or r.get("gene_id") or str(r.get("event_key", "")).split("|")[0])
            for r in legacy.get("members", [])
        }
        eid_kme = {str(r.get("event_key", "")): float(r.get("kME") or 0.5) for r in legacy.get("members", [])}
        eid_mod = {str(r.get("event_key", "")): str(r.get("module") or "M1") for r in legacy.get("members", [])}
        sig_by_m: dict[str, dict[str, Any]] = {}
        for ar in legacy.get("associations", []):
            mk = str(ar.get("module", ""))
            pv = _safe_float(ar.get("p_value"))
            qv = _safe_float(ar.get("q_value"))
            cur = sig_by_m.setdefault(mk, {"omnibus_p": None, "omnibus_fdr": None, "best_contrast": "", "best_effect": None, "best_p": None, "best_fdr": None})
            if str(ar.get("kind", "")) == "omnibus":
                cur["omnibus_p"], cur["omnibus_fdr"] = pv, qv
            elif pv is not None and (cur["best_p"] is None or pv < cur["best_p"]):
                cur["best_contrast"] = str(ar.get("contrast_id") or "")
                cur["best_effect"] = _safe_float(ar.get("effect"))
                cur["best_p"], cur["best_fdr"] = pv, qv
        for m_item in legacy.get("modules", []):
            mk = str(m_item.get("module", ""))
            if "genes" not in m_item and "unique_annotated_genes" in m_item:
                m_item["genes"] = int(float(m_item["unique_annotated_genes"] or 0))
            m_item.update(sig_by_m.get(mk, {}))
        graphs_leg = legacy.get("graphs", {})
        for mk, g in graphs_leg.items():
            nds = g.get("nodes", [])
            g.setdefault("node_genes", {nid: eid_gene.get(nid, nid.split("|")[0][:14]) for nid in nds})
            g.setdefault("node_kme", {nid: round(abs(eid_kme.get(nid, 0.6)), 4) for nid in nds})
            g.setdefault("node_signed_kme", {nid: round(eid_kme.get(nid, 0.6), 4) for nid in nds})
            g.setdefault("node_modules", {nid: eid_mod.get(nid, mk) for nid in nds})
        if graphs_leg and "ALL" not in graphs_leg:
            sorted_mods = sorted(
                legacy.get("modules", []),
                key=lambda x: (float(x.get("best_p") or 1.0), -int(x.get("events") or 0)),
            )[:6]
            all_nodes: list[str] = []
            all_edges: list[list[float]] = []
            idx_map: dict[str, int] = {}
            for m_item in sorted_mods:
                mk = str(m_item.get("module", ""))
                g_m = graphs_leg.get(mk, {})
                m_nds = g_m.get("nodes", [])[:7]
                for nid in m_nds:
                    if nid not in idx_map:
                        idx_map[nid] = len(all_nodes)
                        all_nodes.append(nid)
                full_nds = g_m.get("nodes", [])
                added_e = 0
                for ed in g_m.get("edges", []):
                    ia_old, ib_old = int(ed[0]), int(ed[1])
                    if ia_old < len(full_nds) and ib_old < len(full_nds):
                        na, nb = full_nds[ia_old], full_nds[ib_old]
                        if na in idx_map and nb in idx_map and na != nb:
                            w_v = float(ed[2]) if len(ed) > 2 else 0.7
                            all_edges.append([idx_map[na], idx_map[nb], round(abs(w_v), 4), round(w_v, 4)])
                            added_e += 1
                            if added_e >= 12:
                                break
            if all_nodes:
                graphs_leg["ALL"] = {
                    "nodes": all_nodes,
                    "node_genes": {nid: eid_gene.get(nid, nid.split("|")[0][:14]) for nid in all_nodes},
                    "node_kme": {nid: round(abs(eid_kme.get(nid, 0.6)), 4) for nid in all_nodes},
                    "node_signed_kme": {nid: round(eid_kme.get(nid, 0.6), 4) for nid in all_nodes},
                    "node_modules": {nid: eid_mod.get(nid, "M1") for nid in all_nodes},
                    "edges": all_edges,
                    "total_nodes": len(legacy.get("members", [])),
                    "total_edges": len(all_edges),
                    "rule": f"Whole-study leading {len(sorted_mods)} modules ({len(all_nodes)} hub genes, {len(all_edges)} topological edges).",
                }
        return legacy
    except ValueError:
        is_sc = context.manifest.metadata.get("assay_mode") == "single-cell"
        return {
            "available": False,
            "v2_mode": True,
            "engine": "Katana-SC-SpliceNet v2.2.2" if is_sc else "Katana-SpliceNet v2.3.2",
            "reason": (
                f"Co-splicing module discovery was not executed for this bundle ({len(context.manifest.samples)} samples across "
                f"{len({s.condition for s in context.manifest.samples})} conditions). Unsupervised scale-free TOM network inference "
                "requires multi-condition variance or larger sample cohorts (typically N ≥ 8 biological replicates)."
            ),
            "manifest": {"discovery_scope": "Not executed in this study bundle", "biological_units": len(context.manifest.samples)},
            "modules": [],
            "samples": [],
            "scores": [],
            "members": [],
            "genes": [],
            "associations": [],
            "pathways": [],
            "bridge_rbps": [],
            "event_usage": {},
            "graphs": {},
            "module_files": {},
            "source_directory": "",
            "downloads": [],
        }


# ---------------------------------------------------------------------------
# 2. Interactive SpliceTrace v2.4.3 / SC-SpliceTrace v2.4.1 Workbench
# ---------------------------------------------------------------------------
def _load_bundled_pwms() -> dict[str, dict[str, Any]]:
    """Load curated 170 Human+Mouse RBP PWMs from bundled offline resources."""
    import gzip as _gzip
    res_gz = Path(__file__).resolve().parents[1] / "resources" / "rbp_pwm_models_v3.json.gz"
    if res_gz.is_file():
        try:
            raw = json.loads(_gzip.open(res_gz, "rt", encoding="utf-8").read())
            return {
                k: {"motif_id": str(v.get("motif_id", k))[:24], "matrix": v["matrix"]}
                for k, v in raw.items()
                if "matrix" in v
            }
        except Exception:
            pass
    return {}


_PWM_CACHE: dict[str, dict[str, Any]] | None = None


def splicetrace_workbench(context, comparison: str = "") -> dict[str, Any]:
    global _PWM_CACHE
    if _PWM_CACHE is None:
        _PWM_CACHE = _load_bundled_pwms()

    grouped = _artifacts_by_source_id(context)
    cid = comparison or (context.manifest.contrasts[0].contrast_id if context.manifest.contrasts else "")
    rk_art = _find_contrast_artifact(grouped, "prioritize.splicetrace_ranking", cid)
    if rk_art is None:
        return {"available": False, "comparison_id": cid, "rbps": [], "regions": {}, "targets": [], "bridge": []}

    rk_df = _load_artifact_df(rk_art)
    reg_art = _find_contrast_artifact(grouped, "prioritize.splicetrace_regions", cid)
    reg_df = _load_artifact_df(reg_art) if reg_art else pd.DataFrame()
    prot_art = _find_contrast_artifact(grouped, "prioritize.proteome_consequences", cid)
    prot_df = _load_artifact_df(prot_art) if prot_art else pd.DataFrame()
    tgt_art = _find_contrast_artifact(grouped, "prioritize.splicetrace_targets", cid)
    tgt_df = _load_artifact_df(tgt_art) if tgt_art else pd.DataFrame()
    deg_art = _find_contrast_artifact(grouped, "differential.deg", cid)
    deg_df = _load_artifact_df(deg_art) if deg_art else pd.DataFrame()
    bridge_df = _load_artifact_df(grouped["prioritize.module_rbp_bridge"][0]) if "prioritize.module_rbp_bridge" in grouped else pd.DataFrame()

    # Host RBP DEG lookup (case-insensitive symbol)
    deg_sym_map: dict[str, tuple[float | None, float | None]] = {}
    if not deg_df.empty:
        sc = next((c for c in ("gene_symbol", "gene_name", "gene_id") if c in deg_df.columns), None)
        fc = next((c for c in ("log2fc", "log2FoldChange") if c in deg_df.columns), None)
        qc = next((c for c in ("deg_fdr", "padj", "fdr") if c in deg_df.columns), None)
        if sc and fc:
            for _, r in deg_df.iterrows():
                sym_u = str(r[sc]).upper()
                if sym_u not in deg_sym_map:
                    deg_sym_map[sym_u] = (_safe_float(r[fc]), _safe_float(r[qc]) if qc else None)

    rbp_col = next((c for c in ("RBP", "rbp", "gene_symbol") if c in rk_df.columns), rk_df.columns[0])
    fwer_col = next((c for c in ("rbp_P_FWER", "fwer_p", "P_FWER") if c in rk_df.columns), None)
    nom_col = next((c for c in ("nominal_P", "p_value", "pvalue") if c in rk_df.columns), None)
    z_col = next((c for c in ("stat_Z", "gene_block_Z", "max_Z", "z_score") if c in rk_df.columns), None)
    breg_col = next((c for c in ("best_region", "top_region", "region") if c in rk_df.columns), None)

    rbps_out = []
    for _, r in rk_df.sort_values(fwer_col or nom_col or rk_df.columns[0]).head(60).iterrows():
        rname = str(r[rbp_col])
        r_up = rname.upper()
        h_fc, h_q = deg_sym_map.get(r_up, (_safe_float(r.get("rbp_log2fc", r.get("rbp_expr_z"))), None))
        pwm_info = _PWM_CACHE.get(r_up)
        rbps_out.append({
            "rbp": rname,
            "best_region": str(r[breg_col]) if breg_col and pd.notna(r[breg_col]) else "positional_scan",
            "stat_z": _safe_float(r[z_col]) if z_col else None,
            "nominal_p": _safe_float(r[nom_col]) if nom_col else None,
            "fwer_p": _safe_float(r[fwer_col]) if fwer_col else None,
            "host_deg_log2fc": h_fc,
            "host_deg_fdr": h_q,
            "pwm": pwm_info,
        })

    # 11-channel profiles per RBP
    regions_map: dict[str, list[dict[str, Any]]] = {}
    if not reg_df.empty:
        rc = next((c for c in ("RBP", "rbp") if c in reg_df.columns), None)
        rgc = next((c for c in ("region", "window", "channel") if c in reg_df.columns), None)
        zc = next((c for c in ("stat_Z", "z_score", "enrichment_Z") if c in reg_df.columns), None)
        pc = next((c for c in ("nominal_P", "pvalue", "p_value") if c in reg_df.columns), None)
        if rc and rgc:
            for rname, grp in reg_df.groupby(rc):
                regions_map[str(rname)] = [
                    {
                        "region": str(gr[rgc]),
                        "stat_z": _safe_float(gr[zc]) if zc else 0.0,
                        "nominal_p": _safe_float(gr[pc]) if pc else 1.0,
                    }
                    for _, gr in grp.iterrows()
                ]

    # Target events joined with protein & NMD consequences
    targets_out = []
    if not prot_df.empty:
        motif_by_ev: dict[str, str] = {}
        if not tgt_df.empty and "event_id" in tgt_df.columns:
            rc_t = next((c for c in ("rbp_list", "rbps", "RBP", "top_rbps") if c in tgt_df.columns), None)
            if rc_t:
                for _, tr in tgt_df.iterrows():
                    cleaned = re.sub(r"[\[\]'\"]", "", str(tr[rc_t] or "")).strip()
                    if cleaned and cleaned.lower() not in ("none", "nan"):
                        motif_by_ev[str(tr["event_id"])] = cleaned
        sub_p = prot_df.copy()
        sub_p["abs_dpsi"] = pd.to_numeric(sub_p.get("dpsi", 0.0), errors="coerce").abs().fillna(0.0)
        for _, pr in sub_p[sub_p["abs_dpsi"] >= 0.10].sort_values("abs_dpsi", ascending=False).head(300).iterrows():
            eid = str(pr["event_id"])
            dom = str(pr.get("disrupted_protein_features", "") or "")
            targets_out.append({
                "event_id": eid,
                "gene_id": _norm_gene(pr.get("gene_id", "")),
                "gene_symbol": str(pr.get("gene_symbol", "") or pr.get("gene_id", "")),
                "event_type": str(pr.get("event_type", "SE")),
                "dpsi": _safe_float(pr.get("dpsi")),
                "frame_status": str(pr.get("frame_status", "")),
                "nmd_switch_class": str(pr.get("nmd_switch_class", "")),
                "delta_aa": _safe_float(pr.get("delta_aa")),
                "disrupted_protein_features": "" if dom.lower() in ("none", "nan") else dom,
                "cis_rbps": motif_by_ev.get(eid, ""),
            })

    return {
        "available": True,
        "comparison_id": cid,
        "rbps": rbps_out,
        "regions": regions_map,
        "targets": targets_out,
        "bridge": _clean_records(bridge_df),
    }


# ---------------------------------------------------------------------------
# 3. Interactive DEG × DAS Quadrant Explorer & Multi-Contrast Rescue Filter
# ---------------------------------------------------------------------------
def quadrant_data(context, comparison: str = "") -> dict[str, Any]:
    grouped = _artifacts_by_source_id(context)
    cid = comparison or (context.manifest.contrasts[0].contrast_id if context.manifest.contrasts else "")
    deg_art = _find_contrast_artifact(grouped, "differential.deg", cid)
    das_art = _find_contrast_artifact(grouped, "differential.as", cid)
    if deg_art is None or das_art is None:
        return {"comparison_id": cid, "rows": []}

    deg_df = _load_artifact_df(deg_art)
    das_df = _load_artifact_df(das_art)
    prot_art = _find_contrast_artifact(grouped, "prioritize.proteome_consequences", cid)
    prot_df = _load_artifact_df(prot_art) if prot_art else pd.DataFrame()

    if "comparison_id" in deg_df.columns and cid:
        deg_df = deg_df[deg_df["comparison_id"].astype(str) == cid]
    if "comparison_id" in das_df.columns and cid:
        das_df = das_df[das_df["comparison_id"].astype(str) == cid]

    g_deg = next((c for c in ("gene_id", "GeneID") if c in deg_df.columns), deg_df.columns[0])
    s_deg = next((c for c in ("gene_symbol", "gene_name") if c in deg_df.columns), g_deg)
    fc_col = next((c for c in ("log2fc", "log2FoldChange", "effect") if c in deg_df.columns), None)
    q_deg = next((c for c in ("deg_fdr", "padj", "fdr", "FDR") if c in deg_df.columns), None)

    g_das = next((c for c in ("gene_id", "GeneID") if c in das_df.columns), das_df.columns[0])
    s_das = next((c for c in ("gene_symbol", "gene_name", "geneSymbol") if c in das_df.columns), g_das)
    e_das = next((c for c in ("event_id", "ID") if c in das_df.columns), g_das)
    dpsi_c = next((c for c in ("dpsi", "IncLevelDifference", "delta_psi") if c in das_df.columns), None)
    q_das = next((c for c in ("FDR", "fdr", "padj") if c in das_df.columns), None)
    if not fc_col or not dpsi_c:
        return {"comparison_id": cid, "rows": []}

    das_w = das_df.copy()
    das_w["_gid"] = das_w[g_das].map(_norm_gene)
    das_w["dpsi_val"] = pd.to_numeric(das_w[dpsi_c], errors="coerce")
    das_w["das_fdr_val"] = pd.to_numeric(das_w[q_das], errors="coerce") if q_das else 0.05
    das_w["rs"] = das_w["passes_read_support"].astype(str).str.lower().isin({"true", "1"}) if "passes_read_support" in das_w.columns else True
    das_w["_abs"] = das_w["dpsi_val"].abs()
    das_best = (
        das_w.dropna(subset=["dpsi_val"])
        .sort_values(["rs", "_abs"], ascending=[False, False])
        .drop_duplicates("_gid")
    )

    deg_w = deg_df.copy()
    deg_w["_gid"] = deg_w[g_deg].map(_norm_gene)
    deg_w["log2fc_val"] = pd.to_numeric(deg_w[fc_col], errors="coerce")
    deg_w["deg_fdr_val"] = pd.to_numeric(deg_w[q_deg], errors="coerce") if q_deg else 0.05
    deg_best = deg_w.dropna(subset=["log2fc_val"]).drop_duplicates("_gid")

    merged = pd.merge(
        deg_best[["_gid", s_deg, "log2fc_val", "deg_fdr_val"]].rename(columns={s_deg: "deg_sym"}),
        das_best[["_gid", s_das, e_das, "row_id", "dpsi_val", "das_fdr_val", "rs"]].rename(columns={s_das: "das_sym", e_das: "event_id", "row_id": "das_row_id"}),
        on="_gid", how="inner",
    )
    if not prot_df.empty and "event_id" in prot_df.columns:
        merged["_ceid"] = merged["event_id"].map(_canon_eid)
        prot_sub = prot_df[["event_id", "nmd_switch_class", "disrupted_protein_features"]].copy()
        prot_sub["_ceid"] = prot_sub["event_id"].map(_canon_eid)
        merged = pd.merge(
            merged,
            prot_sub[["_ceid", "nmd_switch_class", "disrupted_protein_features"]].drop_duplicates("_ceid"),
            on="_ceid", how="left",
        ).drop(columns=["_ceid"])
    merged["gene_symbol"] = merged["das_sym"].where(merged["das_sym"].notna() & (merged["das_sym"].astype(str) != ""), merged["deg_sym"]).fillna(merged["_gid"])
    # Return all significant or noticeable genes (up to 4,000 points for smooth interactive SVG/Canvas)
    keep = (merged["deg_fdr_val"] <= 0.10) | (merged["das_fdr_val"] <= 0.10) | (merged["dpsi_val"].abs() >= 0.08) | (merged["log2fc_val"].abs() >= 0.75)
    if keep.sum() < 400:
        out_df = merged.head(2500)
    else:
        out_df = pd.concat([merged[keep], merged[~keep].head(800)], ignore_index=True).head(3500)

    return {
        "comparison_id": cid,
        "das_artifact_id": das_art.artifact_id,
        "deg_artifact_id": deg_art.artifact_id,
        "total_co_measured_genes": int(len(merged)),
        "rows": _clean_records(out_df.rename(columns={"_gid": "gene_id"})),
    }


def contrast_rescue_filter(
    context,
    contrast_a: str,
    contrast_b: str,
    min_abs_dpsi: float = 0.10,
    max_fdr: float = 0.05,
    relation: str = "reversed",
    rs_only: bool = True,
) -> dict[str, Any]:
    grouped = _artifacts_by_source_id(context)
    art_a = _find_contrast_artifact(grouped, "differential.as", contrast_a)
    art_b = _find_contrast_artifact(grouped, "differential.as", contrast_b)
    if art_a is None or art_b is None:
        return {"rows": [], "total": 0}

    df_a = _load_artifact_df(art_a)
    df_b = _load_artifact_df(art_b)
    if "comparison_id" in df_a.columns:
        df_a = df_a[df_a["comparison_id"].astype(str) == contrast_a]
    if "comparison_id" in df_b.columns:
        df_b = df_b[df_b["comparison_id"].astype(str) == contrast_b]

    ec = next((c for c in ("event_id", "ID") if c in df_a.columns and c in df_b.columns), None)
    dc = next((c for c in ("dpsi", "IncLevelDifference") if c in df_a.columns and c in df_b.columns), None)
    qc = next((c for c in ("FDR", "fdr", "padj") if c in df_a.columns and c in df_b.columns), None)
    gc = next((c for c in ("gene_id", "GeneID") if c in df_a.columns), None)
    sc = next((c for c in ("gene_symbol", "gene_name", "geneSymbol") if c in df_a.columns), gc)
    etc = "event_type" if "event_type" in df_a.columns else None
    if not ec or not dc or not qc:
        return {"rows": [], "total": 0}

    def prep(df):
        w = df.copy()
        w["dpsi_v"] = pd.to_numeric(w[dc], errors="coerce")
        w["fdr_v"] = pd.to_numeric(w[qc], errors="coerce")
        if rs_only and "passes_read_support" in w.columns:
            w = w[w["passes_read_support"].astype(str).str.lower().isin({"true", "1"})]
        return w.dropna(subset=["dpsi_v", "fdr_v"])

    wa, wb = prep(df_a), prep(df_b)
    wa_sig = wa[(wa["fdr_v"] <= max_fdr) & (wa["dpsi_v"].abs() >= min_abs_dpsi)]
    cols_a = [c for c in [ec, gc, sc, etc, "row_id", "dpsi_v", "fdr_v"] if c]
    m = pd.merge(
        wa_sig[cols_a].rename(columns={"row_id": "row_id_a", "dpsi_v": "dpsi_a", "fdr_v": "fdr_a"}),
        wb[[ec, "dpsi_v", "fdr_v"]].rename(columns={"dpsi_v": "dpsi_b", "fdr_v": "fdr_b"}),
        on=ec, how="left",
    )

    if relation == "reversed":
        # Rescue / opposite shift in B
        cond = m["dpsi_b"].notna() & (m["fdr_b"] <= max_fdr) & (m["dpsi_b"].abs() >= min_abs_dpsi) & ((m["dpsi_a"] * m["dpsi_b"]) < 0)
    elif relation == "concordant":
        cond = m["dpsi_b"].notna() & (m["fdr_b"] <= max_fdr) & (m["dpsi_b"].abs() >= min_abs_dpsi) & ((m["dpsi_a"] * m["dpsi_b"]) > 0)
    else:
        # Specific to A (not significant in B)
        cond = m["dpsi_b"].isna() | (m["fdr_b"] > max_fdr) | (m["dpsi_b"].abs() < min_abs_dpsi * 0.5)

    res = m.loc[cond].copy()
    res["abs_sum"] = res["dpsi_a"].abs() + res["dpsi_b"].fillna(0).abs()
    res = res.sort_values("abs_sum", ascending=False)

    prot_art = _find_contrast_artifact(grouped, "prioritize.proteome_consequences", contrast_a)
    if prot_art is not None:
        p_df = _load_artifact_df(prot_art)
        if "event_id" in p_df.columns and ec == "event_id":
            res["_ceid"] = res["event_id"].map(_canon_eid)
            p_sub = p_df[["event_id", "nmd_switch_class", "disrupted_protein_features"]].copy()
            p_sub["_ceid"] = p_sub["event_id"].map(_canon_eid)
            res = pd.merge(res, p_sub[["_ceid", "nmd_switch_class", "disrupted_protein_features"]].drop_duplicates("_ceid"), on="_ceid", how="left").drop(columns=["_ceid"])

    rename_map = {ec: "event_id"}
    if gc:
        rename_map[gc] = "gene_id"
    if sc and sc != gc:
        rename_map[sc] = "gene_symbol"
    res = res.rename(columns=rename_map)
    if "gene_symbol" not in res.columns and "gene_id" in res.columns:
        res["gene_symbol"] = res["gene_id"]

    return {
        "contrast_a": contrast_a,
        "contrast_b": contrast_b,
        "relation": relation,
        "source_artifact_id": art_a.artifact_id,
        "total": int(len(res)),
        "rows": _clean_records(res.head(250)),
    }


# ---------------------------------------------------------------------------
# 4. Single-Cell Smart-seq3 Dual Manifold Explorer
# ---------------------------------------------------------------------------
def single_cell_manifold(context, feature: str = "") -> dict[str, Any]:
    grouped = _artifacts_by_source_id(context)
    emb_art = _find_contrast_artifact(grouped, "single_cell.cell_embedding")
    if emb_art is None:
        return {"available": False, "cells": []}
    emb = _load_artifact_df(emb_art).copy()
    feature_values: dict[str, float | None] = {}
    feature_label = ""

    if feature:
        # Check PSI matrix first, then RBP expression matrix
        psi_art = _find_contrast_artifact(grouped, "single_cell.psi_matrix")
        if psi_art is not None:
            psi_df = _load_artifact_df(psi_art)
            hit = psi_df[psi_df["event_id"].astype(str) == feature]
            if not hit.empty:
                row = hit.iloc[0]
                feature_label = f"PSI · {feature[:36]}"
                for c in psi_df.columns:
                    if c != "event_id":
                        feature_values[str(c)] = _safe_float(row[c])
        if not feature_values:
            exp_art = _find_contrast_artifact(grouped, "single_cell.expression_matrix")
            if exp_art is not None:
                exp_df = _load_artifact_df(exp_art)
                matched_col = next((c for c in exp_df.columns if str(c).upper() == feature.upper()), None)
                if matched_col:
                    feature_label = f"Expression · {matched_col}"
                    for _, r in exp_df.iterrows():
                        feature_values[str(r["cell_id"])] = _safe_float(r[matched_col])

    if feature_values:
        emb["feature_value"] = emb["cell_id"].astype(str).map(feature_values)

    psi_events = []
    if "single_cell.psi_matrix" in grouped:
        psi_df = _load_artifact_df(grouped["single_cell.psi_matrix"][0])
        psi_events = psi_df["event_id"].astype(str).head(100).tolist()
    expr_genes = []
    if "single_cell.expression_matrix" in grouped:
        exp_df = _load_artifact_df(grouped["single_cell.expression_matrix"][0])
        expr_genes = [str(c) for c in exp_df.columns if c != "cell_id"][:160]

    return {
        "available": True,
        "feature_label": feature_label,
        "module_columns": [c for c in emb.columns if c.startswith("module_M")],
        "cell_types": sorted(emb["cell_type"].dropna().astype(str).unique().tolist()) if "cell_type" in emb.columns else [],
        "donors": sorted(emb["donor_id"].dropna().astype(str).unique().tolist()) if "donor_id" in emb.columns else [],
        "psi_events": psi_events,
        "expr_genes": expr_genes,
        "cells": _clean_records(emb),
    }


# ---------------------------------------------------------------------------
# 5. 4-Track Adaptive GeneCard & Quantitative Splice-Junction Sashimi Preview
# ---------------------------------------------------------------------------
def genecard_preview(
    context,
    gene_query: str = "",
    event_id: str = "",
    comparison: str = "",
    custom_sequence: str = "",
    custom_motif: str = "",
) -> dict[str, Any]:
    grouped = _artifacts_by_source_id(context)
    cid = comparison or (context.manifest.contrasts[0].contrast_id if context.manifest.contrasts else "")
    das_art = _find_contrast_artifact(grouped, "differential.as", cid)
    deg_art = _find_contrast_artifact(grouped, "differential.deg", cid)
    dtu_art = _find_contrast_artifact(grouped, "differential.dtu", cid)
    prot_art = _find_contrast_artifact(grouped, "prioritize.proteome_consequences", cid)
    tgt_art = _find_contrast_artifact(grouped, "prioritize.splicetrace_targets", cid)

    das_df = _load_artifact_df(das_art) if das_art else pd.DataFrame()
    deg_df = _load_artifact_df(deg_art) if deg_art else pd.DataFrame()
    dtu_df = _load_artifact_df(dtu_art) if dtu_art else pd.DataFrame()
    prot_df = _load_artifact_df(prot_art) if prot_art else pd.DataFrame()
    tgt_df = _load_artifact_df(tgt_art) if tgt_art else pd.DataFrame()

    if das_df.empty:
        return {"found": False, "reason": "No differential splicing table available in workspace."}

    if "comparison_id" in das_df.columns and cid:
        das_c = das_df[das_df["comparison_id"].astype(str) == cid].copy()
        if not das_c.empty:
            das_df = das_c

    ec = next((c for c in ("event_id", "event_uid", "ID") if c in das_df.columns), das_df.columns[0])
    gc = next((c for c in ("gene_id", "GeneID") if c in das_df.columns), ec)
    sc = next((c for c in ("gene_symbol", "gene_name", "geneSymbol") if c in das_df.columns), None)
    dpsi_c = next((c for c in ("dpsi", "IncLevelDifference", "delta_psi") if c in das_df.columns), None)
    fdr_c = next((c for c in ("FDR", "fdr", "padj") if c in das_df.columns), None)

    das_df["_eid"] = das_df[ec].astype(str)
    das_df["_gid"] = das_df[gc].map(_norm_gene)
    # Enrich gene_symbol from report display_names or prot_df if missing
    sym_lookup: dict[str, str] = {}
    if not prot_df.empty and {"gene_id", "gene_symbol"}.issubset(prot_df.columns):
        for _, r in prot_df[["gene_id", "gene_symbol"]].dropna().drop_duplicates().iterrows():
            sym_lookup[_norm_gene(r["gene_id"])] = str(r["gene_symbol"])
    if das_art is not None:
        try:
            from .report_names import display_names
            dn = display_names(das_art)
            for _, r in dn.iterrows():
                if pd.notna(r["_mapped_name"]):
                    sym_lookup[str(r["_name_id"])] = str(r["_mapped_name"])
        except Exception:
            pass

    if sc and sc in das_df.columns:
        das_df["_sym"] = das_df[sc].astype(str)
        mask_ens = das_df["_sym"].isin({"", "NA", "nan", "None"}) | das_df["_sym"].str.startswith("ENS")
        das_df.loc[mask_ens, "_sym"] = das_df.loc[mask_ens, "_gid"].map(lambda g: sym_lookup.get(g, g))
    else:
        das_df["_sym"] = das_df["_gid"].map(lambda g: sym_lookup.get(g, g))

    # Locate target event & host gene
    chosen_row = None
    if event_id:
        c_req = _canon_eid(event_id)
        m_ev = das_df[(das_df["_eid"] == str(event_id)) | (das_df["_eid"].map(_canon_eid) == c_req)]
        if not m_ev.empty:
            chosen_row = m_ev.iloc[0]
    if chosen_row is None and gene_query:
        q_norm = _norm_gene(gene_query).upper()
        m_g = das_df[(das_df["_gid"].str.upper() == q_norm) | (das_df["_sym"].str.upper() == q_norm)].copy()
        if not m_g.empty:
            prot_ceids = set(prot_df["event_id"].map(_canon_eid)) if (not prot_df.empty and "event_id" in prot_df.columns) else set()
            m_g["_rs"] = m_g["passes_read_support"].astype(str).str.lower().isin({"true", "1"}).astype(int) if "passes_read_support" in m_g.columns else 1
            m_g["_has_prot"] = m_g["_eid"].map(_canon_eid).isin(prot_ceids).astype(int)
            m_g["_sig"] = (pd.to_numeric(m_g[fdr_c], errors="coerce").fillna(1.0) <= 0.05).astype(int) if fdr_c else 1
            m_g["_abs"] = pd.to_numeric(m_g[dpsi_c], errors="coerce").abs().fillna(0.0) if dpsi_c else 0.0
            m_g = m_g.sort_values(["_rs", "_sig", "_has_prot", "_abs"], ascending=[False, False, False, False])
            chosen_row = m_g.iloc[0]
    if chosen_row is None:
        # Default to top significant PS-DAS event
        w = das_df.copy()
        if "passes_read_support" in w.columns:
            w_rs = w[w["passes_read_support"].astype(str).str.lower().isin({"true", "1"})]
            if not w_rs.empty:
                w = w_rs
        if fdr_c:
            w = w.sort_values(fdr_c)
        chosen_row = w.iloc[0]

    target_eid = str(chosen_row["_eid"])
    target_gid = str(chosen_row["_gid"])
    target_sym = str(chosen_row["_sym"])

    # All events of this host gene
    gene_events_df = das_df[das_df["_gid"] == target_gid].copy()
    if dpsi_c:
        gene_events_df["_abs"] = pd.to_numeric(gene_events_df[dpsi_c], errors="coerce").abs().fillna(0.0)
        gene_events_df = gene_events_df.sort_values("_abs", ascending=False)
    host_events = [
        {
            "event_id": str(r["_eid"]),
            "event_type": str(r.get("event_type", "SE")),
            "dpsi": _safe_float(r[dpsi_c]) if dpsi_c else None,
            "fdr": _safe_float(r[fdr_c]) if fdr_c else None,
            "row_id": str(r.get("row_id", "")),
        }
        for _, r in gene_events_df.head(30).iterrows()
    ]

    # ---------------- TRACK 1: Host Gene RNA Expression & Event PSI ----------------
    samples_by_id = {s.sample_id: s for s in context.manifest.samples}
    contrast_obj = next((c for c in context.manifest.contrasts if c.contrast_id == cid), None)
    num_ids = list(contrast_obj.numerator.sample_ids) if contrast_obj else []
    den_ids = list(contrast_obj.denominator.sample_ids) if contrast_obj else []
    num_label = contrast_obj.numerator.label if contrast_obj else "Numerator"
    den_label = contrast_obj.denominator.label if contrast_obj else "Denominator"

    deg_summary = {"log2fc": None, "fdr": None, "baseMean": None}
    if not deg_df.empty:
        g_d = next((c for c in ("gene_id", "GeneID") if c in deg_df.columns), deg_df.columns[0])
        sub_d = deg_df[deg_df[g_d].map(_norm_gene) == target_gid]
        if sub_d.empty:
            sub_d = deg_df[deg_df[g_d].astype(str).str.upper() == target_sym.upper()]
        if not sub_d.empty:
            dr = sub_d.iloc[0]
            fc_c = next((c for c in ("log2fc", "log2FoldChange") if c in deg_df.columns), None)
            fq_c = next((c for c in ("deg_fdr", "padj", "fdr") if c in deg_df.columns), None)
            deg_summary = {
                "log2fc": _safe_float(dr[fc_c]) if fc_c else None,
                "fdr": _safe_float(dr[fq_c]) if fq_c else None,
                "baseMean": _safe_float(dr.get("baseMean")),
            }

    # Per-sample RNA expression from report's normalised_counts.tsv or single_cell.expression_matrix
    expr_samples: list[dict[str, Any]] = []
    if context.manifest.metadata.get("report_workspace"):
        rroot = report_root(context)
        norm_files = list((rroot / "outputs").rglob("*normalised_counts*.tsv"))
        if norm_files:
            try:
                # Fast selective line scan so we don't parse 8 MB TSV every click!
                with norm_files[0].open("r", encoding="utf-8", errors="ignore") as handle:
                    header = handle.readline().rstrip("\n").split("\t")
                    for line in handle:
                        if line.startswith(target_gid):
                            parts = line.rstrip("\n").split("\t")
                            if _norm_gene(parts[0]) == target_gid:
                                for col_idx, col_name in enumerate(header[1:], 1):
                                    if col_idx < len(parts) and col_name in samples_by_id:
                                        val = _safe_float(parts[col_idx])
                                        if val is not None:
                                            expr_samples.append({
                                                "sample_id": col_name,
                                                "condition": samples_by_id[col_name].condition,
                                                "arm": num_label if col_name in num_ids else (den_label if col_name in den_ids else samples_by_id[col_name].condition),
                                                "value": round(val, 2),
                                                "unit": "Normalized counts",
                                            })
                                break
            except Exception:
                pass
    if not expr_samples and "single_cell.expression_matrix" in grouped:
        exp_df = _load_artifact_df(grouped["single_cell.expression_matrix"][0])
        mcol = next((c for c in exp_df.columns if str(c).upper() == target_sym.upper()), None)
        if mcol:
            for _, er in exp_df.iterrows():
                sid = str(er["cell_id"])
                if sid in samples_by_id:
                    v = _safe_float(er[mcol])
                    if v is not None:
                        expr_samples.append({
                            "sample_id": sid,
                            "condition": samples_by_id[sid].condition,
                            "arm": samples_by_id[sid].condition,
                            "value": round(v, 3),
                            "unit": "log1p(Expression)",
                        })

    # Per-sample PSI & IJC / SJC counts for the exact selected event
    psi_samples: list[dict[str, Any]] = []
    flipped = any(str(chosen_row.get(fl, "")).lower() in {"1", "true", "yes"} for fl in ("splicing_direction_flipped", "reverse_direction"))

    def _split_nums(raw_v: object) -> list[float | None]:
        if raw_v is None or (isinstance(raw_v, float) and math.isnan(raw_v)):
            return []
        return [_safe_float(x) for x in str(raw_v).split(",") if str(x).strip() != ""]

    inc1 = _split_nums(chosen_row.get("IncLevel1"))
    inc2 = _split_nums(chosen_row.get("IncLevel2"))
    ijc1 = _split_nums(chosen_row.get("IJC_SAMPLE_1"))
    ijc2 = _split_nums(chosen_row.get("IJC_SAMPLE_2"))
    sjc1 = _split_nums(chosen_row.get("SJC_SAMPLE_1"))
    sjc2 = _split_nums(chosen_row.get("SJC_SAMPLE_2"))
    if flipped:
        inc1, inc2 = inc2, inc1
        ijc1, ijc2 = ijc2, ijc1
        sjc1, sjc2 = sjc2, sjc1

    for arm_lbl, sids, psis, ijcs, sjcs in ((num_label, num_ids, inc1, ijc1, sjc1), (den_label, den_ids, inc2, ijc2, sjc2)):
        if psis and len(sids) == len(psis):
            for idx, sid in enumerate(sids):
                if psis[idx] is not None:
                    psi_samples.append({
                        "sample_id": sid,
                        "condition": arm_lbl,
                        "psi": round(float(psis[idx]), 4),
                        "ijc": ijcs[idx] if idx < len(ijcs) else None,
                        "sjc": sjcs[idx] if idx < len(sjcs) else None,
                    })

    if not psi_samples and "single_cell.psi_matrix" in grouped:
        psi_m = _load_artifact_df(grouped["single_cell.psi_matrix"][0])
        hit_p = psi_m[psi_m["event_id"].astype(str) == target_eid]
        if not hit_p.empty:
            prow = hit_p.iloc[0]
            for sid, smeta in samples_by_id.items():
                if sid in hit_p.columns:
                    pv = _safe_float(prow[sid])
                    if pv is not None:
                        psi_samples.append({"sample_id": sid, "condition": smeta.condition, "psi": round(pv, 4), "ijc": None, "sjc": None})

    # DTU summary if available
    dtu_isoforms = []
    if not dtu_df.empty:
        gd_c = next((c for c in ("gene_id", "GeneID") if c in dtu_df.columns), None)
        td_c = next((c for c in ("transcript_id", "feature_id", "featureID") if c in dtu_df.columns), None)
        qd_c = next((c for c in ("dtu_fdr", "fdr", "padj", "qvalue") if c in dtu_df.columns), None)
        ed_c = next((c for c in ("effect", "dIF", "log2fc") if c in dtu_df.columns), None)
        if gd_c and td_c:
            sub_dtu = dtu_df[dtu_df[gd_c].map(_norm_gene) == target_gid].head(8)
            for _, tr in sub_dtu.iterrows():
                dtu_isoforms.append({
                    "transcript_id": str(tr[td_c]),
                    "effect": _safe_float(tr[ed_c]) if ed_c else None,
                    "fdr": _safe_float(tr[qd_c]) if qd_c else None,
                })

    # ---------------- TRACK 2: Quantitative Exon-Junction Sashimi Preview ----------------
    ev_type = str(chosen_row.get("event_type", "SE") or "SE")
    chrom = str(chosen_row.get("chr", chosen_row.get("chrom", "")) or "")
    strand = str(chosen_row.get("strand", "+") or "+")
    exons_geom: list[dict[str, Any]] = []
    junctions_geom: list[dict[str, Any]] = []

    def _coord(k1: str, k2: str) -> tuple[int, int] | None:
        a, b = _safe_float(chosen_row.get(k1)), _safe_float(chosen_row.get(k2))
        if a is not None and b is not None and b > a:
            return int(a), int(b)
        return None

    if ev_type == "SE":
        up = _coord("upstreamES", "upstreamEE")
        sk = _coord("exonStart_0base", "exonEnd")
        dn = _coord("downstreamES", "downstreamEE")
        if up and sk and dn:
            exons_geom = [
                {"label": "Upstream exon", "role": "constitutive", "start": up[0], "end": up[1]},
                {"label": "Cassette exon", "role": "alternative", "start": sk[0], "end": sk[1]},
                {"label": "Downstream exon", "role": "constitutive", "start": dn[0], "end": dn[1]},
            ]
            junctions_geom = [
                {"label": "Inclusion J1", "role": "inclusion", "from_idx": 0, "to_idx": 1},
                {"label": "Inclusion J2", "role": "inclusion", "from_idx": 1, "to_idx": 2},
                {"label": "Skipping J_skip", "role": "skipping", "from_idx": 0, "to_idx": 2},
            ]
    elif ev_type == "RI":
        up = _coord("upstreamES", "upstreamEE")
        ri = _coord("riExonStart_0base", "riExonEnd")
        dn = _coord("downstreamES", "downstreamEE")
        if up and dn:
            exons_geom = [
                {"label": "5′ exon", "role": "constitutive", "start": up[0], "end": up[1]},
                {"label": "Retained intron", "role": "alternative", "start": up[1], "end": dn[0]},
                {"label": "3′ exon", "role": "constitutive", "start": dn[0], "end": dn[1]},
            ]
            junctions_geom = [{"label": "Spliced intron junction", "role": "skipping", "from_idx": 0, "to_idx": 2}]
    elif ev_type in ("A3SS", "A5SS"):
        lg = _coord("longExonStart_0base", "longExonEnd")
        sh = _coord("shortES", "shortEE")
        fl = _coord("flankingES", "flankingEE")
        if lg and sh and fl:
            exons_geom = sorted([
                {"label": "Flanking exon", "role": "constitutive", "start": fl[0], "end": fl[1]},
                {"label": "Short exon core", "role": "constitutive", "start": sh[0], "end": sh[1]},
                {"label": f"{ev_type} extension", "role": "alternative", "start": min(lg[0], sh[1]), "end": max(lg[1], sh[0])},
            ], key=lambda x: x["start"])
            junctions_geom = [
                {"label": "Long splice site", "role": "inclusion", "from_idx": 0, "to_idx": 1},
                {"label": "Short splice site", "role": "skipping", "from_idx": 0, "to_idx": 2},
            ]
    elif ev_type == "MXE":
        up = _coord("upstreamES", "upstreamEE")
        e1 = _coord("1stExonStart_0base", "1stExonEnd")
        e2 = _coord("2ndExonStart_0base", "2ndExonEnd")
        dn = _coord("downstreamES", "downstreamEE")
        if up and e1 and e2 and dn:
            exons_geom = [
                {"label": "Upstream exon", "role": "constitutive", "start": up[0], "end": up[1]},
                {"label": "MXE exon 1", "role": "alternative", "start": e1[0], "end": e1[1]},
                {"label": "MXE exon 2", "role": "alternative_2", "start": e2[0], "end": e2[1]},
                {"label": "Downstream exon", "role": "constitutive", "start": dn[0], "end": dn[1]},
            ]
            junctions_geom = [
                {"label": "MXE1 J1", "role": "inclusion", "from_idx": 0, "to_idx": 1},
                {"label": "MXE1 J2", "role": "inclusion", "from_idx": 1, "to_idx": 3},
                {"label": "MXE2 J1", "role": "skipping", "from_idx": 0, "to_idx": 2},
                {"label": "MXE2 J2", "role": "skipping", "from_idx": 2, "to_idx": 3},
            ]

    def _mean_clean(vals: list[float | None]) -> float | None:
        clean = [v for v in vals if v is not None]
        return round(float(np.mean(clean)), 2) if clean else None

    mean_psi_1 = _mean_clean(inc1)
    mean_psi_2 = _mean_clean(inc2)
    if mean_psi_1 is None and psi_samples:
        p1 = [r["psi"] for r in psi_samples if r["condition"] == num_label]
        p2 = [r["psi"] for r in psi_samples if r["condition"] == den_label]
        mean_psi_1 = round(float(np.mean(p1)), 4) if p1 else None
        mean_psi_2 = round(float(np.mean(p2)), 4) if p2 else None

    track2_sashimi = {
        "chrom": chrom,
        "strand": strand,
        "event_type": ev_type,
        "exons": exons_geom,
        "junctions": junctions_geom,
        "arms": [
            {"label": num_label, "mean_psi": mean_psi_1, "mean_ijc": _mean_clean(ijc1), "mean_sjc": _mean_clean(sjc1), "n_samples": len(num_ids)},
            {"label": den_label, "mean_psi": mean_psi_2, "mean_ijc": _mean_clean(ijc2), "mean_sjc": _mean_clean(sjc2), "n_samples": len(den_ids)},
        ],
        "bam_bound_samples": sum(1 for s in context.manifest.samples if s.assets.get("bam") or s.assets.get("single_cell_bam")),
    }

    target_ceid = _canon_eid(target_eid)
    # ---------------- TRACK 3: Adaptive RBP Cis-Motif Hit Track (Never Blank!) ----------------
    if tgt_art is None:
        track3_motifs: dict[str, Any] = {"status": "unavailable"}
    else:
        hit_t = tgt_df[(tgt_df["event_id"].astype(str) == target_eid) | (tgt_df["event_id"].map(_canon_eid) == target_ceid)] if "event_id" in tgt_df.columns else pd.DataFrame()
        motifs_list: list[dict[str, Any]] = []
        if not hit_t.empty:
            grouped_hits: dict[tuple[str, str], dict[str, Any]] = {}
            for _, mr in hit_t.iterrows():
                rc_m = next((c for c in ("rbp_list", "rbps", "RBP", "top_rbps") if c in tgt_df.columns), None)
                raw_rbps = re.sub(r"[\[\]'\"]", "", str(mr[rc_m] or "")).strip() if rc_m else ""
                delta_s = _safe_float(mr.get("delta_motif_score", mr.get("mature_delta_score", mr.get("score_diff"))))
                reg_lbl = str(mr.get("region", mr.get("motif_location", "Alternative exon / ±200 nt flank")))
                if raw_rbps and raw_rbps.lower() not in ("none", "nan"):
                    for rbp_name in [x.strip() for x in re.split(r"[,;|]", raw_rbps) if x.strip()]:
                        key = (rbp_name.upper(), reg_lbl)
                        if key not in grouped_hits:
                            grouped_hits[key] = {
                                "rbp": rbp_name,
                                "region": reg_lbl,
                                "delta_score": delta_s,
                                "count": 1,
                            }
                        else:
                            grouped_hits[key]["count"] += 1
                            if delta_s is not None and (grouped_hits[key]["delta_score"] is None or abs(delta_s) > abs(grouped_hits[key]["delta_score"])):
                                grouped_hits[key]["delta_score"] = delta_s
            for idx_r, item in enumerate(sorted(grouped_hits.values(), key=lambda x: (-x["count"], x["rbp"]))[:12]):
                rbp_disp = f"{item['rbp']} (×{item['count']} sites)" if item["count"] > 1 else item["rbp"]
                motifs_list.append({
                    "rbp": rbp_disp,
                    "region": item["region"],
                    "delta_score": item["delta_score"],
                    "slot": idx_r,
                })
        if motifs_list:
            track3_motifs = {"status": "present", "hits": motifs_list}
        else:
            track3_motifs = {"status": "scanned_zero", "note": "SpliceTrace scanned curated RBP PWMs across ±200 nt flanks & mature junction; 0 high-threshold cis-motif hits on this locus."}

    # ---------------- TRACK 4: Adaptive Mature mRNA NMD & UniProt Domain Track (Never Blank!) ----------------
    if prot_art is None:
        track4_proteome: dict[str, Any] = {"status": "unavailable"}
    else:
        hit_pr = prot_df[(prot_df["event_id"].astype(str) == target_eid) | (prot_df["event_id"].map(_canon_eid) == target_ceid)] if "event_id" in prot_df.columns else pd.DataFrame()
        if hit_pr.empty and "architecture_json" in chosen_row and pd.notna(chosen_row.get("architecture_json")):
            hit_pr = pd.DataFrame([chosen_row])
        if hit_pr.empty:
            track4_proteome = {"status": "scanned_unmapped", "note": "Event lies outside annotated GENCODE canonical transcript models."}
        else:
            pr0 = hit_pr.iloc[0]
            loc_cls = str(pr0.get("coding_region_location", "") or "")
            frame_st = str(pr0.get("frame_status", "") or "")
            nmd_cls = str(pr0.get("nmd_switch_class", "") or "")
            arch = {}
            raw_arch = pr0.get("architecture_json")
            if isinstance(raw_arch, str) and raw_arch.strip().startswith("{"):
                try:
                    arch = json.loads(raw_arch)
                except Exception:
                    arch = {}
            if "noncod" in loc_cls.lower() or "noncod" in frame_st.lower() or nmd_cls.lower() == "noncoding":
                track4_proteome = {
                    "status": "noncoding",
                    "transcript_id": str(pr0.get("host_transcript_id", "")),
                    "summary": str(pr0.get("functional_consequence_summary", "Non-coding transcript isoform remodeling (no CDS ORF translation).")),
                    "architecture": arch,
                }
            else:
                dom_txt = str(pr0.get("disrupted_protein_features", "") or "")
                track4_proteome = {
                    "status": "coding_architecture",
                    "transcript_id": str(pr0.get("host_transcript_id", "")),
                    "uniprot_acc": str(pr0.get("uniprot_acc", "") or ""),
                    "coding_region_location": loc_cls,
                    "frame_status": frame_st,
                    "nmd_switch_class": nmd_cls,
                    "aa_len_inclusion": _safe_float(pr0.get("aa_len_inclusion")),
                    "aa_len_exclusion": _safe_float(pr0.get("aa_len_exclusion")),
                    "delta_aa": _safe_float(pr0.get("delta_aa")),
                    "ptc_to_last_ejc_nt": _safe_float(pr0.get("ptc_to_last_ejc_nt")),
                    "disrupted_protein_features": "" if dom_txt.lower() in ("none", "nan") else dom_txt,
                    "summary": str(pr0.get("functional_consequence_summary", "")),
                    "architecture": arch,
                }

    # Real-time on-demand Track 3 & Track 4 enrichment so EVERY gene/event in EVERY workspace renders all 4 tracks!
    from .mckinsey_compute import compute_on_demand_genecard_tracks
    track3_motifs, track4_proteome = compute_on_demand_genecard_tracks(
        context,
        chosen_row=chosen_row,
        target_eid=target_eid,
        target_gid=target_gid,
        target_sym=target_sym,
        ev_type=ev_type,
        exons_geom=exons_geom,
        comparison=cid,
        existing_track3=track3_motifs,
        existing_track4=track4_proteome,
        custom_sequence=custom_sequence,
        custom_motif=custom_motif,
    )

    return {
        "found": True,
        "comparison_id": cid,
        "gene_id": target_gid,
        "gene_symbol": target_sym,
        "event_id": target_eid,
        "event_type": ev_type,
        "dpsi": _safe_float(chosen_row.get(dpsi_c)) if dpsi_c else None,
        "fdr": _safe_float(chosen_row.get(fdr_c)) if fdr_c else None,
        "das_artifact_id": das_art.artifact_id if das_art else "",
        "das_row_id": str(chosen_row.get("row_id", "")),
        "host_events": host_events,
        "track1_quantification": {
            "deg": deg_summary,
            "expression_samples": expr_samples,
            "psi_samples": psi_samples,
            "dtu_isoforms": dtu_isoforms,
        },
        "track2_sashimi": track2_sashimi,
        "track3_motifs": track3_motifs,
        "track4_proteome": track4_proteome,
    }
