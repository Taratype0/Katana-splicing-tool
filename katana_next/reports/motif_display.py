"""Select report motifs from recorded native statistics without fitting tests.

The FDR is attached to one recorded window/full-sequence test in the complete
comparison family. It is not a newly computed, whole-model FDR. Nonsignificant
models and their original statistics remain in the complete output files.
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd

from .sources import digest


NATIVE_TEST_FAMILY = "Native complete comparison window/full-sequence BH family; not model-level FDR"
FDR_COLUMN = "q_new_window_family"
P_COLUMN = "raw_P_new_paired_swaps"


def significant_motif_ranking(ranking, *, threshold=.05, limit=10):
    """Return up to ``limit`` distinct eligible models, keeping native row values."""
    if not 0 < threshold < 1 or not isinstance(limit, int) or limit < 1:
        raise ValueError("Motif display requires a valid FDR threshold and positive limit")
    required = {"model_id", FDR_COLUMN, P_COLUMN}
    if not required <= set(ranking):
        raise ValueError("Motif display requires the complete recorded model ranking and paired P/FDR")
    if ranking.model_id.isna().any() or ranking.model_id.astype(str).duplicated().any():
        raise ValueError("Motif display requires one recorded test per distinct model")
    for column in (FDR_COLUMN, P_COLUMN):
        values = pd.to_numeric(ranking.loc[ranking[column].notna(), column], errors="raise")
        if not np.isfinite(values).all() or not values.between(0., 1.).all():
            raise ValueError("Recorded native motif P/FDR must be finite in [0, 1] or missing")
    eligible = ranking.loc[ranking[FDR_COLUMN].le(threshold) & ranking[P_COLUMN].notna()]
    return eligible.sort_values([FDR_COLUMN, P_COLUMN, "model_id"], kind="stable").head(limit).copy()


def motif_display_summary(ranking, comparison_id, *, threshold=.05, limit=10):
    """A compact scoped status, including a truthful zero-significance outcome."""
    selected = significant_motif_ranking(ranking, threshold=threshold, limit=limit)
    valid = ranking[FDR_COLUMN].notna() & ranking[P_COLUMN].notna()
    count = int((ranking[FDR_COLUMN].le(threshold) & ranking[P_COLUMN].notna()).sum())
    status = "significant" if count else "no_significant_motif" if valid.any() else "not_estimable"
    if status == "no_significant_motif":
        notice = (f"No motif has a recorded native test with FDR ≤ {threshold:g}. "
                  "No nonsignificant sequence panels are shown; all models and original statistics remain in the outputs.")
    elif status == "not_estimable":
        notice = ("No native motif test is estimable for this comparison. "
                  "All recorded models and source outputs are retained; no significance is inferred.")
    else:
        notice = (f"{count} models have at least one recorded test with FDR ≤ {threshold:g}; "
                  f"showing {len(selected)} distinct models ordered by that test's FDR, paired P and model ID.")
    return {"schema": "katana.motif_significance_display/v1", "chapter": 7,
            "contrast_id": str(comparison_id), "status": status, "threshold": float(threshold),
            "models": len(ranking), "estimable_models": int(valid.sum()),
            "significant_models": count, "displayed_models": len(selected),
            "min_native_fdr": float(ranking.loc[valid, FDR_COLUMN].min()) if valid.any() else None,
            "main_model_ids": selected.model_id.astype(str).tolist(),
            "native_test_family": NATIVE_TEST_FAMILY, "model_level_fdr_created": False,
            "notice": notice}


def apply_motif_significance_display(manifest, report_root, *, threshold=.05, limit=10):
    """Return an independent display-only manifest for an existing portable report.

    Native ranking files and all source panels are retained. This helper only
    changes which existing sequence cards are in the body and adds scoped
    summaries. Missing significant cards are rejected rather than hidden.
    """
    root = Path(report_root).resolve()
    result = deepcopy(manifest)
    summaries, source_hashes = [], {}
    candidates = [r for r in result.get("complete_results", [])
                  if r.get("name", "").startswith("all-motif-ranked-recorded-tests | ")]
    seen = set()
    for table in candidates:
        cid = str(table.get("contrast_id", ""))
        if not cid or cid in seen:
            raise ValueError("Motif report requires one complete ranking per comparison")
        seen.add(cid)
        path = (root / table["source"]).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError("Complete motif ranking must be an existing bounded report member")
        source_hashes[table["source"]] = digest(path)
        ranking = pd.read_csv(path, sep="\t", dtype={"model_id": str}, float_precision="round_trip")
        if "comparison_id" not in ranking or not ranking.comparison_id.astype(str).eq(cid).all():
            raise ValueError("Complete motif ranking belongs to another comparison")
        summary = motif_display_summary(ranking, cid, threshold=threshold, limit=limit)
        summary.update(ranking_path=table["source"], ranking_sha256=source_hashes[table["source"]])
        selected = set(summary["main_model_ids"])
        panels = [p for p in result.get("panels", []) if p.get("contrast_id") == cid
                  and p.get("metadata", {}).get("motif_view_kind") == "top_motif"]
        actual = [str(p["metadata"].get("model_id")) for p in panels]
        if len(actual) != len(set(actual)) or not selected <= set(actual):
            raise ValueError("Existing report omits or duplicates an eligible significant motif card")
        for panel in panels:
            token = str(panel["metadata"]["model_id"])
            eligible = token in selected
            panel["output_only"] = not eligible
            panel["metadata"].update(main_display_rule="native selected-test FDR <= threshold; max 10 distinct models",
                                     main_display_fdr_threshold=float(threshold),
                                     main_display_eligible=eligible)
        if summary["status"] != "significant":
            for panel in result.get("panels", []):
                if (panel.get("contrast_id") == cid and
                        panel.get("metadata", {}).get("component") == "motiftrace_sequence"):
                    panel["output_only"] = True
        summaries.append(summary)
    motif_cids = {str(p["contrast_id"]) for p in result.get("panels", [])
                  if p.get("metadata", {}).get("motif_view_kind") == "top_motif"}
    if not motif_cids <= seen:
        raise ValueError("Sequence cards without their complete recorded ranking cannot be significance-filtered")
    for source, sha in source_hashes.items():
        if digest(root / source) != sha:
            raise RuntimeError("Native motif ranking changed during display selection")
    result["motif_summaries"] = summaries
    return result, {"schema": "katana.motif_significance_display_acceptance/v1", "passed": True,
                    "comparisons": len(summaries), "comparison_summaries": summaries,
                    "main_sequence_cards": sum(s["displayed_models"] for s in summaries),
                    "native_ranking_sha256": source_hashes, "native_statistics_unchanged": True,
                    "new_scans": 0, "new_statistical_tests": 0}
