"""Positive, comparison-complete MotifTrace calculation/report acceptance.

The caller authenticates the collection and task contracts first. This gate
checks the full frozen scientific family, preserved numbers and actual Top 10
presentation. It never rescans sequences, refits tests or treats an absent
``unavailable`` message as proof that a comparison was rendered.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .model import file_sha256, read_json


WORKFLOW = "prioritize.motiftrace_consequences"


def _inside(root, relative):
    relative = Path(relative)
    path = (Path(root).resolve() / relative).resolve()
    if relative.is_absolute() or not path.is_relative_to(Path(root).resolve()) or not path.is_file():
        raise ValueError("MotifTrace report member is absent or escapes the report")
    return path


def _table(path):
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t", float_precision="round_trip", dtype={"model_id": str, "event_id": str})


def _empty_selection(task, receipt, receipts, transported, cid):
    """Recheck the actual full AS source; a claimed zero count is insufficient."""
    from katana_next.science.das_selection import selected_event_mask
    selection = receipt.get("selection") or {}
    expected = {"selected_rows": 0, "event_set": "stable_significant", "contrast_ids": [cid]}
    request = task.get("selection") or {}
    if (receipt.get("status") != "not_applicable" or
        request.get("analysis_id") != "differential.as" or request.get("contrast_ids") != [cid] or
        request.get("event_set") != "stable_significant" or
        any(selection.get(key) != value for key, value in expected.items()) or
        any(receipt.get("execution_evidence", {}).get(key) != selection.get(key)
            for key in ("source_artifact_id", "source_rows", *expected))):
        raise ValueError("MotifTrace N/A requires the exact zero RS foreground selection evidence")
    source_id = selection.get("source_artifact_id")
    sources = [a for a in receipt.get("artifacts", []) if a.get("artifact_id") == source_id
               and a.get("analysis_id") == "differential.as"]
    if len(sources) != 1:
        raise ValueError("MotifTrace zero foreground has no unique full differential-AS source")
    source = sources[0]
    owners = [r for r in receipts.values() if r.get("status") in {"completed", "reused"}
              and any(item.get("path") == source.get("path") for item in r.get("outputs", []))]
    if not owners:
        raise ValueError("MotifTrace zero foreground source is not a completed transported task output")
    frame = _table(transported(owners[0], source["path"]))
    if "comparison_id" not in frame or cid not in source.get("contrast_ids", source.get("parameters", {}).get("contrast_ids", [])):
        raise ValueError("MotifTrace zero foreground source does not establish its comparison")
    frame = frame.loc[frame.comparison_id.astype(str).eq(cid)]
    if len(frame) != selection.get("source_rows") or selected_event_mask(frame, "stable_significant").any():
        raise ValueError("MotifTrace N/A contradicts the actual measured source and RS foreground")
    audits = [a for a in receipt.get("generated_artifacts", []) if a.get("entity_type") == "eligibility_audit"]
    if len(audits) != 1:
        raise ValueError("MotifTrace zero foreground lacks its completed eligibility audit")
    audit = _table(transported(receipt, audits[0]["path"]))
    if len(audit) != 1:
        raise ValueError("MotifTrace zero foreground audit must contain its exact one selection")
    row = audit.iloc[0]
    for key in ("source_artifact_id", "source_rows", "selected_rows", "event_set"):
        if row.get(key) != selection.get(key):
            raise ValueError("MotifTrace eligibility audit differs from the actual selection")
    if list(row.get("contrast_ids", [])) != [cid]:
        raise ValueError("MotifTrace eligibility audit belongs to another comparison")


def _family(receipt, transported, cid):
    from katana_next.reports.motif_views import SEQUENCE_FILES, STAT_FILES
    if receipt.get("status") not in {"completed", "reused"} or receipt.get("job", {}).get("status") not in {"complete", "cache_hit"}:
        raise ValueError("MotifTrace requires an actually completed scientific workflow")
    markers = [r for r in receipt.get("outputs", []) if Path(r["path"]).name == "COMPLETE.json"]
    if len(markers) != 1:
        raise ValueError("MotifTrace requires exactly one complete unsplit family per comparison")
    marker = markers[0]
    complete_path = transported(receipt, marker["path"])
    complete = read_json(complete_path)
    members = complete.get("members", {})
    if not complete.get("complete") or not set((*SEQUENCE_FILES, *STAT_FILES)) <= set(members):
        raise ValueError("MotifTrace complete family omits required sequence/statistical outputs")
    original = Path(marker["path"]).parent
    paths = {}
    for name, record in members.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("MotifTrace complete family has an unbounded member")
        path = transported(receipt, str(original / relative))
        if path.stat().st_size != record.get("bytes") or file_sha256(path) != record.get("sha256"):
            raise ValueError("MotifTrace full scientific family changed or was truncated")
        paths[name] = path
    scientific = read_json(paths[SEQUENCE_FILES[0]])
    models = scientific.get("models", {})
    if (scientific.get("comparison_id") != cid or not isinstance(models, dict) or not models or
        complete.get("statistical_family_split", scientific.get("statistical_family_split", scientific.get("provenance", {}).get("statistical_family_split"))) is not False):
        raise ValueError("MotifTrace requires its complete declared model family and exact comparison")
    selected = _table(paths["FROZEN_SELECTED_EVENTS.tsv"])
    if selected.empty or selected.event_id.duplicated().any():
        raise ValueError("A completed MotifTrace family requires nonempty exact selected events")
    if "comparison_id" in selected and not selected.comparison_id.eq(cid).all():
        raise ValueError("Frozen MotifTrace events include another comparison")
    tests = _table(paths["NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz"])
    expected_columns = {"model_id", "event_type", "direction", "region", "method", "window_nt", "window_start0",
                        "raw_P_new_paired_swaps", "q_new_window_family", "mean_paired_score_difference"}
    if not expected_columns <= set(tests):
        if not tests.empty:
            raise ValueError("Native MotifTrace test family lacks explicit identity/statistical columns")
        from katana_next.reports.motif_views import verify_empty_completed_statistics
        verify_empty_completed_statistics(complete_path.parent, members)
        tests = tests.reindex(columns=sorted(set(tests) | expected_columns))
    if not set(tests.model_id.astype(str)) <= set(models):
        raise ValueError("Native MotifTrace tests contain undeclared models")
    if not tests.empty and set(tests.model_id.astype(str)) != set(models):
        raise ValueError("Native MotifTrace tests dropped part of the declared model family")
    for column in ("raw_P_new_paired_swaps", "q_new_window_family"):
        values = pd.to_numeric(tests.loc[tests[column].notna(), column], errors="raise")
        if not np.isfinite(values).all() or not values.between(0, 1).all():
            raise ValueError("Native MotifTrace P/FDR must be finite in [0, 1] or missing")
    return marker, complete_path, complete, paths, scientific, tests


def _report_family(report, report_dir, cid, marker, complete_path, complete, paths, scientific, tests):
    model_ids = set(scientific["models"])
    expected_count = min(10, len(model_ids))
    cards = [p for p in report.get("panels", []) if p.get("contrast_id") == cid and
             p.get("metadata", {}).get("component") == "motiftrace_sequence" and
             p.get("metadata", {}).get("motif_view_kind") == "top_motif" and not p.get("output_only")]
    if len(cards) != expected_count:
        raise ValueError("Every completed MotifTrace comparison requires its own min(10, all models) main motif cards")
    if len({p.get("id") for p in cards}) != expected_count:
        raise ValueError("MotifTrace main cards have duplicate panel identities")
    complete_sha = file_sha256(complete_path)
    ranking_paths = set()
    test_sha = file_sha256(paths["NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz"])
    for card in cards:
        meta = card["metadata"]
        if (meta.get("complete_manifest_sha256") != complete_sha or meta.get("statistical_family_split") is not False or
            meta.get("das_scope", {}).get("recorded_event_set") != "stable_significant" or card.get("policy") is not None):
            raise ValueError("MotifTrace main cards changed their completed family or fixed RS scope")
        image = _inside(report_dir, card.get("image", ""))
        if file_sha256(image) != meta.get("motif_image_sha256"):
            raise ValueError("MotifTrace main motif image is missing or changed")
        links = card.get("data_links", [])
        if any(link.get("contrast_id") != cid or link.get("policy") != "stable_significant" for link in links):
            raise ValueError("MotifTrace main card links mix comparison or measurement scopes")
        ranked = [link for link in links if link.get("label") == "All ranked models and native test statistics"]
        native = [link for link in links if link.get("role") == "full_results" and link.get("source_sha256") == test_sha]
        if len(ranked) != 1 or not native:
            raise ValueError("MotifTrace cards lack complete ranked models or the unchanged full native test family")
        ranking_paths.add(_inside(report_dir, ranked[0]["path"]))
        if any(file_sha256(_inside(report_dir, link["path"])) != link.get("portable_sha256") for link in links):
            raise ValueError("MotifTrace card data link changed after report construction")
        if any(file_sha256(_inside(report_dir, link["path"])) != test_sha for link in native):
            raise ValueError("MotifTrace native P/FDR rows were replaced in the report")
    if len(ranking_paths) != 1:
        raise ValueError("A comparison's MotifTrace cards must share one complete declared ranking")
    ranking = _table(next(iter(ranking_paths)))
    if (len(ranking) != len(model_ids) or ranking.model_id.duplicated().any() or set(ranking.model_id) != model_ids or
        not ranking.comparison_id.eq(cid).all() or not ranking.report_policy.eq("stable_significant").all()):
        raise ValueError("MotifTrace complete model ranking is missing, truncated or belongs to another comparison")
    keys = ["q_new_window_family", "raw_P_new_paired_swaps", "event_type", "direction", "region", "method", "window_nt", "window_start0"]
    valid = tests.loc[tests.q_new_window_family.notna() & tests.raw_P_new_paired_swaps.notna()]
    selected = {str(token): part.sort_values(keys, kind="stable").iloc[0] for token, part in valid.groupby("model_id", sort=False)}
    for row in ranking.to_dict("records"):
        native = selected.get(str(row["model_id"]))
        if native is None:
            if pd.notna(row["q_new_window_family"]) or pd.notna(row["raw_P_new_paired_swaps"]) or row["significance"] != "not_estimable":
                raise ValueError("MotifTrace report invented statistics for a non-estimable recorded model")
            continue
        for column in (*keys, "mean_paired_score_difference"):
            if not (pd.isna(row[column]) and pd.isna(native[column])) and row[column] != native[column]:
                raise ValueError("MotifTrace report P/FDR/effect must come from the same actual selected native test row")
        expected = "recorded_test_FDR_le_0.05" if row["q_new_window_family"] <= .05 else "not_significant"
        if row["significance"] != expected:
            raise ValueError("MotifTrace report changed a native test's significance label")
    ordered = ranking.sort_values(["q_new_window_family", "raw_P_new_paired_swaps", "model_id"], kind="stable", na_position="last")
    if ranking.model_id.tolist() != ordered.model_id.tolist() or ranking.display_rank.tolist() != list(range(1, len(ranking) + 1)):
        raise ValueError("MotifTrace ranking changed the declared native FDR/P/model order")
    expected_top = ranking.head(expected_count)
    actual_top = sorted((p["metadata"].get("display_rank"), p["metadata"].get("model_id")) for p in cards)
    if actual_top != list(zip(expected_top.display_rank, expected_top.model_id)):
        raise ValueError("MotifTrace main cards omit or substitute a ranked model")
    # Every frozen output remains a complete downloadable source, even when it
    # is not selected for the report body. Never certify only observed Top 10.
    members = {name: record["sha256"] for name, record in complete["members"].items()}
    members["COMPLETE.json"] = complete_sha
    for name, expected_sha in members.items():
        copies = [item for item in report.get("complete_results", []) if item.get("contrast_id") == cid and
                  item.get("name") == f"Complete native MotifTrace output | {name} | {cid}"]
        if len(copies) != 1 or file_sha256(_inside(report_dir, copies[0]["source"])) != expected_sha:
            raise ValueError("MotifTrace report download omits part of the complete scientific family: " + name)


def required_motif_report_gaps(root, plan, report_manifest, report_dir, *, project, receipts, transported):
    """Return typed gaps for all planned comparisons, including omitted tasks."""
    expected = plan.get("contrast_ids", [])
    declared = {c["contrast_id"] for c in project.get("contrasts", [])}
    gaps = []
    if not expected or len(expected) != len(set(expected)) or not set(expected) <= declared:
        return [{"component": "motiftrace_sequence", "reason": "Complete MotifTrace acceptance requires unique declared planned comparisons"}]
    tasks = [t for t in plan.get("tasks", []) if t.get("workflow_id") == WORKFLOW]
    for task in tasks:
        scope = task.get("scope", {}).get("contrast_ids", [])
        if len(scope) != 1 or scope[0] not in expected:
            gaps.append({"component": "motiftrace_sequence", "task_id": task["id"],
                         "reason": "Each MotifTrace task must cover exactly one selected declared comparison"})
    for cid in expected:
        matches = [t for t in tasks if t.get("scope", {}).get("contrast_ids") == [cid]]
        try:
            if len(matches) != 1:
                raise ValueError("Each selected comparison requires exactly one independent complete MotifTrace task")
            task = matches[0]
            receipt = receipts[task["id"]]
            if receipt.get("status") == "not_applicable":
                _empty_selection(task, receipt, receipts, transported, cid)
                if any(p.get("contrast_id") == cid and p.get("metadata", {}).get("motif_view_kind") == "top_motif"
                       for p in report_manifest.get("panels", [])):
                    raise ValueError("An empty RS foreground cannot have computed motif cards")
            else:
                family = _family(receipt, transported, cid)
                _report_family(report_manifest, report_dir, cid, *family)
        except (ValueError, OSError, KeyError, TypeError, AttributeError) as exc:
            gaps.append({"component": "motiftrace_sequence", "contrast_id": cid, "reason": str(exc)})
    return gaps
