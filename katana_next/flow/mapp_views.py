"""Read-only native MAEI views; never relabel its statistics as calibration."""
from __future__ import annotations

import json
from pathlib import Path, PureWindowsPath
import re


def build_mapp_views(root, plan, receipts, transport):
    import numpy as np
    import pandas as pd
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap
    from .model import file_sha256, read_json, write_json
    from .runtime import canonical_hash
    root = Path(root).resolve()
    project = read_json(root / "project.json")
    translated = {}
    for item in transport:
        relative = Path(item["collection_path"])
        path = (root / relative).resolve()
        original = item["source"]["path"]
        if (relative.is_absolute() or '\\' in item["collection_path"] or PureWindowsPath(str(relative)).drive or '..' in relative.parts
                or not path.is_relative_to(root)):
            raise ValueError("MAPP transported member escapes its declared collection")
        if original in translated:
            raise ValueError("MAPP transported member has a duplicate source identity")
        translated[original] = path
    families, not_applicable = {}, []
    for task in plan["tasks"]:
        if task.get("adapter") != "mapp.maei":
            continue
        _, receipt = receipts[task["id"]]
        if receipt.get("status") not in {"completed", "reused", "not_applicable"}:
            raise ValueError("MAPP views require an actual completed native fit")
        if receipt.get("receipt_sha256") != canonical_hash({k: v for k, v in receipt.items() if k != "receipt_sha256"}):
            raise ValueError("MAPP native completion proof changed")

        def verified(original):
            members = [row for row in receipt.get("outputs", []) if row["path"] == str(original)]
            if len(members) != 1 or str(original) not in translated:
                raise ValueError("MAPP result is not a verified transported member of its actual task")
            path = translated[str(original)]
            if file_sha256(path) != members[0]["sha256"]:
                raise ValueError("MAPP native numerical output changed")
            return path

        if receipt.get("status") == "not_applicable":
            from .acceptance import verify_mapp_not_applicable
            ids = task.get("scope", {}).get("contrast_ids", [])
            contrast = next((c for c in project.get("contrasts", []) if ids == [c["contrast_id"]]), None)
            if contrast is None:
                raise ValueError("MAPP non-applicability has no exact declared contrast")
            sample_ids = contrast["numerator"]["sample_ids"] + contrast["denominator"]["sample_ids"]
            proof = verify_mapp_not_applicable(task, receipt, verified, expected_sample_ids=sample_ids)
            not_applicable.append((task, proof))
            continue
        if receipt.get("status") not in {"completed", "reused"} or receipt.get("execution_evidence", {}).get("native_fitting_performed") is not True:
            raise ValueError("MAPP views require an actual completed native fit")
        members = [r for r in receipt["outputs"] if Path(r["path"]).name == "AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv"]
        if len(members) != 1 or members[0]["path"] not in translated:
            raise ValueError("MAPP native significance table is not a verified transported member")
        member = members[0]
        path = translated[member["path"]]
        if file_sha256(path) != member["sha256"]:
            raise ValueError("MAPP native numerical output changed")
        contrasts = task.get("scope", {}).get("contrast_ids", [])
        if len(contrasts) != 1:
            raise ValueError("MAPP report views require one exact comparison per fit")
        frame = pd.read_csv(path, sep="\t", index_col=0)
        frame.index = frame.index.astype(str)
        if frame.empty or "combined.Zscore" not in frame or not any(str(c).startswith("activities_") for c in frame):
            raise ValueError("Native MAEI table is missing its actual sample activity or Z statistics")
        ranking = pd.to_numeric(frame["combined.Zscore"], errors="coerce")
        # Preview selection only. Native full tables and correction families
        # remain separate; no global P/FDR is calculated across windows.
        preview = frame.loc[ranking.nlargest(20).index].copy()
        preview.insert(0, "native_model_id", preview.index)
        preview.insert(0, "fit_task_id", task["id"])
        preview.insert(0, "site", task.get("scope", {}).get("site", "declared native window"))
        families.setdefault(contrasts[0], []).append((task, path, preview))
    panels, tables, files, unavailable = [], [], [], []
    for task, proof in not_applicable:
        cid = proof["comparison_id"]
        tables.append({"chapter": 7, "name": f"Native MAEI eligibility audit | {cid} | {task['id']}",
                       "path": proof["verified_eligibility_table"], "analysis_id": "prioritize.mapp_native",
                       "contrast_id": cid})
        unavailable.append({"chapter": 7, "title": "Native MAEI not applicable | " + cid + " | " + task["id"],
                            "status": "not_applicable", "task_id": task["id"],
                            "reason": "The complete measured SE/transcript eligibility audit found no exact mapped included/skipped TPM targets for this native window. No native model or P value was fabricated.",
                            "metadata": {"component": "mapp_native", "reason_code": "no_exact_mapped_se",
                                         "measured_events": proof["measured_events"], "eligible_events": 0}})
    cmap = LinearSegmentedColormap.from_list("katana-maei", ["#2863A5", "#FFFFFF", "#C63F4A"])
    for cid, entries in families.items():
        directory = root / "presentation" / "mapp" / (re.sub(r"[^A-Za-z0-9_.-]", "_", cid)[:70] + "-" + canonical_hash(cid)[:12])
        directory.mkdir(parents=True, exist_ok=True)
        frame = pd.concat([item[2] for item in entries], ignore_index=True)
        preview = frame.nlargest(20, "combined.Zscore").copy()
        preview_path = directory / "native_display_rows.tsv"
        preview.to_csv(preview_path, sep="\t", index=False)
        links = [{"role": "plot_data", "label": "Exact plotted native rows", "path": str(preview_path),
                  "contrast_id": cid, "policy": None, "unit": "native motif/window/sample activity"}]
        catalog = [{"task_id": task["id"], "path": path.relative_to(root).as_posix(),
                    "sha256": file_sha256(path), "scope": task["scope"]} for task, path, _ in entries]
        index_path = write_json(directory / "complete_native_fits.json", {"comparison_id": cid, "fits": catalog,
            "correction": "Native MAEI per-window/sample Bonferroni; no pooled cross-window FDR",
            "calibrated_motiftrace_parallel": False})
        links.append({"role": "full_results", "label": "Every completed native fit and original statistics",
                      "path": str(index_path), "source_collection_root": str(root),
                      "contrast_id": cid, "policy": None, "unit": "native motif/window fit"})
        tables.extend({"chapter": 7, "name": f"Native MAPP MAEI | {cid} | {task['id']}", "path": str(path),
                       "analysis_id": "prioritize.mapp_native", "contrast_id": cid} for task, path, _ in entries)
        activity_cols = [c for c in preview if str(c).startswith("activities_")]
        values = preview[activity_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
        labels = [str(site) + " | " + str(model) for site, model in zip(preview.site, preview.native_model_id)]
        finite = np.abs(values[np.isfinite(values)])
        bound = float(finite.max()) if len(finite) else 1.0
        fig = Figure(figsize=(10, max(4, .28 * len(preview))))
        FigureCanvasAgg(fig)
        ax = fig.subplots()
        image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=-bound or -1, vmax=bound or 1, aspect="auto")
        ax.set_xticks(range(len(activity_cols)), [c.removeprefix("activities_") for c in activity_cols], rotation=45, ha="right")
        ax.set_yticks(range(len(preview)), labels, fontsize=7)
        ax.set_title("Native MAEI motif activity | " + cid)
        fig.colorbar(image, ax=ax, label="Original native sample activity")
        fig.tight_layout()
        path = directory / "native_activity_heatmap.png"
        fig.savefig(path, dpi=150); fig.clear()
        panels.append({"id": "native-maei-" + cid, "chapter": 7, "contrast_id": cid,
                       "title": "Native MAPP MAEI motif activities | " + cid, "path": str(path),
                       "data_links": links,
                       "metadata": {"analysis_id": "prioritize.mapp_native", "component": "mapp_native",
                                    "fit_task_ids": [task["id"] for task, _, _ in entries],
                                    "fit_receipt_sha256": {task["id"]: receipts[task["id"]][1]["receipt_sha256"]
                                                           for task, _, _ in entries},
                                    "native_fit_index_sha256": file_sha256(index_path),
                                    "native_plot_data_sha256": file_sha256(preview_path),
                                    "native_activity_image_sha256": file_sha256(path),
                                    "calibrated_motiftrace_parallel": False},
                       "caption": "Original native sample activities. Displayed rows use combined Z ordering within native windows; this is not a treatment-effect test or calibrated MotifTrace/MAPP parallel evidence. All native fits remain downloadable."})
        files.extend([preview_path, index_path, path])
    if not panels and not unavailable:
        return {"extra_views": [], "files": []}
    return {"extra_views": [{"project_id": project["project_id"], "source_project_manifest_sha256": file_sha256(root / "project.json"),
            "status": "complete", "panels": panels, "tables": tables, "unavailable": unavailable}], "files": files}
