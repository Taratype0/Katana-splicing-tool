"""Positive per-network checks for the fixed report, including silent omissions.

This is a read-only acceptance gate. A missing panel cannot be hidden by simply
omitting its unavailable entry. Empty statistical families are authenticated
against the completed native products, never replaced with invented results.
"""
from pathlib import Path

from .model import read_json, file_sha256


def _inside(root, value):
    value = Path(value)
    path = (root / value).resolve()
    if value.is_absolute() or not path.is_relative_to(root):
        raise ValueError("Network report proof escapes its collection or report")
    return path


def _role(metadata):
    if metadata.get("network_view_role"):
        return metadata["network_view_role"]
    if metadata.get("flow_stage") == "network.fit" and metadata.get("flow_artifact_key") == "dendrogram":
        return "dendrogram"
    if metadata.get("flow_stage") == "network.bootstrap" and metadata.get("flow_artifact_key") == "confidence":
        return "bootstrap_confidence"
    return None


def network_report_gaps(collection, plan, report_manifest, report_dir):
    """Require each planned network's own completed figures and module browser."""
    import pandas as pd
    from .presentation import _network_task_identity, validate_network_view_proof

    root, report_dir = Path(collection).resolve(), Path(report_dir).resolve()
    fits = [task for task in plan["tasks"] if task.get("adapter") == "network.fit"]
    gaps = []

    def missing(reason, *, task=None, role=None):
        gaps.append({"component": "network", "reason": reason,
                     **({"task_id": task["id"], "contrast_ids": task.get("scope", {}).get("contrast_ids", [])} if task else {}),
                     **({"network_view_role": role} if role else {})})

    if not fits:
        missing("The fixed report has no planned network fit")
        return gaps
    try:
        project = read_json(root / "project.json")
        descriptor = read_json(root / "report-project.json")
        exported = read_json(root / "artifacts.json")["artifacts"]
        assets = [*exported, *descriptor.get("extra_assets", [])]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        missing("Network collection catalogue is absent or invalid: " + str(exc))
        return gaps
    expected_contrasts = plan.get("contrast_ids", [c["contrast_id"] for c in project.get("contrasts", [])])
    scope = plan.get("configuration", {}).get("network", {}).get("scope", "contrast")
    if scope == "contrast":
        actual = [task.get("parameters", {}).get("contrast_ids") for task in fits]
        if any(not isinstance(value, list) for value in actual) or sorted(actual) != sorted([[cid] for cid in expected_contrasts]):
            missing("The selected comparisons do not each have one independently planned network fit")

    for fit in fits:
        try:
            identity = _network_task_identity(fit, project)
            if identity is None:
                raise ValueError("A complete generic flow network requires an explicit sample/comparison identity")
            cid = identity["contrast_ids"][0] if identity["network_scope"] == "contrast" else "all"
            scoped = [a for a in assets if a.get("metadata", {}).get("network_identity") == identity]
            proof_names = {a.get("metadata", {}).get("completion_manifest") for a in scoped
                           if a.get("metadata", {}).get("completion_scope") == "flow_verified_network_complete"}
            proof_names.discard(None)
            if len(proof_names) != 1:
                raise ValueError("This network lacks its unique completed four-stage presentation proof")
            proof_path = _inside(root, next(iter(proof_names)))
            proof = read_json(proof_path)
            chosen = [c for c in project["contrasts"] if c["contrast_id"] in identity["contrast_ids"]]
            chosen.sort(key=lambda row: identity["contrast_ids"].index(row["contrast_id"]))
            _, _, precomputed, _ = validate_network_view_proof(proof_path, [], chosen)
            if precomputed["provenance"]["source_task_id"] != fit["id"]:
                raise ValueError("Network views are bound to another fit task")
            score_files = [member for member in proof["files"]
                           if member["task_id"] == fit["id"]
                           and Path(member["source_path"]).name == "module_program_scores.tsv"]
            if len(score_files) != 1:
                raise ValueError("Network proof lacks its one exact complete module activity matrix")
            scores = pd.read_csv(_inside(root, score_files[0]["path"]), sep="\t", index_col=0)
            selected = pd.read_csv(precomputed["selected_tests"], sep="\t")
            tested = pd.read_csv(precomputed["tested"], sep="\t")
            required = {"all_module_sizes", "all_module_event_types", "dendrogram", "bootstrap_confidence", "module_browser"}
            if len(scores.columns):
                required.add("all_module_activity")
            if len(selected):
                required.add("matched_module_effects")
            if len(tested):
                required.add("module_enrichment_overview")
            # Zero modules/zero testable terms are valid recorded outcomes. The
            # required overview/browser and complete tables still expose them.
            for role in sorted(required):
                candidates = [a for a in scoped if _role(a.get("metadata", {})) == role]
                if len(candidates) != 1:
                    raise ValueError(f"Network lacks one exact {role} source asset")
                source = candidates[0]
                metadata = source.get("metadata", {})
                source_path = _inside(root, source["path"])
                if not source_path.is_file():
                    raise ValueError(f"Network {role} source file is absent")
                records = report_manifest.get("complete_results", []) if role == "module_browser" else report_manifest.get("panels", [])
                matches = [r for r in records if r.get("chapter") == 6 and r.get("contrast_id") == cid
                           and not r.get("output_only")
                           and (r.get("metadata", {}).get("source_artifact_id") or r.get("source_artifact_id")) == source["artifact_id"]]
                if len(matches) != 1:
                    raise ValueError(f"Network {role} is absent or duplicated in its exact comparison report")
                rendered = matches[0]
                rmeta = rendered.get("metadata", {})
                if rmeta.get("network_identity") != identity or _role(rmeta) != role:
                    raise ValueError(f"Network {role} report identity differs from its completed fit")
                target = _inside(report_dir, rendered["source" if role == "module_browser" else "image"])
                if role == "module_browser":
                    if rmeta.get("source_sha256") != file_sha256(source_path) or file_sha256(target) != rmeta.get("portable_sha256"):
                        raise ValueError("Network browser source or portable bytes changed")
                else:
                    if file_sha256(target) != file_sha256(source_path):
                        raise ValueError(f"Network {role} report image differs from its completed source")
                    plot_sources = [link for link in metadata.get("data_links", []) if link.get("role") == "plot_data"]
                    plotted = [link for link in rendered.get("data_links", []) if link.get("role") == "plot_data"]
                    if len(plot_sources) != 1 or len(plotted) != 1:
                        raise ValueError(f"Network {role} lacks its exact plotted data")
                    original_data = _inside(root, plot_sources[0]["path"])
                    digest = file_sha256(original_data)
                    if (plotted[0].get("contrast_id") != cid or plotted[0].get("source_sha256") != digest
                            or file_sha256(_inside(report_dir, plotted[0]["path"])) != digest):
                        raise ValueError(f"Network {role} plotted data changed or belongs to another comparison")
                    if role == "module_enrichment_overview":
                        displayed = pd.read_csv(original_data, sep="\t")
                        if not pd.to_numeric(displayed["FDR"], errors="coerce").between(0, .05, inclusive="both").all():
                            raise ValueError("Network report enrichment includes nonsignificant or invalid FDR terms")
        except (ValueError, OSError, KeyError, TypeError) as exc:
            missing(str(exc), task=fit)
    return gaps
