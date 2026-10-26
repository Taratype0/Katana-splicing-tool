"""Typed calculation and presentation checks for the approved complete report.

Native MAEI fitting and importing a previously calibrated parallel collection
are different contracts. A caption mentioning MAPP cannot establish either one.
All paths below are authenticated inside the portable collection; original
compute directories need not remain mounted when the report is generated.
"""
from __future__ import annotations

from pathlib import Path

from .model import file_sha256, read_json
from .runtime import canonical_hash


def _inside(root, relative):
    relative = Path(relative)
    path = (root / relative).resolve()
    if relative.is_absolute() or not path.is_relative_to(root):
        raise ValueError("Report calculation proof escapes its verified collection")
    return path


def verify_mapp_not_applicable(task, receipt, transported, *, expected_sample_ids=None):
    """Accept only the explicit, complete zero-exact-SE eligibility result.

    ``transported`` authenticates a member against this task's own output
    ledger. Missing references, failed fits and absent tables are never N/A.
    """
    import pandas as pd
    evidence = receipt.get("execution_evidence", {})
    if receipt.get("status") != "not_applicable" or evidence.get("native_fitting_performed") is not False:
        raise ValueError("MAPP non-applicability requires an explicit unfitted task receipt")
    proof = read_json(transported(evidence.get("applicability_proof", "")))
    if (proof.get("schema") != "katana.flow.mapp_applicability/v1" or
        proof.get("status") != "not_applicable" or proof.get("reason_code") != "no_exact_mapped_se" or
        proof.get("eligible_events") != 0 or proof.get("complete_measured_universe") is not True or
        proof.get("native_fitting_performed") is not False):
        raise ValueError("MAPP non-applicability lacks the complete zero-eligible measured universe")
    scope = task.get("scope", {})
    contrasts = scope.get("contrast_ids", [])
    if len(contrasts) != 1 or proof.get("comparison_id") != contrasts[0]:
        raise ValueError("MAPP non-applicability belongs to a different comparison")
    for key in ("site", "window_start", "window_end"):
        if key in scope and proof.get(key) != scope[key]:
            raise ValueError("MAPP non-applicability belongs to a different native window")
    if expected_sample_ids is not None and (len(proof.get("sample_ids", [])) != len(expected_sample_ids) or
                                            set(proof.get("sample_ids", [])) != set(expected_sample_ids)):
        raise ValueError("MAPP non-applicability belongs to a different biological sample design")
    table = transported(proof.get("eligibility_table", ""))
    if file_sha256(table) != proof.get("eligibility_table_sha256"):
        raise ValueError("MAPP full eligibility audit changed")
    frame = pd.read_csv(table, sep="\t", keep_default_na=False)
    allowed = {"no_exact_annotated_included_or_skipped_transcript", "outside_native_exon_length_range",
               "native_min_total_TPM_in_every_sample", "native_uninformative_exon_in_all_samples"}
    if (not {"event_id", "gene_id", "status"} <= set(frame) or
        len(frame) != proof.get("measured_events") or frame["event_id"].duplicated().any() or
        not set(frame["status"]) <= allowed):
        raise ValueError("MAPP zero-eligible proof does not match its complete measured event audit")
    return {**proof, "verified_eligibility_table": str(table)}


def _component(gap):
    metadata = gap.get("metadata") or {}
    component = gap.get("component") or metadata.get("component")
    if component:
        return "mapp_parallel" if component == "mapp" else component
    # Exact identifiers preserve compatibility with older explicit disclosures.
    # Arbitrary caption words, gene names and resource paths are never gates.
    ids = {"rich-mapp-not-computed": "mapp_parallel",
           "rich-mapp_parallel-not-computed": "mapp_parallel",
           "rich-network-bootstrap-not-completed": "network_bootstrap",
           "rich-chapter-6-not-computed": "network",
           "rich-chapter-7-not-computed": "motiftrace_sequence"}
    titles = {"Transcript 5′ / 3′ coverage bias": "gene_body_coverage",
              "Bootstrap module confidence": "network_bootstrap",
              "Original network dendrogram": "network_dendrogram"}
    return ids.get(gap.get("id"), titles.get(gap.get("title")))


def required_report_gaps(collection, plan, report_manifest, report_dir, *, collection_receipt=None):
    """Return concrete missing required calculations/views, never word matches.

    The collection validator checks every planned task and transported byte.
    Additional native checks bind each fit to its exact comparison/window and
    to the actual activity panel. Calibrated parallel imports are conditional
    on a declared asset or explicit import task and remain separate evidence.
    """
    from .collect import validate_collection
    root, report_dir = Path(collection).resolve(), Path(report_dir).resolve()
    collection_receipt = collection_receipt or validate_collection(root)
    project = read_json(root / collection_receipt.get("project_manifest", "project.json"))
    proofs = {row["task_id"]: row for row in collection_receipt["task_receipts"]}
    receipts = {}
    for task in plan["tasks"]:
        record = proofs.get(task["id"])
        if record is None:
            raise ValueError("Report collection lacks a planned task receipt")
        path = _inside(root, record["collection_path"])
        value = read_json(path)
        if (file_sha256(path) != record["sha256"] or
            value.get("receipt_sha256") != canonical_hash({k: v for k, v in value.items() if k != "receipt_sha256"}) or
            value.get("plan_id") != plan["plan_id"] or value.get("task_id") != task["id"] or
            value.get("task_contract_sha256") != canonical_hash(task) or
            value.get("source_project_manifest_sha256") != plan["project_manifest_sha256"]):
            raise ValueError("Report task proof does not match its sealed calculation contract")
        receipts[task["id"]] = value
    transport = {row["source"]["path"]: row for row in read_json(root / "transport-provenance.json")["products"]}

    def transported(receipt, original):
        owned = next((item for item in receipt.get("outputs", []) if item["path"] == str(original)), None)
        copied = transport.get(str(original))
        if owned is None or copied is None or copied["source"]["sha256"] != owned["sha256"]:
            raise ValueError("Required result is not an unchanged transported member of its actual task")
        path = _inside(root, copied["collection_path"])
        if file_sha256(path) != owned["sha256"]:
            raise ValueError("Required native result changed after collection")
        return path

    gaps, native_families, native_not_applicable = [], {}, set()

    def missing(component, reason, **details):
        gaps.append({"component": component, "reason": reason, **details})

    from .selection import analysis_selection, component_enabled
    enabled = analysis_selection(plan.get("configuration", {}))
    native_tasks = [task for task in plan["tasks"] if task.get("adapter") == "mapp.maei"]
    if enabled["mapp"] and not native_tasks:
        missing("mapp_native", "The complete profile has no planned native MAEI calculation")
    for task in native_tasks:
        receipt = receipts[task["id"]]
        try:
            if receipt.get("status") == "not_applicable":
                ids = task.get("scope", {}).get("contrast_ids", [])
                contrast = next((c for c in project.get("contrasts", []) if ids == [c["contrast_id"]]), None)
                if contrast is None:
                    raise ValueError("MAPP non-applicability has no exact declared contrast")
                sample_ids = contrast["numerator"]["sample_ids"] + contrast["denominator"]["sample_ids"]
                verify_mapp_not_applicable(task, receipt, lambda original: transported(receipt, original),
                                          expected_sample_ids=sample_ids)
                native_not_applicable.add(task["id"])
                continue
            if (receipt.get("status") not in {"completed", "reused"} or
                receipt.get("execution_evidence", {}).get("native_fitting_performed") is not True):
                raise ValueError("Native MAEI fitting has not actually completed")
            members = [item for item in receipt.get("outputs", [])
                       if Path(item["path"]).name == "AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv"]
            if len(members) != 1:
                raise ValueError("Native MAEI fitting lacks its unique unchanged native significance table")
            native_path = transported(receipt, members[0]["path"])
            import numpy as np
            import pandas as pd
            frame = pd.read_csv(native_path, sep="\t", index_col=0)
            activities = [column for column in frame if str(column).startswith("activities_")]
            if (frame.empty or "combined.Zscore" not in frame or not activities or
                not np.isfinite(frame[activities].apply(pd.to_numeric, errors="coerce").to_numpy()).any() or
                not pd.to_numeric(frame["combined.Zscore"], errors="coerce").notna().any()):
                raise ValueError("Native MAEI result lacks actual sample activities/Z statistics")
            contrasts = task.get("scope", {}).get("contrast_ids", [])
            if len(contrasts) != 1:
                raise ValueError("Each native MAEI report fit must have one exact comparison")
            native_families.setdefault(contrasts[0], {})[task["id"]] = (task, receipt, native_path)
        except (ValueError, OSError, KeyError) as exc:
            missing("mapp_native", str(exc), task_id=task["id"])

    for cid, entries in native_families.items():
        panels = [panel for panel in report_manifest.get("panels", [])
                  if panel.get("metadata", {}).get("component") == "mapp_native" and panel.get("contrast_id") == cid]
        try:
            if len(panels) != 1:
                raise ValueError("Every native MAEI comparison requires its actual sample-activity panel")
            panel, expected = panels[0], set(entries)
            metadata = panel["metadata"]
            if (set(metadata.get("fit_task_ids", [])) != expected or
                metadata.get("fit_receipt_sha256") != {key: value[1]["receipt_sha256"] for key, value in entries.items()} or
                metadata.get("calibrated_motiftrace_parallel") is not False):
                raise ValueError("Native activity panel does not bind every completed window fit")
            image = _inside(report_dir, panel.get("image", ""))
            if not image.is_file() or file_sha256(image) != metadata.get("native_activity_image_sha256"):
                raise ValueError("Native MAEI activity panel is absent or changed")
            links = panel.get("data_links", [])
            for role, key in (("plot_data", "native_plot_data_sha256"), ("full_results", "native_fit_index_sha256")):
                matches = [link for link in links if link.get("role") == role and
                           link.get("source_sha256") == metadata.get(key)]
                if len(matches) != 1 or file_sha256(_inside(report_dir, matches[0]["path"])) != metadata.get(key):
                    raise ValueError("Native MAEI panel lacks its exact data and complete native-fit index")
                if role == "full_results":
                    index = read_json(_inside(report_dir, matches[0]["path"]))
                    rows = index.get("fits", [])
                    if (index.get("comparison_id") != cid or index.get("calibrated_motiftrace_parallel") is not False or
                        len(rows) != len(expected) or {row["task_id"] for row in rows} != expected):
                        raise ValueError("Native-fit index does not retain every completed fit")
                    for row in rows:
                        task, receipt, native_path = entries[row["task_id"]]
                        if (row.get("sha256") != file_sha256(native_path) or row.get("scope") != task.get("scope") or
                            _inside(root, row["path"]) != native_path):
                            raise ValueError("Native-fit index changed the numerical family or window identity")
            if report_manifest.get("specification_version") == "1.1":
                from katana_next.reports.mapp_catalogs import validate_portable_mapp_panel
                validate_portable_mapp_panel(panel, report_root=report_dir)
        except (ValueError, OSError, KeyError) as exc:
            missing("mapp_native", str(exc), contrast_id=cid)

    parallel_assets = [asset for asset in project.get("assets", []) if asset.get("kind") == "motiftrace_parallel_results"]
    requested_id = plan.get("configuration", {}).get("parameters", {}).get("prioritize.motiftrace_parallel", {}).get("parallel_bundle_asset_id")
    parallel_tasks = [task for task in plan["tasks"] if task.get("workflow_id") == "prioritize.motiftrace_parallel"]
    parallel_required = enabled["mapp"] and bool(parallel_assets or requested_id or parallel_tasks)
    if parallel_required and not parallel_tasks:
        missing("mapp_parallel", "An explicitly supplied/requested calibrated parallel collection has no planned import")
    for task in parallel_tasks:
        receipt = receipts[task["id"]]
        try:
            if receipt.get("status") not in {"completed", "reused"} or receipt.get("job", {}).get("status") not in {"complete", "cache_hit"}:
                raise ValueError("The explicitly requested calibrated parallel import has not actually completed")
            browsers = [item for item in receipt.get("generated_artifacts", [])
                        if item.get("analysis_id") == "prioritize.motiftrace_parallel" and
                        item.get("entity_type") == "motiftrace_parallel_report"]
            if len(browsers) != 1:
                raise ValueError("The requested calibrated parallel collection lacks its verified result browser")
            transported(receipt, browsers[0]["path"])
        except (ValueError, OSError, KeyError) as exc:
            missing("mapp_parallel", str(exc), task_id=task["id"])

    selection_components = {"clip": "prioritize.rbp_mechanism", "motif_enrichment": "prioritize.motiftrace_consequences",
                            "motiftrace_sequence": "prioritize.motiftrace_consequences"}
    required = {"gene_body_coverage", "network", "network_bootstrap", "network_dendrogram",
                "rbp_identity", "motif_enrichment", "motiftrace_sequence", "clip", "mapp_native"}
    from .presentation import network_bootstrap_enabled
    if not network_bootstrap_enabled(plan):
        required.remove("network_bootstrap")
    interactions_required = bool(plan.get("configuration", {}).get("interactions") or
                                 len(plan.get("configuration", {}).get("design_factors", [])) >= 2)
    if interactions_required:
        required.add("interaction")
    if parallel_required:
        required.add("mapp_parallel")
    for gap in report_manifest.get("unavailable", []):
        component = _component(gap)
        if not component_enabled(component, enabled):
            continue
        if component == "mapp_native" and (gap.get("task_id") in native_not_applicable or
                                            len(native_not_applicable) == len(native_tasks)):
            continue
        workflow = selection_components.get(component)
        affected = [task for task in plan["tasks"] if task.get("workflow_id") == workflow] if workflow else []
        if affected and all(receipts[task["id"]].get("status") == "not_applicable" for task in affected):
            continue
        if component in required or gap.get("required") is True:
            gaps.append({**gap, "component": component})
    from .network_acceptance import network_report_gaps
    if enabled["network"]:
        gaps.extend(network_report_gaps(root, plan, report_manifest, report_dir))
    from .motif_acceptance import required_motif_report_gaps
    if enabled["motiftrace"]:
        gaps.extend(required_motif_report_gaps(
            root, plan, report_manifest, report_dir,
            project=project, receipts=receipts, transported=transported,
        ))
    return gaps
