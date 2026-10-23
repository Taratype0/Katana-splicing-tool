"""Bind completed native sequencing-QC plots and exact data to chapter one."""
from __future__ import annotations
import copy
from pathlib import Path


def build_qc_views(root, plan, receipts, transport):
    from .model import file_sha256, read_json, write_json
    from .runtime import canonical_hash
    root = Path(root).resolve()
    project = read_json(root / "project.json")
    translated = {item["source"]["path"]: root / item["collection_path"] for item in transport}
    panels, tables, source_proofs = [], [], []
    for task in plan["tasks"]:
        if task.get("adapter") != "qc.sequencing":
            continue
        _, receipt = receipts[task["id"]]
        if receipt.get("status") not in {"completed", "reused"}:
            raise ValueError("Sequencing QC views require a completed task")
        if receipt.get("receipt_sha256") != canonical_hash({k: v for k, v in receipt.items() if k != "receipt_sha256"}):
            raise ValueError("Sequencing QC completion proof changed")
        outputs = {item["path"]: item for item in receipt["outputs"]}

        def verified(original):
            original = str(original)
            if original not in outputs or original not in translated:
                raise ValueError("A sequencing QC view member was not transported from its completed task")
            path = translated[original]
            if file_sha256(path) != outputs[original]["sha256"]:
                raise ValueError("Sequencing QC scientific bytes changed after transport")
            return path

        evidence = receipt.get("execution_evidence", {})
        index = read_json(verified(evidence["qc_view_index"]))
        if index.get("schema") != "katana.flow.native_sequencing_qc/v1" or index.get("status") != "completed":
            raise ValueError("Unrecognised sequencing QC view provenance")
        for artifact in index["artifacts"]:
            path = verified(artifact["path"])
            phase = artifact.get("metadata", {}).get("qc_phase", "native")
            title = phase.replace("_", " ").capitalize() + " | " + artifact["name"]
            if path.suffix == ".png":
                metadata = copy.deepcopy(artifact.get("metadata", {}))
                links = []
                for link in metadata.get("data_links", []):
                    links.append({**link, "path": str(verified(link["path"]))})
                if {link.get("role") for link in links} != {"plot_data", "full_results", "methods"}:
                    raise ValueError("A native sequencing QC figure lost its exact data/method links")
                for export in metadata.get("exports", {}).values():
                    verified(export)
                panels.append({"id": "sequencing-" + phase + "-" + path.stem,
                    "chapter": 1, "title": title, "path": str(path), "contrast_id": "all",
                    "caption": metadata.get("caption", "Original native sequencing QC; no upstream rerun."),
                    "data_links": links, "metadata": {"qc_phase": phase, "native_QC_values_preserved": True}})
            elif path.suffix in {".tsv", ".csv", ".html"}:
                tables.append({"chapter": 1, "name": title, "path": str(path),
                               "analysis_id": "quality.reads"})
        for source in index["sources"]:
            path = verified(source["transported_native_path"])
            if file_sha256(path) != source["sha256"]:
                raise ValueError("Sequencing QC native source binding changed")
            tables.append({"chapter": 1, "kind": "file", "name": "Original native QC | " + source["biological_sample_id"] + " | " + path.name,
                           "path": str(path), "analysis_id": "quality.reads"})
            source_proofs.append({**source, "collection_path": path.relative_to(root).as_posix()})
    if not panels:
        return {"extra_views": [], "files": []}
    directory = root / "presentation" / "sequencing-qc"
    directory.mkdir(parents=True, exist_ok=True)
    proof = write_json(directory / "native-source-bindings.json", {"sources": source_proofs,
        "upstream_jobs": 0, "gene_body_coverage_scope": "separate_required_measured_workflow"})
    return {"extra_views": [{"project_id": project["project_id"],
        "source_project_manifest_sha256": file_sha256(root / "project.json"),
        "status": "complete", "panels": panels, "tables": tables,
        "audit": {"native_sources": str(proof), "upstream_jobs": 0}}], "files": [proof]}
