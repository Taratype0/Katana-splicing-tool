"""Import a genuine upstream bundle using explicit IDs and sample-vector order."""
from __future__ import annotations

import copy
import csv
import io
import os
import shutil
from pathlib import Path, PurePosixPath

from .model import file_sha256, read_json, write_json


def _upstream_metadata(source, paths):
    """Recover declarations only from the ledger-verified source map.

    Older exporters retained this file but omitted its declarations from their
    manifest. Declared labels remain declarations, never inferred validation.
    """
    source_map = paths.get("provenance/source_map.json")
    recorded = read_json(source_map) if source_map else {}
    result = {}
    for key in ("reference", "reference_declarations", "gene_exclusions", "assay",
                "assay_validation", "design", "quantification_method"):
        candidates = [item[key] for item in (source, recorded) if item.get(key) not in (None, {}, "")]
        if candidates and any(candidate != candidates[0] for candidate in candidates[1:]):
            raise ValueError(f"Conflicting upstream {key} declarations")
        if candidates:
            result[key] = copy.deepcopy(candidates[0])
    for key in ("reference", "reference_declarations", "gene_exclusions"):
        if key in result and not isinstance(result[key], dict):
            raise ValueError(f"Invalid upstream {key} declarations")
    builds = {str(value["genome_build"]) for key in ("reference", "reference_declarations")
              if isinstance(value := result.get(key), dict) and value.get("genome_build")}
    if len(builds) > 1:
        raise ValueError("Conflicting upstream genome_build declarations")
    return result, next(iter(builds), "")


def import_bundle(bundle, project_template, outdir, *, owned_link_root=None, analysis_configuration=None):
    from katana_next.domain.models import ProjectManifest
    from .selection import resource_kinds
    selected_resources = resource_kinds(analysis_configuration or {})
    optional_resources = resource_kinds({})
    skipped_resources = []
    bundle = Path(bundle).resolve()
    template_path = Path(project_template).resolve()
    template = read_json(template_path)
    source = read_json(bundle / "manifest.json")
    if source.get("schema") != "katana-pipeline-bundle/v1" or source.get("status") != "complete":
        raise ValueError("Complete mode requires a complete genuine upstream bundle")
    out = Path(outdir).resolve()
    # Managed handoffs may link their own immutable checkpoint files within one
    # filesystem. Borrowed user inputs retain the independent-copy behavior.
    owned = Path(owned_link_root).resolve(strict=True) if owned_link_root is not None else None
    if owned is not None and (not bundle.is_relative_to(owned) or not out.is_relative_to(owned)):
        raise ValueError("Linked bundle handoff must remain inside its owned stage")
    if out.exists():
        raise ValueError("Import destination already exists; completed inputs are preserved")
    artifacts = source.get("artifacts", [])
    if not artifacts:
        raise ValueError("Bundle has no native artifact ledger")
    paths = {}
    for item in artifacts:
        if not str(item.get("status", "")).startswith("available"):
            continue
        relative = PurePosixPath(item["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or ":" in str(relative):
            raise ValueError("Unsafe upstream bundle member")
        path = bundle / relative
        if not path.resolve().is_relative_to(bundle) or not path.is_file():
            raise ValueError(f"Missing or escaping native bundle member: {relative}")
        if path.stat().st_size != item.get("bytes") or file_sha256(path) != item.get("sha256"):
            raise ValueError(f"Upstream bundle member failed SHA256 validation: {relative}")
        if str(relative) in paths:
            raise ValueError("Duplicate upstream bundle member")
        paths[str(relative)] = path
    upstream_metadata, genome_build = _upstream_metadata(source, paths)
    project_metadata = copy.deepcopy(template.get("metadata", {}))
    releases = {value["annotation_release"] for key in ("reference", "reference_declarations")
                if isinstance(value := upstream_metadata.get(key), dict) and value.get("annotation_release")}
    if len(releases) > 1:
        raise ValueError("Conflicting upstream annotation_release declarations")
    if releases:
        release = next(iter(releases))
        if project_metadata.get("annotation_release") not in (None, "", release):
            raise ValueError("Template annotation_release conflicts with verified upstream declaration")
        project_metadata["annotation_release"] = release
    groups, samples, sample_ids = {}, [], set()
    for row in source["samples"]:
        sid, group = row["sample_id"], row["condition"]
        if sid in sample_ids:
            raise ValueError("Duplicate biological sample ID in upstream bundle")
        sample_ids.add(sid)
        groups.setdefault(group, []).append(sid)
        attributes = copy.deepcopy(template.get("sample_attributes", {}).get(sid, {}))
        for key in ("source", "assay"):
            if row.get(key):
                if key in attributes and attributes[key] != row[key]:
                    raise ValueError(f"Template {key} conflicts with verified upstream sample {sid}")
                attributes[key] = row[key]
        samples.append({"sample_id": sid, "biological_sample_id": row.get("subject") or sid,
                        "condition": group, "batch": row.get("batch") or None, "attributes": attributes, "assets": {}})
    by_sample = {s["sample_id"]: s for s in samples}
    comparisons, contrast_index = [], {}
    project_metadata["comparison_designs"] = {}
    for row in source["comparisons"]:
        cid, numerator, denominator = row["comparison_id"], row["treatment"], row["control"]
        if cid in contrast_index or numerator not in groups or denominator not in groups:
            raise ValueError("Invalid declared comparison in upstream bundle")
        comparison = {"contrast_id": cid, "label": template.get("contrast_labels", {}).get(cid, cid),
                      "numerator": {"label": numerator, "sample_ids": groups[numerator]},
                      "denominator": {"label": denominator, "sample_ids": groups[denominator]}}
        comparisons.append(comparison)
        contrast_index[cid] = row
        project_metadata["comparison_designs"][cid] = {
            "paired": row.get("paired", False), "formula": row.get("design", "~ condition"),
            "blocking_column": "biological_sample_id" if row.get("paired", False) else None}
    annotation, assets = {"genome_build": genome_build}, []
    normalized_members = {}
    stage_r_contrasts = {item.get("comparison_id") for item in artifacts
                         if item.get("kind") == "stager" and item.get("relative_path") in paths}
    chosen_mode = template.get("rmats_mode", "JCEC")
    if chosen_mode not in {"JC", "JCEC"}:
        raise ValueError("Explicit JC or JCEC mode is required")
    for index, item in enumerate(artifacts):
        relative = item.get("relative_path")
        if relative not in paths:
            continue
        kind = item["kind"]
        if kind in {"gtf", "genome_fasta", "transcript_fasta"}:
            annotation[kind] = relative
        if kind in {"bam", "bai", "quant"} and item.get("sample_id") in by_sample:
            by_sample[item["sample_id"]]["assets"][{"quant": "quant_sf"}.get(kind, kind)] = relative
        mapped = {"rmats": "rmats_result", "deg": "deg_result", "dtu": "dtu_result", "stager": "dtu_result", "dexseq": "dexseq_result",
                  "normalized_counts": "normalized_counts", "raw_counts": "raw_counts", "qc": "sequencing_qc"}.get(kind)
        if kind == "dtu" and item.get("comparison_id") in stage_r_contrasts:
            # Keep the underlying DEXSeq transcript test as a separate source.
            # The completed stageR transcript-confirmation table is the DTU
            # inference; counting both would duplicate each transcript.
            mapped = "dtu_native_test_result"
        native_exon = (kind == "auxiliary" and relative.startswith("dexseq_exon/")
                       and PurePosixPath(relative).name.startswith("DEXSeqResults."))
        projection = None
        asset_relative = relative
        if kind == "dexseq" or native_exon:
            # The upstream ledger historically called perGeneQValue 'dexseq'
            # and kept the actual exon-bin results as an auxiliary native CSV.
            # Preserve both native files; only exon rows enter exon inference.
            text = paths[relative].read_text(encoding="utf-8-sig")
            if not text.strip():
                raise ValueError("Empty native DEXSeq table: " + relative)
            delimiter = "\t" if "\t" in text.splitlines()[0] else ","
            rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
            if (not rows or not {"groupID", "padj"}.issubset(rows[0]) or
                    len(rows[0]) != len(set(rows[0])) or any(len(row) != len(rows[0]) for row in rows[1:])):
                raise ValueError("Invalid native DEXSeq table: " + relative)
            exon_level = "featureID" in rows[0]
            if native_exon and not exon_level:
                raise ValueError("Native exon results lack featureID: " + relative)
            mapped = "dexseq_result" if exon_level else "dexseq_gene_result"
            projection = {"native_path": relative, "native_sha256": item["sha256"],
                          "unit": "exon_bin" if exon_level else "gene",
                          "native_delimiter": delimiter, "native_statistics_recomputed": False}
            if exon_level and delimiter != "\t":
                asset_relative = f"adapters/native-{index}.dexseq-exons.tsv"
                buffer = io.StringIO(newline="")
                csv.writer(buffer, delimiter="\t", lineterminator="\n").writerows(rows)
                normalized_members[asset_relative] = buffer.getvalue().encode("utf-8")
        if not mapped:
            continue
        cid = item.get("comparison_id")
        metadata = {"upstream_sha256": item["sha256"], "upstream_artifact": copy.deepcopy(item)}
        if projection:
            metadata["native_table_projection"] = projection
            if item.get("native_fold_change_direction"):
                metadata["native_fold_change_direction"] = item["native_fold_change_direction"]
        if kind == "qc":
            metadata.update({key: item[key] for key in ("sample_id", "mate", "method", "qc_type", "qc_phase") if key in item})
        if cid:
            if cid not in contrast_index:
                raise ValueError("Native result has an undeclared comparison")
            row = contrast_index[cid]
            metadata["contrast_id"] = cid
            if kind == "rmats":
                mode = item.get("mode") or ("JCEC" if ".JCEC." in relative else "JC")
                if mode != chosen_mode:
                    continue
                left, right = row["rmats_group1"], row["rmats_group2"]
                if {left, right} != {row["treatment"], row["control"]}:
                    raise ValueError("Native rMATS direction is incompatible with declared contrast")
                for group in (left, right):
                    order = source.get("rmats_sample_order", {}).get(group)
                    if order is None or len(order) != len(set(order)) or set(order) != set(groups[group]):
                        raise ValueError("Native rMATS vector order must cover exactly the declared biological samples")
                metadata.update(event_type=item.get("event_type") or PurePosixPath(relative).name.split(".")[0],
                                mode=mode, reverse_direction=left != row["treatment"],
                                numerator_sample_ids=source["rmats_sample_order"][row["treatment"]],
                                denominator_sample_ids=source["rmats_sample_order"][row["control"]])
            elif kind == "deg":
                metadata["reverse_direction"] = row["deg_target"] != row["treatment"]
        assets.append({"asset_id": f"native-{index}", "kind": mapped, "path": asset_relative, "required": True, "metadata": metadata})
    for asset in template.get("assets", []):
        if asset.get("kind") in optional_resources - selected_resources:
            skipped_resources.append(asset["asset_id"])
            continue
        asset = copy.deepcopy(asset)
        path = Path(asset["path"])
        path = path if path.is_absolute() else template_path.parent / path
        if not path.is_file():
            raise ValueError(f"Required analysis resource unavailable: {asset['asset_id']}")
        relative = "resources/" + str(len(assets)) + "/" + path.name
        asset["path"] = relative
        paths[relative] = path.resolve()
        assets.append(asset)
    payload = {"schema_version": "2.0", "project_id": template.get("project_id") or "bundle-" + file_sha256(bundle / "manifest.json")[:20],
               "name": template.get("name", source.get("project_name", "Katana study")),
               "input_root": "bundle", "output_root": "results", "annotation": annotation,
               "samples": samples, "contrasts": comparisons, "assets": assets,
               "metadata": {**project_metadata, "upstream_bundle_sha256": file_sha256(bundle / "manifest.json"),
                            "native_statistics_preserved": True, "rmats_mode": chosen_mode,
                            "upstream_metadata": upstream_metadata}}
    ProjectManifest.model_validate(payload)
    out.mkdir(parents=True)
    # A single server-side copy creates a relocatable project. BAMs never need
    # downloading to the reporting machine; they stay in the execution project.
    for relative, path in paths.items():
        destination = out / "bundle" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if owned is not None and path.resolve().is_relative_to(owned) and relative not in normalized_members:
            try:
                os.link(path, destination)
            except OSError:
                shutil.copyfile(path, destination)
        else:
            shutil.copyfile(path, destination)
        if file_sha256(path) != file_sha256(destination):
            raise ValueError("Imported native bytes changed while copying")
    for relative, content in normalized_members.items():
        destination = out / "bundle" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    write_json(out / "project.json", payload)
    write_json(out / "import.receipt.json", {"schema": "katana.flow.import/v1", "status": "completed",
               "source_manifest_sha256": file_sha256(bundle / "manifest.json"), "samples": len(samples),
               "comparisons": len(comparisons), "files_verified": len(paths), "mode": chosen_mode,
               "source_statistics_refitted": False,
               "not_selected_resource_ids": skipped_resources,
               "format_projections": [{"path": relative, "sha256": file_sha256(out / "bundle" / relative)}
                                      for relative in normalized_members]})
    return out / "project.json"
