"""Reuse the native MultiQC report and its machine-readable measured data."""
import json
from pathlib import Path
import re
import shutil


def completed_multiqc(manifest):
    assets = [a for a in manifest.assets if a.kind == "sequencing_qc" and
              (a.metadata.get("qc_type") == "multiqc" or a.path.name in {"multiqc_report.html", "multiqc_data.json"})]
    if not assets:
        return None
    html = [a for a in assets if a.path.suffix.lower() == ".html"]
    data = [a for a in assets if a.path.name == "multiqc_data.json"]
    if len(html) != 1 or len(data) != 1:
        raise ValueError("Supply one matched MultiQC HTML report and multiqc_data.json for this study")
    document = json.loads(data[0].path.read_text(encoding="utf-8-sig"))
    report_text = html[0].path.read_text(encoding="utf-8-sig")
    if "multiqc" not in report_text.lower() or "<html" not in report_text.lower():
        raise ValueError("Invalid MultiQC HTML report")
    # Original module sample names are retained, including read mate suffixes.
    sample_names = set()
    modules = []
    for name, values in document.get("report_saved_raw_data", {}).items():
        if isinstance(values, dict) and values and name.lower() not in {"software_versions", "custom_content"}:
            sample_names.update(str(k) for k in values)
            modules.append(name)
    for values in document.get("report_general_stats_data", []):
        if isinstance(values, dict):
            sample_names.update(str(k) for k in values)
    declared = {s.sample_id for s in manifest.samples}
    mapping = {}
    for name in sample_names:
        candidates = [sid for sid in declared if name == sid or re.match(re.escape(sid) + r"[_. -]", name)]
        if candidates:
            longest = max(map(len, candidates))
            best = [sid for sid in candidates if len(sid) == longest]
            if len(best) == 1:
                mapping[name] = best[0]
    if not modules or not mapping:
        raise ValueError("MultiQC contains no measured QC mapped to this study; version-only reports are insufficient")
    missing = sorted(declared - set(mapping.values()))
    if missing:
        raise ValueError("MultiQC is missing declared samples: " + ", ".join(missing))
    return html[0], data[0], {"sample_mapping": mapping, "modules": sorted(modules),
        "multiqc_version": document.get("config_version", document.get("multiqc_version", "recorded in native report")),
        "coverage_available": any("qualimap" in x.lower() or "gene_body" in x.lower() for x in modules),
        "raw_reads_or_bams_read": False, "upstream_jobs": 0}


def reuse_multiqc(manifest, output_dir):
    from .model import file_sha256, write_json
    selected = completed_multiqc(manifest)
    if selected is None:
        return None
    html, data, audit = selected
    out = Path(output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Existing sequencing QC outputs are preserved")
    out.mkdir(parents=True, exist_ok=True)
    sources = []
    for asset, name in ((html, "multiqc_report.html"), (data, "multiqc_data.json")):
        target = out / name
        shutil.copyfile(asset.path, target)
        digest = file_sha256(target)
        if digest != file_sha256(asset.path):
            raise ValueError("MultiQC source changed during transport")
        sources.append({"source_path": str(asset.path), "transported_native_path": str(target),
                        "sha256": digest, "biological_sample_id": "all samples"})
    proof = write_json(out / "multiqc-reuse.json", audit)
    links = [{"role": role, "label": label, "path": str(path), "contrast_id": "all",
              "policy": None, "unit": "native MultiQC sample/module"}
             for role, label, path in (("plot_data", "Original MultiQC data", out / "multiqc_data.json"),
                                       ("full_results", "Complete MultiQC numerical results", out / "multiqc_data.json"),
                                       ("methods", "QC source and sample mapping", proof))]
    index = write_json(out / "sequencing-qc-views.json", {
        "schema": "katana.flow.native_sequencing_qc/v1", "status": "completed",
        "sources": sources, "upstream_jobs": 0,
        "artifacts": [{"name": "MultiQC sequencing quality", "path": str(out / "multiqc_report.html"),
                       "metadata": {"native_multiqc": True, "qc_phase": "multiqc", "data_links": links,
                                    "caption": "Original MultiQC results. Only measured modules are shown; no FASTQ or BAM scan was repeated."}}]})
    outputs = [{"path": str(p), "size_bytes": p.stat().st_size, "sha256": file_sha256(p),
                "kind": "file", "name": p.name, "analysis_id": "flow.qc.sequencing"}
               for p in out.iterdir() if p.is_file()]
    return {"schema": "katana.flow.stage/v1", "adapter": "qc.sequencing", "status": "complete",
            "outputs": outputs, "evidence": {**audit, "qc_view_index": str(index), "native_multiqc": True}}
