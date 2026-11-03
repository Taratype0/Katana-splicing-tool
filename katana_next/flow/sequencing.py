"""Transport and display native sequencing QC; never launch upstream tools."""
from __future__ import annotations

import json
import csv
import math
from pathlib import Path
import re
import shutil

from .model import file_sha256, read_json, write_json

BAM_SOURCES = {"bam", "genome_bam", "transcriptome_bam"}


def _project_sources(spec):
    from katana_next.infrastructure.manifest_store import ManifestStore
    path = spec.get("inputs", {}).get("project")
    if not path:
        raise ValueError("Sequencing QC requires the exact verified project manifest")
    _, manifest = ManifestStore().load(path)
    sources = [asset for asset in manifest.assets if asset.kind == "sequencing_qc"]
    if not sources:
        raise ValueError("No completed native sequencing QC assets were imported into this project")
    return manifest, sources


def _source_type(path, metadata):
    name = Path(path).name
    if name == "strandedness.summary.tsv":
        if metadata.get("qc_type", "strandedness") != "strandedness" or metadata.get("qc_phase", "library") != "library":
            raise ValueError("Conflicting native strandedness QC declaration")
        return "strandedness", "library"
    if name.endswith(("_fastqc.zip", "_fastqc.html")):
        phase = metadata.get("qc_phase")
        if phase not in {"raw", "trimmed", "bam_records"}:
            phase = "trimmed" if re.search(r"_val_[12]_fastqc|_trimmed_fastqc", name) else "raw"
        return "fastqc", phase
    if name.endswith("_trimming_report.txt"):
        return "trim_galore", "trimming"
    if name.endswith("Log.final.out"):
        return "star", "alignment"
    raise ValueError("Unsupported native sequencing QC file: " + name)


def _sample(path, metadata, sample_ids):
    declared = metadata.get("sample_id")
    if declared:
        if declared not in sample_ids:
            raise ValueError("Sequencing QC refers to an undeclared biological sample: " + str(declared))
        return declared
    name = Path(path).name
    matches = [sid for sid in sample_ids if name == sid or name.startswith(sid + "_")
               or name.startswith(sid + ".") or name.startswith(sid + "Log.final.out")]
    if not matches:
        raise ValueError("Native sequencing QC has no unambiguous sample identity: " + name)
    longest = max(map(len, matches))
    matches = [sid for sid in matches if len(sid) == longest]
    if len(matches) != 1:
        raise ValueError("Ambiguous sequencing QC sample mapping: " + name)
    return matches[0]


def _source_samples(path, metadata, sample_ids, kind):
    if kind != "strandedness":
        return {_sample(path, metadata, sample_ids)}
    columns = {"sample", "requested", "strandedness", "method", "informative_forward_fraction",
               "informative_reverse_fraction", "failed_to_determine_fraction", "usable_reads_sampled"}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)) or not columns <= set(reader.fieldnames):
            raise ValueError("Invalid native strandedness summary columns")
        rows = list(reader)
    seen = set()
    for row in rows:
        sid = row.get("sample")
        if None in row or any(row.get(key) is None for key in columns):
            raise ValueError("Malformed native strandedness summary row")
        if sid not in sample_ids or sid in seen:
            raise ValueError("Unknown or duplicate sample in native strandedness summary: " + str(sid))
        seen.add(sid)
        resolved = row["strandedness"]
        if resolved not in {"forward", "reverse", "unstranded"}:
            raise ValueError("Unresolved native strandedness for " + sid)
        if row["method"] == "user_declared_not_inferred":
            if row["requested"] != resolved:
                raise ValueError("Conflicting declared strandedness for " + sid)
            if any(row[key] for key in columns - {"sample", "requested", "strandedness", "method"}):
                raise ValueError("Declared strandedness must not claim inferred measurements: " + sid)
        elif row["method"] == "rseqc_infer_experiment" and row["requested"] == "auto":
            for key in ("informative_forward_fraction", "informative_reverse_fraction", "failed_to_determine_fraction"):
                try:
                    value = float(row[key])
                except ValueError as error:
                    raise ValueError("Invalid native strandedness fraction for " + sid) from error
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("Invalid native strandedness fraction for " + sid)
            if not row["usable_reads_sampled"].isdigit() or int(row["usable_reads_sampled"]) <= 0:
                raise ValueError("Invalid native strandedness read count for " + sid)
        else:
            raise ValueError("Unsupported native strandedness method for " + sid)
    if seen != sample_ids:
        raise ValueError("Native strandedness summary is missing samples: " + ", ".join(sorted(sample_ids - seen)))
    return seen


def readiness(spec):
    if "project" in spec.get("readiness_deferred_inputs", []):
        return []
    try:
        manifest, sources = _project_sources(spec)
        sample_ids = {sample.sample_id for sample in manifest.samples}
        # Genome-BAM input does not create fictional raw FASTQ or STAR QC.
        # Legacy projects with no source declaration use the complete cohort.
        expected_samples = {sample.sample_id for sample in manifest.samples
                            if sample.attributes.get("source", "").casefold() not in BAM_SOURCES}
        covered = {}
        for asset in sources:
            if not asset.path.is_file() or not asset.path.stat().st_size:
                raise ValueError("Missing or empty completed native QC asset: " + asset.asset_id)
            kind, phase = _source_type(asset.path, asset.metadata)
            covered_samples = _source_samples(asset.path, asset.metadata, sample_ids, kind)
            if not (kind == "fastqc" and asset.path.name.endswith(".html")):
                covered.setdefault(phase, set()).update(covered_samples)
            expected_hash = asset.metadata.get("upstream_sha256")
            if expected_hash and file_sha256(asset.path) != expected_hash:
                raise ValueError("Completed native QC source SHA256 changed: " + asset.asset_id)
        bam_samples = sample_ids - expected_samples
        for phase, samples in covered.items():
            eligible = sample_ids if phase == "library" else (bam_samples if phase == "bam_records" else expected_samples)
            missing = eligible - samples
            if missing:
                raise ValueError("Native sequencing QC phase " + phase + " is missing declared eligible samples: " + ", ".join(sorted(missing)))
        return []
    except (ValueError, OSError, KeyError) as error:
        return [str(error)]


def execute_stage(spec, output_dir):
    from katana_next.reports.qc import build_qc_assets, _identity
    issues = readiness(spec)
    if issues:
        raise ValueError("; ".join(issues))
    out = Path(output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Existing sequencing QC outputs are preserved")
    out.mkdir(parents=True, exist_ok=True)
    manifest, sources = _project_sources(spec)
    sample_ids = {sample.sample_id for sample in manifest.samples}
    phases, source_ledger, duplicate_keys = {}, [], {}
    for index, asset in enumerate(sources):
        kind, phase = _source_type(asset.path, asset.metadata)
        covered_samples = _source_samples(asset.path, asset.metadata, sample_ids, kind)
        sid = "all samples" if kind == "strandedness" else next(iter(covered_samples))
        before = asset.path.stat()
        digest = file_sha256(asset.path)
        destination = out / "native" / phase / f"{index:06d}" / asset.path.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(asset.path, destination)
        after = asset.path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns) or file_sha256(destination) != digest:
            raise ValueError("Native QC source changed during verified transport")
        if asset.metadata.get("upstream_sha256") and asset.metadata["upstream_sha256"] != digest:
            raise ValueError("Native QC source no longer matches upstream provenance")
        row = {"relative_path": destination.relative_to(out / "native" / phase).as_posix(),
               "recorded_size": before.st_size, "sha256": digest, "sample_id": sid,
               "biological_sample_id": sid, "asset_id": asset.asset_id,
               "native_id": asset.metadata.get("source_native_id"), "source_path": str(asset.path),
               "transported_native_path": str(destination), "qc_type": kind, "qc_phase": phase,
               "source_sample_ids": sorted(covered_samples)}
        if asset.path.name.endswith("_fastqc.zip"):
            mate = _identity(asset.path, sid)[1]
            duplicate_keys.setdefault((phase, sid, mate), []).append(row)
        phases.setdefault(phase, []).append(row)
        source_ledger.append(row)
    # Distinct technical read files retain their true biological mapping while
    # a unique display label keeps the native FastQC pivot lossless.
    for rows in duplicate_keys.values():
        if len(rows) > 1:
            for row in rows:
                row["sample_id"] += " | " + row["asset_id"]
                row["technical_file_display_label"] = True
    artifacts, phase_receipts = [], {}
    for phase, rows in phases.items():
        native = out / "native" / phase
        write_json(native / "QC_FETCH_MANIFEST.json", {"files": rows})
        if phase == "library":
            # A cohort summary is a measured/declaration table, not a per-read
            # FastQC report. Preserve its exact bytes in the existing QC links.
            for row in rows:
                artifacts.append({"name": "Resolved library strandedness (all samples)",
                    "path": row["transported_native_path"], "kind": "table", "status": "complete",
                    "analysis_id": "quality.reads", "metadata": {"qc_phase": phase,
                    "original_native_QC": True, "native_QC_values_preserved": True,
                    "sample_ids": row["source_sample_ids"], "upstream_jobs": 0}})
            phase_receipts[phase] = {"status": "completed", "source_files": len(rows),
                "samples": sorted(sample_ids), "native_values_preserved": True, "upstream_jobs": 0}
            continue
        rendered = build_qc_assets(native, out / "rendered" / phase)
        phase_receipts[phase] = rendered["receipt"]
        for artifact in rendered["artifacts"]:
            if phase == "bam_records":
                metadata = artifact.setdefault("metadata", {})
                metadata["caption"] = (
                    metadata.get("caption", "").replace("FASTQ file", "supplied BAM file") +
                    " FastQC measures records retained in the supplied BAM; this is not original raw-read "
                    "or trimming QC and cannot recover reads absent from that BAM.")
            artifact.setdefault("metadata", {}).update(qc_phase=phase, original_native_QC=True,
                upstream_jobs=0, biological_sample_ids=sorted({row["biological_sample_id"] for row in rows}))
            artifacts.append(artifact)
    figures = [artifact for artifact in artifacts if Path(artifact["path"]).suffix == ".png"]
    if not figures:
        raise ValueError("Supplied sequencing QC contains no renderable measured QC data")
    # Gene-body coverage belongs to its separately scheduled measured analysis.
    # The existing parser's coverage placeholder is deliberately not promoted
    # into this independent sequencing-QC stage or its final report views.
    index = write_json(out / "sequencing-qc-views.json", {
        "schema": "katana.flow.native_sequencing_qc/v1", "status": "completed",
        "artifacts": artifacts, "sources": source_ledger, "phase_receipts": phase_receipts,
        "gene_body_coverage_scope": "separate_required_measured_workflow",
        "upstream_jobs": 0, "raw_reads_or_bams_read": False})
    paths = [path for path in out.rglob("*") if path.is_file()]
    outputs = [{"path": str(path), "size_bytes": path.stat().st_size, "sha256": file_sha256(path),
                "kind": "file", "name": path.name, "analysis_id": "flow.qc.sequencing"} for path in paths]
    evidence = {"qc_view_index": str(index), "native_qc_files": len(sources),
        "native_qc_figures": len(figures), "phases": sorted(phases), "upstream_jobs": 0,
        "raw_reads_or_bams_read": False, "native_sources_preserved": True,
        "expected_eligible_samples": sorted(sample.sample_id for sample in manifest.samples
            if sample.attributes.get("source", "").casefold() not in BAM_SOURCES),
        "gene_body_coverage_scope": "separate_required_measured_workflow"}
    return {"schema": "katana.flow.stage/v1", "adapter": "qc.sequencing", "status": "complete",
            "outputs": outputs, "evidence": evidence}
