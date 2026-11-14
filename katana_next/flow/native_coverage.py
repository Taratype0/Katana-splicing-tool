"""Select completed coverage profiles without opening a BAM or running QC."""
from pathlib import Path
import json

from katana_next.science.coverage_bias import read_gene_body_coverage, validate_coverage_frame


def multiqc_qualimap_profiles(asset, declared):
    """Read Qualimap's saved raw curve, not its rescaled plot or general statistics."""
    import pandas as pd
    document = json.loads(asset.path.read_text(encoding="utf-8-sig"))
    raw = document.get("report_saved_raw_data", {}).get("qualimap_rnaseq_cov_hist", {})
    if not raw:
        return None
    rows = []
    sample_map = asset.metadata.get("sample_map", {})
    for original, values in raw.items():
        sid = sample_map.get(original, original)
        # Coverage is one profile per biological sample, not per read mate.
        if sid not in declared:
            raise ValueError("MultiQC Qualimap curve needs exact sample identity or sample_map: " + original)
        if not isinstance(values, dict):
            raise ValueError("Unsupported MultiQC Qualimap raw curve: " + original)
        rows.extend({"sample_id": sid, "percentile": position, "coverage": value,
                     "axis_orientation": "source transcript position; use native Qualimap bias for strand-oriented end metrics"}
                    for position, value in values.items())
    return validate_coverage_frame(pd.DataFrame(rows))


def profile_selection(manifest):
    """One source per declared sample, preferring Qualimap over existing RSeQC."""
    declared = {s.sample_id for s in manifest.samples}
    candidates = {}
    frames = {}
    for asset in manifest.assets:
        if asset.kind == "sequencing_qc" and asset.path.name == "multiqc_data.json":
            frame = multiqc_qualimap_profiles(asset, declared)
            if frame is None:
                continue
            priority = 1
        elif asset.kind == "gene_body_coverage":
            frame = read_gene_body_coverage(asset.path, asset.metadata.get("sample_map"),
                                            sample_id=asset.metadata.get("sample_id"))
            qualimap = frame.axis_orientation.str.startswith("source transcript position").all()
            priority = 0 if qualimap else 2
        else:
            continue
        unknown = set(frame.sample_id) - declared
        if unknown:
            raise ValueError("Coverage profile contains undeclared samples: " + ", ".join(sorted(unknown)))
        frames[asset.asset_id] = frame
        for sid in frame.sample_id.unique():
            candidates.setdefault(sid, []).append((priority, asset.asset_id))
    selected = {}
    for sid, choices in candidates.items():
        priority = min(x[0] for x in choices)
        best = sorted(aid for rank, aid in choices if rank == priority)
        if len(best) != 1:
            raise ValueError(f"Multiple coverage profiles of the same priority for {sid}: {best}")
        selected[sid] = best[0]
    return selected, frames


def bind_native_sample_metadata(asset, sample_ids):
    """Resolve explicitly supplied per-sample QC directories against the sample sheet."""
    metadata = asset.setdefault("metadata", {})
    if asset.get("kind") not in {"gene_body_coverage", "qualimap_rnaseq_metrics"}:
        return
    if metadata.get("sample_id"):
        if metadata["sample_id"] not in sample_ids:
            raise ValueError("Unknown coverage sample: " + str(metadata["sample_id"]))
        return
    relative = metadata.get("qc_relative_path", "")
    candidates = set(Path(relative).parts[:-1]) & set(sample_ids)
    if len(candidates) != 1:
        raise ValueError("Put native QC in a directory named exactly as its sample ID, or declare sample_id: " + relative)
    metadata["sample_id"] = candidates.pop()


def stage_native_qc(source, destination):
    """Copy only small completed QC products, never recursively copy reads/BAMs."""
    import shutil
    source, destination = Path(source), Path(destination)
    names = {"coverage_profile_along_genes_(total).txt": "gene_body_coverage",
             "rnaseq_qc_results.txt": "qualimap_rnaseq_metrics",
             "multiqc_report.html": "sequencing_qc", "multiqc_data.json": "sequencing_qc"}
    assets = []
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        kind = names.get(path.name)
        if kind is None:
            continue
        relative = path.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        metadata = {"qc_relative_path": relative.as_posix(), "source_path": str(path.resolve()),
                    "reused_native_result": True}
        if path.name.startswith("multiqc_"):
            metadata.update(qc_type="multiqc", qc_phase="multiqc")
        assets.append({"asset_id": "native-qc-" + str(len(assets)), "kind": kind,
                       "path": target.relative_to(destination.parent).as_posix(), "metadata": metadata})
    if not assets:
        raise ValueError("--qc_results contains no supported Qualimap or MultiQC results")
    return assets
