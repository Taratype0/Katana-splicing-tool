"""Strict sequence/coordinate-bound external evidence, never causal inference."""
import hashlib
import json
import gzip

import numpy as np
import pandas as pd


def read_binding_peaks(path):
    """BED comments are whole lines; '#' inside an experimental peak ID is data."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        skip = [i for i, line in enumerate(stream) if line.startswith(("#", "track ", "browser ")) or not line.strip()]
    return pd.read_csv(path, sep="\t", header=None, skiprows=skip, compression="infer")


def species_name(value):
    return str(value).replace("_", " ").strip()


def verified_resource(asset, *, species, genome_build=None):
    metadata = asset.metadata
    required = {"source_url", "version", "sha256", "species"}
    if not required <= metadata.keys() or any(not str(metadata[k]).strip() for k in required):
        raise ValueError(f"{asset.asset_id}: source_url/version/sha256/species are required")
    if species_name(metadata["species"]) != species_name(species):
        raise ValueError(f"{asset.asset_id}: species mismatch; orthology is not binding evidence")
    if genome_build is not None and metadata.get("genome_build") != genome_build:
        raise ValueError(f"{asset.asset_id}: genome build mismatch; no implicit liftOver")
    with asset.path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    if digest != metadata["sha256"]:
        raise ValueError(f"{asset.asset_id}: SHA256 mismatch; resource changed")
    return {"asset_id": asset.asset_id, "path": str(asset.path), **metadata}


def protein_feature_support(features, resources, *, species):
    records, audits = [], []
    by_sequence = {}
    for asset in resources:
        audits.append(verified_resource(asset, species=species))
        entry = json.loads(asset.path.read_text(encoding="utf-8"))
        if species_name(entry.get("organism", {}).get("scientificName")) != species_name(species):
            raise ValueError("UniProt entry organism differs from the declared resource species")
        sequence = entry.get("sequence", {}).get("value", "")
        audit = entry.get("entryAudit", {})
        if not sequence or not audit.get("sequenceVersion") or not audit.get("entryVersion"):
            raise ValueError("UniProt JSON requires sequence, sequenceVersion and entryVersion")
        key = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        by_sequence.setdefault(key, []).append((entry, asset))
    for row in features.to_dict("records"):
        matches = by_sequence.get(row.get("protein_sha256"), []) if row.get("protein_status") == "validated_annotated_ORF" else []
        if not matches:
            records.append({"transcript_id": row["transcript_id"], "mapping_status": "not_evaluable_no_exact_sequence_match",
                            "evidence_type": "external_protein_annotation", "feature_type": "", "description": "",
                            "not_evaluable_reason": "no_exact_reference_sequence" if row.get("protein_status") == "validated_annotated_ORF" else "no_validated_complete_ORF"})
        for entry, asset in matches:
            found = False
            for feature in entry.get("features", []):
                if feature.get("type") not in {"Domain", "Region", "Motif", "Active site", "Binding site", "Signal", "Transmembrane"}:
                    continue
                location = feature.get("location", {})
                start, end = location.get("start", {}), location.get("end", {})
                if start.get("modifier", "EXACT") != "EXACT" or end.get("modifier", "EXACT") != "EXACT":
                    continue
                a, b = start.get("value"), end.get("value")
                if not isinstance(a, int) or not isinstance(b, int) or not 1 <= a <= b <= len(entry["sequence"]["value"]):
                    raise ValueError("Invalid exact UniProt protein feature coordinates")
                found = True
                records.append({"transcript_id": row["transcript_id"], "mapping_status": "exact_protein_sequence",
                    "evidence_type": "external_protein_annotation", "protein_accession": entry["primaryAccession"],
                    "protein_sequence_version": entry["entryAudit"]["sequenceVersion"], "entry_version": entry["entryAudit"]["entryVersion"],
                    "feature_type": feature["type"], "start_1based": a, "end_inclusive": b,
                    "description": feature.get("description", ""), "evidence_codes": json.dumps(feature.get("evidences", [])),
                    "source_url": asset.metadata["source_url"], "resource_sha256": asset.metadata["sha256"]})
            if not found:
                records.append({"transcript_id": row["transcript_id"], "mapping_status": "exact_sequence_no_supported_feature",
                                "evidence_type": "external_protein_annotation", "feature_type": "", "description": "",
                                "protein_accession": entry["primaryAccession"], "source_url": asset.metadata["source_url"],
                                "resource_sha256": asset.metadata["sha256"], "entry_version":entry["entryAudit"]["entryVersion"],
                                "protein_sequence_version":entry["entryAudit"]["sequenceVersion"]})
    return pd.DataFrame(records), pd.DataFrame(audits)


def binding_support(events, assets, *, species, genome_build, window_bp=100):
    from katana_cosplice.motif_scanner import _event_windows
    records, audits = [], []
    for asset in assets:
        audit = verified_resource(asset, species=species, genome_build=genome_build)
        if asset.metadata.get("coordinate_system") != "0-based-half-open":
            raise ValueError("Binding peaks must explicitly declare 0-based-half-open BED coordinates")
        for field in ("assay", "biosample", "target", "accession", "peak_calling"):
            if not asset.metadata.get(field):
                raise ValueError(f"Binding resource requires {field}; motifs cannot substitute for binding")
        if asset.metadata["assay"] not in {"eCLIP", "iCLIP", "PAR-CLIP", "HITS-CLIP"}:
            raise ValueError("Only explicitly declared experimental CLIP evidence is supported")
        audits.append(audit)
        peaks = read_binding_peaks(asset.path)
        if peaks.shape[1] < 6:
            raise ValueError("Binding peaks require BED6 or narrowPeak including strand")
        for col in (1, 2):
            peaks[col] = pd.to_numeric(peaks[col], errors="raise")
        if ((peaks[1] < 0) | (peaks[2] <= peaks[1]) | (peaks[1] != peaks[1].round()) | (peaks[2] != peaks[2].round())).any():
            raise ValueError("Invalid BED intervals")
        if not peaks[5].isin(["+", "-"]).all():
            raise ValueError("Unstranded peaks cannot supply strand-specific binding support")
        grouped = {(str(chrom), str(strand)): group for (chrom, strand), group in peaks.groupby([0, 5])}
        for row in events.to_dict("records"):
            from katana_cosplice.matrix import EVENT_COORDINATES
            cols = EVENT_COORDINATES.get(row.get("event_type"), [])
            if cols and all(pd.notna(row.get(c)) for c in cols):
                values = [float(row[c]) for c in cols]
                if all(np.isfinite(v) and v >= 0 and v == int(v) for v in values):
                    row["event_coordinates"] = "|".join(str(int(v)) for v in values)
            windows = _event_windows(pd.Series(row), window_bp)
            candidates = grouped.get((str(row.get("chr")), str(row.get("strand"))))
            supported = bool(windows) and row.get("strand") in {"+", "-"}
            hits = set()
            if supported and candidates is not None:
                for a, b in windows:
                    hits.update(candidates.index[(candidates[1] < b) & (candidates[2] > a)])
            records.append({"row_id": f"{row['row_id']}:{asset.asset_id}", "source_row_id": str(row["row_id"]),
                "resource_id": asset.asset_id,
                "event_id": row.get("event_id", row.get("event_key", "")), "gene_id": row.get("gene_id", ""),
                "module": row.get("module", ""), "comparison_id": row.get("comparison_id", ""),
                "RBP": asset.metadata["target"], "biosample": asset.metadata["biosample"],
                "assay": asset.metadata["assay"], "accession": asset.metadata["accession"],
                "overlapping_peaks": len(hits) if supported else None,
                "overlapping_peak_records": json.dumps([
                    {"bed_row_0based": int(i), "name": str(peaks.at[i, 3]), "chr": str(peaks.at[i, 0]),
                     "start_0based": int(peaks.at[i, 1]), "end_exclusive": int(peaks.at[i, 2]),
                     "strand": str(peaks.at[i, 5])} for i in sorted(hits)]),
                "queried_windows_0based_half_open": json.dumps(windows),
                "status": "overlap_in_reference_biosample" if hits else "no_overlap_in_resource" if supported else "not_evaluable",
                "not_evaluable_reason": "" if supported else (
                    "missing_or_unresolved_strand" if row.get("strand") not in {"+", "-"} else "missing_boundary_coordinates"),
                "evidence_type": "experimental_binding_in_reference_biosample_not_user_condition",
                "source_url": asset.metadata["source_url"], "resource_sha256": asset.metadata["sha256"]})
    return pd.DataFrame(records), pd.DataFrame(audits)
