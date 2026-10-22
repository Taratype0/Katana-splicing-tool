"""Read original SpliceSeq counts without redefining its native PSI."""
from dataclasses import dataclass
import csv
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class NativeSpliceSeqMeasurements:
    raw_psi: pd.DataFrame
    supported_psi: pd.DataFrame
    included_counts: pd.DataFrame
    excluded_counts: pd.DataFrame
    sample_audit: pd.DataFrame


def _number(value, *, integer=False):
    if value in {"", "\\N", "NULL", "NA"}:
        return np.nan
    result = float(value)
    if not np.isfinite(result) or (integer and (result < 0 or not result.is_integer())):
        raise ValueError("Invalid native count or numeric measurement: " + repr(value))
    return result


def read_native_spliceseq_counts(paths, event_map, samples, *, minimum_support=10):
    """Retain every declared event/sample; low support or absent data remain NA.

    Accept the unmodified `SELECT s.*,a.*` export, including its two numeric
    sample_id columns. Join biological samples through the explicit native
    sample name, never through the database-local numeric ID.
    """
    required = {"event_key", "native_as_id", "exact_full_path_identity"}
    if not required <= set(event_map) or not {"sample_id", "patient_id", "condition"} <= set(samples):
        raise ValueError("Explicit native event and biological sample mapping required")
    if minimum_support < 1 or not float(minimum_support).is_integer():
        raise ValueError("Read support must be a positive integer")
    if event_map.event_key.isna().any() or event_map.event_key.duplicated().any():
        raise ValueError("Frozen event keys must be nonmissing and unique")
    if samples[["sample_id", "patient_id", "condition"]].isna().any().any():
        raise ValueError("Missing biological identity")
    if samples.sample_id.duplicated().any() or samples.duplicated(["patient_id", "condition"]).any():
        raise ValueError("One declared library per patient and tissue scope required")
    ids = pd.to_numeric(event_map.native_as_id, errors="raise")
    if ids.dropna().duplicated().any() or (ids.dropna() % 1).ne(0).any():
        raise ValueError("Native event IDs must be unique integers")
    keys = event_map.event_key.astype(str).tolist()
    sample_names = samples.sample_id.astype(str).tolist()
    native_lookup = {int(v): i for i, v in enumerate(ids) if pd.notna(v)}
    sample_lookup = {s: i for i, s in enumerate(sample_names)}
    shape = (len(keys), len(sample_names))
    psi, included, excluded = [np.full(shape, np.nan) for _ in range(3)]
    seen = set()
    raw_rows = np.zeros(len(sample_names), dtype=int)
    for path in paths:
        with Path(path).open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            duplicate_headers = {h for h in header if header.count(h) > 1}
            if duplicate_headers - {"sample_id"}:
                raise ValueError("Ambiguous duplicate native columns")
            name = "name" if "name" in header else "sample_name"
            fields = [name, "as_id", "psi", "included_counts", "excluded_counts"]
            if not set(fields) <= set(header):
                raise ValueError("Missing original native export fields")
            positions = {f: header.index(f) for f in fields}
            duplicate_ids = [i for i, h in enumerate(header) if h == "sample_id"]
            for line, row in enumerate(reader, 2):
                if len(row) != len(header):
                    raise ValueError(f"Malformed native export row {line}")
                if len(duplicate_ids) > 1 and len({row[i] for i in duplicate_ids}) != 1:
                    raise ValueError("Native sample/count join numeric IDs disagree")
                sample = row[positions[name]]
                if sample not in sample_lookup:
                    raise ValueError("Undeclared native sample name: " + sample)
                j = sample_lookup[sample]
                raw_rows[j] += 1
                event_id = _number(row[positions["as_id"]], integer=True)
                if not np.isfinite(event_id):
                    raise ValueError("Missing native event ID")
                if int(event_id) not in native_lookup:
                    continue
                i = native_lookup[int(event_id)]
                if (i, j) in seen:
                    raise ValueError("Duplicate frozen event/sample measurement")
                seen.add((i, j))
                value = _number(row[positions["psi"]])
                if np.isfinite(value) and not 0 <= value <= 1:
                    raise ValueError("Native PSI outside [0, 1]")
                psi[i, j] = value
                included[i, j] = _number(row[positions["included_counts"]], integer=True)
                excluded[i, j] = _number(row[positions["excluded_counts"]], integer=True)
    support = included + excluded
    passed = np.isfinite(support) & (support >= minimum_support)
    supported = np.where(passed & np.isfinite(psi), psi, np.nan)
    frame = lambda a: pd.DataFrame(a, index=pd.Index(keys, name="event_key"), columns=sample_names)
    audit = samples[["sample_id", "patient_id", "condition"]].copy()
    audit["native_export_rows"] = raw_rows
    audit["frozen_events_with_native_psi"] = np.isfinite(psi).sum(axis=0)
    audit["frozen_events_supported"] = np.isfinite(supported).sum(axis=0)
    audit["low_or_missing_support_cells"] = (~passed).sum(axis=0)
    audit["library_observed_in_export"] = raw_rows > 0
    audit["minimum_included_plus_excluded_counts"] = minimum_support
    audit["psi_recomputed_from_counts"] = False
    return NativeSpliceSeqMeasurements(frame(psi), frame(supported), frame(included), frame(excluded), audit)
