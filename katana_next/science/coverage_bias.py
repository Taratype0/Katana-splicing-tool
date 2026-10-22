"""Read RSeQC gene-body profiles; descriptive bias is not a sample-failure rule."""
from pathlib import Path
import re

import numpy as np
import pandas as pd


def read_qualimap_summary(path: Path, *, sample_id):
    if not sample_id or not str(sample_id).strip():
        raise ValueError("Native Qualimap metrics require an explicit sample_id")
    text = path.read_text(encoding="utf-8")
    keys = {"qualimap_five_prime_bias": "5' bias", "qualimap_three_prime_bias": "3' bias",
            "qualimap_five_to_three_bias": "5'-3' bias", "bias_window_bp": "5'-3' bias region size",
            "bias_top_transcripts": "5'-3' bias number of top transcripts"}
    result = {"row_id": str(sample_id), "sample_id": str(sample_id)}
    for name, key in keys.items():
        match = re.search(r"^\s*" + re.escape(key) + r"\s*=\s*([\d.]+)\s*$", text, re.M)
        if not match:
            raise ValueError(f"Missing native Qualimap field {key}: {path}")
        value = float(match.group(1))
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"Invalid native Qualimap metric {key}")
        result[name] = value
    for name, key in (("qualimap_protocol", "protocol"), ("qualimap_gtf", "gff file")):
        match = re.search(r"^\s*" + re.escape(key) + r"\s*=\s*(.+)$", text, re.M)
        if not match:
            raise ValueError(f"Missing Qualimap input provenance {key}")
        result[name] = match.group(1).strip()
    result["metric_semantics"] = "Native median-derived terminal-window statistics; not ratios of aggregate percentile curves"
    return result


def read_gene_body_coverage(path: Path, sample_map: dict | None = None, *, sample_id=None) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t")
    if list(raw.columns) == ["#Transcript position", "Transcript coverage profile"]:
        if not sample_id:
            raise ValueError("Qualimap profiles require an explicit sample_id; filenames do not define sample identity")
        frame = raw.rename(columns={raw.columns[0]: "percentile", raw.columns[1]: "coverage"})
        frame.insert(0, "sample_id", sample_id)
        frame["axis_orientation"] = "source transcript position; use native Qualimap bias for strand-oriented end metrics"
    elif {"sample_id", "percentile", "coverage"} <= set(raw):
        frame = raw[["sample_id", "percentile", "coverage"]].copy()
    else:
        if len(raw.columns) < 11 or str(raw.columns[0]).casefold() not in {"percentile", "sample", "sample_id"}:
            raise ValueError("Expected RSeQC geneBodyCoverage.txt or sample_id/percentile/coverage TSV")
        frame = raw.rename(columns={raw.columns[0]: "sample_id"}).melt(
            id_vars="sample_id", var_name="percentile", value_name="coverage"
        )
    if "axis_orientation" not in frame:
        frame["axis_orientation"] = "5-prime to 3-prime"
    if frame.empty or frame.sample_id.isna().any() or frame.sample_id.astype(str).str.strip().eq("").any():
        raise ValueError("Coverage rows need explicit sample IDs")
    frame["sample_id"] = frame.sample_id.astype(str)
    if sample_map:
        frame["sample_id"] = frame.sample_id.map(lambda sample: sample_map.get(sample, sample))
    if frame.sample_id.isna().any() or frame.sample_id.astype(str).str.strip().eq("").any():
        raise ValueError("Coverage sample mapping needs nonempty sample IDs")
    for column in ("percentile", "coverage"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if not np.isfinite(frame[["percentile", "coverage"]].to_numpy()).all():
        raise ValueError("Nonfinite gene-body coverage")
    if (frame.coverage.to_numpy() < 0).any() or ((frame.percentile.to_numpy() < 0) | (frame.percentile.to_numpy() > 100)).any():
        raise ValueError("Coverage must be nonnegative; percentiles must be between 0 and 100")
    if frame.duplicated(["sample_id", "percentile"]).any():
        raise ValueError("Duplicate sample/percentile values; do not merge profiles implicitly")
    for sample, group in frame.groupby("sample_id"):
        if len(group) < 10 or group.percentile.min() > 5 or group.percentile.max() < 95:
            raise ValueError(f"Incomplete 5-prime to 3-prime profile: {sample}")
    return frame.sort_values(["sample_id", "percentile"]).reset_index(drop=True)


def summarize_coverage(profiles: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if metadata.sample_id.duplicated().any():
        raise ValueError("Coverage metadata needs one row per sample")
    unknown = set(profiles.sample_id) - set(metadata.sample_id)
    if unknown:
        raise ValueError(f"Coverage samples missing from metadata: {sorted(unknown)[:10]}; declare sample_map if needed")
    curves = profiles.merge(metadata, on="sample_id", validate="many_to_one")
    means = curves.groupby("sample_id").coverage.transform("mean")
    curves["relative_coverage"] = curves.coverage.div(means.replace(0, np.nan))
    rows = []
    for sample, group in curves.groupby("sample_id", sort=False):
        five = group.loc[group.percentile.le(20), "coverage"].mean()
        three = group.loc[group.percentile.ge(80), "coverage"].mean()
        oriented = "axis_orientation" not in group or group.axis_orientation.eq("5-prime to 3-prime").all()
        ratio = three / five if five > 0 else np.nan
        rows.append({
            "row_id": sample, "sample_id": sample,
            "mean_coverage": group.coverage.mean(), "five_prime_mean": five, "three_prime_mean": three,
            "three_to_five_ratio": ratio if oriented else np.nan,
            "position_end_to_start_ratio": ratio,
            "end_ratio_semantics": "last/first 20-percent windows; not native Qualimap median 100-bp bias",
            "profile_status": "measurable" if group.coverage.sum() > 0 else "no_coverage",
            "interpretation": "Descriptive only; no automatic pass/fail or correction",
        })
    curves["row_id"] = curves.sample_id + "::" + curves.percentile.astype(str)
    metrics = pd.DataFrame(rows).merge(metadata, on="sample_id", validate="one_to_one")
    return curves, metrics
