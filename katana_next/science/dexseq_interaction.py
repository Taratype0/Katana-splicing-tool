"""Explicit two-factor exon-usage interaction; never infer design from filenames."""
from pathlib import Path
import pandas as pd


def read_design(path, reference_a, reference_b, sample_ids=None):
    path = Path(path)
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"sample_id", "factor_a", "factor_b", "count_file"}
    if not required <= set(frame):
        raise ValueError(f"Design requires columns: {sorted(required)}")
    if frame[list(required)].eq("").any().any() or frame.sample_id.duplicated().any():
        raise ValueError("Design requires nonempty values and unique sample IDs")
    if "biological_sample_id" in frame and frame.biological_sample_id.duplicated().any():
        raise ValueError("Repeated biological units require a paired model; this workflow is independent-sample only")
    if sample_ids is not None and not set(frame.sample_id) <= set(sample_ids):
        raise ValueError("Design contains samples outside the project manifest")
    for column, reference in [("factor_a", reference_a), ("factor_b", reference_b)]:
        if frame[column].nunique() != 2 or reference not in set(frame[column]):
            raise ValueError(f"{column} requires exactly two levels and an explicit valid reference")
    groups = frame.groupby(["factor_a", "factor_b"]).size()
    if len(groups) != 4 or groups.min() < 2:
        raise ValueError("A complete 2x2 design with at least two independent samples per cell is required")
    paths = []
    for value in frame.count_file:
        candidate = Path(value)
        candidate = candidate if candidate.is_absolute() else path.parent / candidate
        candidate = candidate.resolve()
        if not candidate.is_file():
            raise FileNotFoundError(candidate)
        paths.append(str(candidate))
    if len(set(paths)) != len(paths):
        raise ValueError("A count file cannot represent multiple independent samples")
    frame["count_file"] = paths
    return frame
