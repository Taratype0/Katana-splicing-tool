from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.audit import load_sample_metadata, resolve_rmats_root
from katana_cosplice.config import CoSpliceConfig
from katana_cosplice.utils import normalize_gene_id, normalize_sample_id, parse_vector


EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
EVENT_COORDINATES: dict[str, list[str]] = {
    "SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "MXE": [
        "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
    ],
}


@dataclass
class MatrixResult:
    psi: pd.DataFrame
    ijc: pd.DataFrame
    sjc: pd.DataFrame
    events: pd.DataFrame
    samples: pd.DataFrame
    sources: pd.DataFrame


def _group_samples(path: Path) -> list[str]:
    raw = path.read_text(encoding="utf-8").replace("\r", "\n").replace("\n", ",")
    samples = [normalize_sample_id(value) for value in raw.split(",") if value.strip()]
    if len(samples) != len(set(samples)):
        duplicated = sorted({sample for sample in samples if samples.count(sample) > 1})
        raise ValueError(f"rMATS bamlist contains duplicate samples: {path}: {duplicated}")
    return samples


def _comparison_groups(name: str, known_groups: list[str]) -> tuple[str, str] | None:
    for left in known_groups:
        prefix = f"{left}-"
        if name.startswith(prefix):
            right = name[len(prefix):]
            if right in known_groups:
                return left, right
    return None


def _event_frame(frame: pd.DataFrame, event_type: str) -> pd.DataFrame:
    coordinate_columns = [column for column in EVENT_COORDINATES[event_type] if column in frame.columns]
    required = ["GeneID", "geneSymbol", "chr", "strand", *coordinate_columns]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{event_type} table is missing event identity columns: {missing}")
    if frame.loc[:, required].isna().any().any():
        bad = frame.loc[:, required].columns[frame.loc[:, required].isna().any()].tolist()
        raise ValueError(f"{event_type} table contains missing event identity values in columns: {bad}")
    output = frame.loc[:, required].copy()
    output["gene_id"] = output["GeneID"].map(normalize_gene_id)
    output["gene_symbol"] = output["geneSymbol"].astype(str).str.strip().str.strip('"')
    output["event_type"] = event_type
    identity = output[["gene_id", "event_type", "chr", "strand", *coordinate_columns]].astype(str)
    output["event_key"] = identity.agg("|".join, axis=1)
    output["event_coordinates"] = output[["chr", "strand", *coordinate_columns]].astype(str).agg("|".join, axis=1)
    keep = ["event_key", "gene_id", "gene_symbol", "event_type", "chr", "strand", "event_coordinates", *coordinate_columns]
    return output.loc[:, keep].drop_duplicates("event_key", keep="first")


def _vector_frame(
    frame: pd.DataFrame,
    event_keys: pd.Series,
    samples: list[str],
    vector_column: str,
) -> pd.DataFrame:
    values = frame[vector_column].map(parse_vector)
    wrong = values.map(len).ne(len(samples))
    if wrong.any():
        example = int(np.flatnonzero(wrong.to_numpy())[0])
        raise ValueError(
            f"{vector_column} vector length does not match bamlist sample count: "
            f"expected {len(samples)}, observed {len(values.iloc[example])}"
        )
    matrix = pd.DataFrame(values.tolist(), index=event_keys, columns=samples, dtype=float)
    matrix = matrix.loc[~matrix.index.duplicated(keep="first")]
    matrix.index.name = "event_key"
    return matrix


def _candidate_sources(
    rmats_root: Path,
    mode: str,
    groups: list[str],
) -> tuple[dict[tuple[str, str], tuple[Path, int]], dict[tuple[str, str], int]]:
    candidates: dict[tuple[str, str], list[tuple[Path, int]]] = {}
    for post_dir in rmats_root.glob("*/rmats_post"):
        pair = _comparison_groups(post_dir.parent.name, groups)
        if pair is None:
            continue
        for event_type in EVENT_TYPES:
            table = post_dir / f"{event_type}.MATS.{mode}.txt"
            if not table.exists():
                continue
            for side, group in enumerate(pair, start=1):
                key = (group, event_type)
                candidates.setdefault(key, []).append((table, side))
    sources = {
        key: sorted(values, key=lambda value: (-value[0].stat().st_size, str(value[0]), value[1]))[0]
        for key, values in candidates.items()
    }
    return sources, {key: len(values) for key, values in candidates.items()}


def _validate_matrix_values(psi: pd.DataFrame, ijc: pd.DataFrame, sjc: pd.DataFrame) -> None:
    if not psi.index.is_unique:
        raise ValueError("Sample-level PSI matrix contains duplicate event_key rows.")
    if not psi.columns.is_unique:
        raise ValueError("Sample-level PSI matrix contains duplicate sample columns.")
    finite_psi = psi.stack(future_stack=True).dropna()
    if not finite_psi.between(0.0, 1.0).all():
        raise ValueError("rMATS IncLevel values must be within [0, 1].")
    for name, matrix in (("IJC", ijc), ("SJC", sjc)):
        finite = matrix.stack(future_stack=True).dropna()
        if not finite.ge(0.0).all():
            raise ValueError(f"rMATS {name} read counts cannot be negative.")


def build_sample_matrices(config: CoSpliceConfig) -> MatrixResult:
    rmats_root = resolve_rmats_root(config.rmats_root)
    samples = load_sample_metadata(config.sample_metadata, config.selected_conditions)
    groups = sorted(samples["condition"].unique().tolist(), key=len, reverse=True)
    bamlist_dir = rmats_root / "bamlist"
    group_samples: dict[str, list[str]] = {}
    assigned_samples: set[str] = set()
    for group in groups:
        path = bamlist_dir / f"{group}_bamlist.txt"
        if not path.exists():
            raise FileNotFoundError(f"Missing rMATS bamlist for condition {group}: {path}")
        group_samples[group] = _group_samples(path)
        observed = set(group_samples[group])
        expected = set(samples.loc[samples["condition"].eq(group), "sample_id"])
        missing = sorted(observed.difference(samples["sample_id"]))
        if missing:
            raise ValueError(f"rMATS bamlist samples are absent from metadata for {group}: {missing}")
        if observed != expected:
            raise ValueError(
                f"rMATS bamlist does not match metadata condition {group}: "
                f"missing_from_bamlist={sorted(expected-observed)}, unexpected={sorted(observed-expected)}"
            )
        overlap = sorted(assigned_samples.intersection(observed))
        if overlap:
            raise ValueError(f"Samples occur in more than one rMATS condition bamlist: {overlap}")
        assigned_samples.update(observed)

    selections, candidate_counts = _candidate_sources(rmats_root, config.rmats_mode, groups)
    missing_sources = [(group, event) for group in groups for event in EVENT_TYPES if (group, event) not in selections]
    if missing_sources:
        raise FileNotFoundError(f"No rMATS source table found for: {missing_sources}")

    path_jobs: dict[Path, list[tuple[str, str, int]]] = {}
    for (group, event_type), (path, side) in selections.items():
        path_jobs.setdefault(path, []).append((group, event_type, side))

    psi_frames: dict[str, list[pd.DataFrame]] = {group: [] for group in groups}
    ijc_frames: dict[str, list[pd.DataFrame]] = {group: [] for group in groups}
    sjc_frames: dict[str, list[pd.DataFrame]] = {group: [] for group in groups}
    event_frames: list[pd.DataFrame] = []
    source_rows: list[dict[str, object]] = []

    for path, jobs in sorted(path_jobs.items(), key=lambda item: str(item[0])):
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        for group, event_type, side in jobs:
            metadata = _event_frame(frame, event_type)
            keys = metadata["event_key"]
            sample_ids = group_samples[group]
            psi_frames[group].append(_vector_frame(frame, keys, sample_ids, f"IncLevel{side}"))
            ijc_frames[group].append(_vector_frame(frame, keys, sample_ids, f"IJC_SAMPLE_{side}"))
            sjc_frames[group].append(_vector_frame(frame, keys, sample_ids, f"SJC_SAMPLE_{side}"))
            event_frames.append(metadata)
            source_rows.append(
                {
                    "condition": group,
                    "event_type": event_type,
                    "comparison": path.parent.parent.name,
                    "side": side,
                    "source_file": str(path),
                    "source_file_size": path.stat().st_size,
                    "candidate_source_count": candidate_counts[(group, event_type)],
                    "selection_rule": "largest_table_then_lexical_path",
                    "n_events": len(frame),
                    "samples": ",".join(sample_ids),
                }
            )

    sample_order = samples["sample_id"].tolist()
    psi_by_group = [pd.concat(psi_frames[group], axis=0) for group in groups]
    ijc_by_group = [pd.concat(ijc_frames[group], axis=0) for group in groups]
    sjc_by_group = [pd.concat(sjc_frames[group], axis=0) for group in groups]
    psi = pd.concat(psi_by_group, axis=1).reindex(columns=sample_order)
    ijc = pd.concat(ijc_by_group, axis=1).reindex(index=psi.index, columns=sample_order)
    sjc = pd.concat(sjc_by_group, axis=1).reindex(index=psi.index, columns=sample_order)
    _validate_matrix_values(psi, ijc, sjc)
    events = pd.concat(event_frames, ignore_index=True).drop_duplicates("event_key", keep="first")
    events = events.set_index("event_key").reindex(psi.index).reset_index()
    sources = pd.DataFrame(source_rows).sort_values(["condition", "event_type"]).reset_index(drop=True)
    return MatrixResult(psi=psi, ijc=ijc, sjc=sjc, events=events, samples=samples, sources=sources)
