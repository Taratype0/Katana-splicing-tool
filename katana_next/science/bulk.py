from __future__ import annotations

import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from katana_next.science.statistics import (
    bh_adjust,
    hypergeometric_overlap,
    safe_spearman,
)
from katana_next.infrastructure.fingerprints import fingerprint_value
from katana_next.services.project_service import ProjectContext


EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
EVENT_COORDINATES = (
    "chr",
    "strand",
    "exonStart_0base",
    "exonEnd",
    "upstreamES",
    "upstreamEE",
    "downstreamES",
    "downstreamEE",
    "riExonStart_0base",
    "riExonEnd",
    "longExonStart_0base",
    "longExonEnd",
    "shortES",
    "shortEE",
    "flankingES",
    "flankingEE",
    "1stExonStart_0base",
    "1stExonEnd",
    "2ndExonStart_0base",
    "2ndExonEnd",
)


@dataclass(frozen=True)
class BulkEventData:
    events: pd.DataFrame
    psi_matrix: pd.DataFrame
    read_matrix: pd.DataFrame
    matrix_audit: pd.DataFrame
    source_cache_status: str = "rebuild"


def clean_identifier(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip().strip('"').strip("'")


def normalize_gene_id(value: object) -> str:
    text = clean_identifier(value)
    match = re.match(r"^(ENS[A-Z]*G\d+)(?:\.\d+)?$", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return text


def optional_series(
    frame: pd.DataFrame, column: str, default: object = ""
) -> pd.Series:
    """Return a Series for optional input columns instead of a scalar default."""
    if column in frame.columns:
        return frame[column]
    return pd.Series(default, index=frame.index, dtype="object")


def parse_number_list(value: object) -> list[float]:
    if value is None:
        return []
    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass
    values: list[float] = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if not token or token.upper() in {"NA", "NAN", "NONE"}:
            continue
        try:
            number = float(token)
        except ValueError:
            continue
        if math.isfinite(number):
            values.append(number)
    return values


def canonical_event_id(row: pd.Series, event_type: str) -> str:
    values = [event_type]
    for column in EVENT_COORDINATES:
        if column not in row.index:
            continue
        value = clean_identifier(row[column])
        if value:
            values.append(f"{column}={value}")
    if len(values) == 1:
        fallback = clean_identifier(row.get("ID", row.name))
        values.append(f"source_id={fallback}")
    return "|".join(values)


def canonical_event_ids(frame: pd.DataFrame, event_type: str) -> pd.Series:
    """Build the same exact IDs as ``canonical_event_id`` without row apply."""
    identifiers = pd.Series(event_type, index=frame.index, dtype="object")
    coordinate_count = pd.Series(0, index=frame.index, dtype="int16")
    for column in EVENT_COORDINATES:
        if column not in frame.columns:
            continue
        values = frame[column].map(clean_identifier)
        present = values.ne("")
        identifiers.loc[present] = (
            identifiers.loc[present]
            + "|"
            + column
            + "="
            + values.loc[present]
        )
        coordinate_count.loc[present] += 1
    fallback = coordinate_count.eq(0)
    if fallback.any():
        source = optional_series(frame, "ID", "").map(clean_identifier)
        source.loc[source.eq("")] = source.index[source.eq("")].astype(str)
        identifiers.loc[fallback] = (
            identifiers.loc[fallback] + "|source_id=" + source.loc[fallback]
        )
    return identifiers


def selected_contrast_ids(
    context: ProjectContext, parameters: dict[str, Any]
) -> list[str]:
    available = [contrast.contrast_id for contrast in context.manifest.contrasts]
    requested = parameters.get("contrast_ids") or available
    if isinstance(requested, str):
        requested = [item.strip() for item in requested.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(available))
    if unknown:
        raise ValueError(f"Unknown contrast IDs: {', '.join(unknown)}")
    return [contrast_id for contrast_id in available if contrast_id in requested]


def rmats_assets(context: ProjectContext, contrast_ids: list[str]) -> list[Any]:
    assets = []
    for asset in context.manifest.assets:
        if asset.kind != "rmats_result" or not asset.path.is_file():
            continue
        contrast_id = str(asset.metadata.get("contrast_id", ""))
        if contrast_id in contrast_ids:
            assets.append(asset)
    return assets


def _contrast(context: ProjectContext, contrast_id: str):
    match = next(
        (item for item in context.manifest.contrasts if item.contrast_id == contrast_id),
        None,
    )
    if match is None:
        raise ValueError(f"Unknown contrast: {contrast_id}")
    return match


def _sample_names(asset: Any, contrast: Any, arm: str, count: int) -> list[str]:
    key = f"{arm}_sample_ids"
    names = asset.metadata.get(key)
    if not isinstance(names, (list, tuple)):
        names = (
            contrast.numerator.sample_ids
            if arm == "numerator"
            else contrast.denominator.sample_ids
        )
    normalized = [str(value) for value in names]
    if len(normalized) < count:
        normalized.extend(
            f"{contrast.contrast_id}:{arm}:replicate_{index + 1}"
            for index in range(len(normalized), count)
        )
    return normalized[:count]


def load_bulk_events(
    context: ProjectContext,
    parameters: dict[str, Any],
) -> BulkEventData:
    contrast_ids = selected_contrast_ids(context, parameters)
    assets = rmats_assets(context, contrast_ids)
    cache_directory = _bulk_event_cache_directory(context, contrast_ids, assets)
    cached = _read_bulk_event_cache(cache_directory)
    if cached is not None:
        return cached
    frames: list[pd.DataFrame] = []
    psi_observations: dict[tuple[str, str], list[float]] = {}
    read_observations: dict[tuple[str, str], list[float]] = {}
    source_observations: dict[tuple[str, str], set[str]] = {}

    for asset in assets:
        contrast_id = str(asset.metadata["contrast_id"])
        event_type = str(asset.metadata.get("event_type") or _event_type_from_name(asset.path.name))
        if event_type not in EVENT_TYPES:
            raise ValueError(f"Unknown rMATS event type for {asset.path}: {event_type}")
        source = pd.read_csv(asset.path, sep="\t", low_memory=False)
        required = {"FDR", "IncLevelDifference", "IncLevel1", "IncLevel2"}
        missing = sorted(required - set(source.columns))
        if missing:
            raise ValueError(f"rMATS file {asset.path.name} is missing columns: {missing}")
        if source.empty:
            continue
        contrast = _contrast(context, contrast_id)
        frame = source.copy()
        frame["event_type"] = event_type
        frame["event_id"] = canonical_event_ids(frame, event_type)
        frame["comparison_id"] = contrast_id
        frame["comparison_label"] = contrast.label
        frame["numerator"] = contrast.numerator.label
        frame["denominator"] = contrast.denominator.label
        frame["gene_id"] = optional_series(frame, "GeneID").map(normalize_gene_id)
        frame["gene_symbol"] = optional_series(frame, "geneSymbol").map(
            clean_identifier
        )
        frame["fdr"] = pd.to_numeric(frame["FDR"], errors="coerce")
        frame["dpsi"] = pd.to_numeric(frame["IncLevelDifference"], errors="coerce")
        reverse_direction = bool(asset.metadata.get("reverse_direction", False))
        frame["reverse_direction"] = reverse_direction
        if reverse_direction:
            frame["dpsi"] = -frame["dpsi"]
            arm1 = "denominator"
        else:
            arm1 = "numerator"
        frame["abs_dpsi"] = frame["dpsi"].abs()
        frame["row_id"] = contrast_id + "::" + frame["event_id"]
        frame["source_path"] = str(asset.path)
        frame["source_row"] = np.arange(len(frame), dtype=int)

        parsed = {
            column: optional_series(frame, column).map(parse_number_list)
            for column in (
                "IncLevel1",
                "IncLevel2",
                "IJC_SAMPLE_1",
                "IJC_SAMPLE_2",
                "SJC_SAMPLE_1",
                "SJC_SAMPLE_2",
            )
        }
        metric_rows = []
        name_cache: dict[tuple[str, int], list[str]] = {}
        rows = zip(
            frame["event_id"].astype(str),
            frame["dpsi"],
            parsed["IncLevel1"],
            parsed["IncLevel2"],
            parsed["IJC_SAMPLE_1"],
            parsed["IJC_SAMPLE_2"],
            parsed["SJC_SAMPLE_1"],
            parsed["SJC_SAMPLE_2"],
            strict=True,
        )
        for event_id, dpsi, psi1, psi2, ijc1, ijc2, sjc1, sjc2 in rows:
            if arm1 == "numerator":
                numerator_psi, denominator_psi = psi1, psi2
                numerator_ijc, denominator_ijc = ijc1, ijc2
                numerator_sjc, denominator_sjc = sjc1, sjc2
            else:
                numerator_psi, denominator_psi = psi2, psi1
                numerator_ijc, denominator_ijc = ijc2, ijc1
                numerator_sjc, denominator_sjc = sjc2, sjc1
            numerator_reads = _paired_sums(numerator_ijc, numerator_sjc)
            denominator_reads = _paired_sums(denominator_ijc, denominator_sjc)
            numerator_key = ("numerator", len(numerator_psi))
            denominator_key = ("denominator", len(denominator_psi))
            if numerator_key not in name_cache:
                name_cache[numerator_key] = _sample_names(
                    asset, contrast, *numerator_key
                )
            if denominator_key not in name_cache:
                name_cache[denominator_key] = _sample_names(
                    asset, contrast, *denominator_key
                )
            numerator_names = name_cache[numerator_key]
            denominator_names = name_cache[denominator_key]
            for sample, psi, reads in zip(
                numerator_names,
                numerator_psi,
                _pad(numerator_reads, len(numerator_psi)),
                strict=True,
            ):
                psi_observations.setdefault((event_id, sample), []).append(psi)
                read_observations.setdefault((event_id, sample), []).append(reads)
                source_observations.setdefault((event_id, sample), set()).add(contrast_id)
            for sample, psi, reads in zip(
                denominator_names,
                denominator_psi,
                _pad(denominator_reads, len(denominator_psi)),
                strict=True,
            ):
                psi_observations.setdefault((event_id, sample), []).append(psi)
                read_observations.setdefault((event_id, sample), []).append(reads)
                source_observations.setdefault((event_id, sample), set()).add(contrast_id)
            metric_rows.append(
                _replicate_metrics(
                    numerator_psi,
                    denominator_psi,
                    numerator_reads,
                    denominator_reads,
                    float(dpsi) if pd.notna(dpsi) else math.nan,
                )
            )
        metrics = pd.DataFrame(metric_rows, index=frame.index)
        frame = pd.concat([frame, metrics], axis=1)
        frames.append(frame)

    if not frames:
        raise ValueError("No compatible rMATS result assets are available for the selected contrasts.")
    events = pd.concat(frames, ignore_index=True, sort=False)
    psi_matrix = _observation_matrix(psi_observations)
    read_matrix = _observation_matrix(read_observations)
    audit_rows = []
    for key, values in psi_observations.items():
        event_id, sample_id = key
        audit_rows.append(
            {
                "row_id": f"{event_id}::{sample_id}",
                "event_id": event_id,
                "sample_id": sample_id,
                "source_comparisons": ", ".join(sorted(source_observations[key])),
                "n_observations": len(values),
                "median_psi": float(np.median(values)),
                "max_absolute_disagreement": float(np.max(values) - np.min(values)),
                "matrix_value_rule": "median across duplicate comparison observations",
            }
        )
    audit = pd.DataFrame(audit_rows)
    result = BulkEventData(
        events=events,
        psi_matrix=psi_matrix,
        read_matrix=read_matrix,
        matrix_audit=audit,
        source_cache_status="rebuild",
    )
    _write_bulk_event_cache(cache_directory, result, contrast_ids)
    return result


def _bulk_event_cache_directory(
    context: ProjectContext, contrast_ids: list[str], assets: list[Any]
) -> Path:
    fingerprints = {
        str(item.path.resolve()): item.model_dump(mode="json")
        for item in context.snapshot.input_fingerprints
    }
    payload = {
        "loader_version": "2.1",
        "contrast_ids": contrast_ids,
        "contrasts": [
            contrast.model_dump(mode="json")
            for contrast in context.manifest.contrasts
            if contrast.contrast_id in contrast_ids
        ],
        "assets": [
            {
                "path": str(asset.path.resolve()),
                "metadata": asset.metadata,
                "fingerprint": fingerprints.get(str(asset.path.resolve())),
            }
            for asset in assets
        ],
    }
    return (
        context.manifest.output_root
        / ".katana"
        / "cache"
        / "bulk_events"
        / fingerprint_value(payload)[:20]
    )


def _read_bulk_event_cache(directory: Path) -> BulkEventData | None:
    paths = {
        "events": directory / "events.parquet",
        "psi": directory / "psi_matrix.parquet",
        "reads": directory / "read_matrix.parquet",
        "audit": directory / "matrix_audit.parquet",
        "metadata": directory / "cache.json",
    }
    if not all(path.is_file() for path in paths.values()):
        return None
    try:
        psi = pd.read_parquet(paths["psi"]).set_index("event_id")
        reads = pd.read_parquet(paths["reads"]).set_index("event_id")
        return BulkEventData(
            events=pd.read_parquet(paths["events"]),
            psi_matrix=psi,
            read_matrix=reads,
            matrix_audit=pd.read_parquet(paths["audit"]),
            source_cache_status="hit",
        )
    except (OSError, ValueError, KeyError):
        return None


def _write_bulk_event_cache(
    directory: Path, result: BulkEventData, contrast_ids: list[str]
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    _atomic_parquet(result.events, directory / "events.parquet")
    _atomic_parquet(
        result.psi_matrix.rename_axis("event_id").reset_index(),
        directory / "psi_matrix.parquet",
    )
    _atomic_parquet(
        result.read_matrix.rename_axis("event_id").reset_index(),
        directory / "read_matrix.parquet",
    )
    _atomic_parquet(result.matrix_audit, directory / "matrix_audit.parquet")
    metadata = directory / "cache.json"
    temporary = metadata.with_name(f"{metadata.name}.{uuid4().hex}.tmp")
    temporary.write_text(
        json.dumps(
            {
                "schema_version": "2.1",
                "contrast_ids": contrast_ids,
                "semantics": "parameter-independent normalized rMATS event and replicate cache",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    os.replace(temporary, metadata)


def _atomic_parquet(frame: pd.DataFrame, target: Path) -> None:
    temporary = target.with_name(f"{target.name}.{uuid4().hex}.tmp")
    frame.to_parquet(temporary, index=False)
    os.replace(temporary, target)


def _event_type_from_name(name: str) -> str:
    upper = name.upper()
    return next((event_type for event_type in EVENT_TYPES if upper.startswith(event_type + ".")), "")


def _paired_sums(left: list[float], right: list[float]) -> list[float]:
    length = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else 0.0)
        + (right[index] if index < len(right) else 0.0)
        for index in range(length)
    ]


def _pad(values: list[float], length: int) -> list[float]:
    return [values[index] if index < len(values) else math.nan for index in range(length)]


def _replicate_metrics(
    numerator_psi: list[float],
    denominator_psi: list[float],
    numerator_reads: list[float],
    denominator_reads: list[float],
    dpsi: float,
) -> dict[str, float | int]:
    differences = [
        left - right
        for left in numerator_psi
        for right in denominator_psi
        if math.isfinite(left) and math.isfinite(right) and left != right
    ]
    expected = dpsi > 0
    pairwise = (
        float(np.mean([(difference > 0) == expected for difference in differences]))
        if differences and math.isfinite(dpsi) and dpsi != 0
        else math.nan
    )
    loo_effects = []
    for group, other, sign in (
        (numerator_psi, denominator_psi, 1.0),
        (denominator_psi, numerator_psi, -1.0),
    ):
        for index in range(len(group)):
            reduced = [value for position, value in enumerate(group) if position != index]
            if reduced and other:
                loo_effects.append(sign * (float(np.mean(reduced)) - float(np.mean(other))))
    loo_consistency = (
        float(np.mean([(effect > 0) == expected for effect in loo_effects if effect != 0]))
        if loo_effects and math.isfinite(dpsi) and dpsi != 0
        else math.nan
    )
    uncertainty = []
    for psi, reads in [
        *zip(numerator_psi, _pad(numerator_reads, len(numerator_psi)), strict=True),
        *zip(denominator_psi, _pad(denominator_reads, len(denominator_psi)), strict=True),
    ]:
        if math.isfinite(psi) and math.isfinite(reads) and reads > 0:
            uncertainty.append(math.sqrt(max(psi * (1 - psi), 0.0) / (reads + 1.0)))
    return {
        "numerator_replicates": len(numerator_psi),
        "denominator_replicates": len(denominator_psi),
        "numerator_mean_reads": float(np.mean(numerator_reads)) if numerator_reads else math.nan,
        "denominator_mean_reads": float(np.mean(denominator_reads)) if denominator_reads else math.nan,
        "replicate_pairwise_direction_consistency": pairwise,
        "leave_one_out_direction_consistency": loo_consistency,
        "max_measurement_se": max(uncertainty) if uncertainty else math.nan,
    }


def _observation_matrix(values: dict[tuple[str, str], list[float]]) -> pd.DataFrame:
    if not values:
        return pd.DataFrame()
    event_ids = sorted({event for event, _ in values})
    sample_ids = sorted({sample for _, sample in values})
    matrix = pd.DataFrame(index=event_ids, columns=sample_ids, dtype=float)
    for (event_id, sample_id), observations in values.items():
        finite = [value for value in observations if math.isfinite(value)]
        matrix.loc[event_id, sample_id] = float(np.median(finite)) if finite else math.nan
    matrix.index.name = "event_id"
    return matrix


STABLE_PRESETS: dict[str, dict[str, float]] = {
    "exploratory": {
        "min_reads": 8.0,
        "min_supported_replicates": 2.0,
        "min_direction_consistency": 0.67,
        "max_measurement_se": 0.18,
    },
    "balanced": {
        "min_reads": 15.0,
        "min_supported_replicates": 2.0,
        "min_direction_consistency": 0.75,
        "max_measurement_se": 0.12,
    },
    "stringent": {
        "min_reads": 20.0,
        "min_supported_replicates": 3.0,
        "min_direction_consistency": 0.80,
        "max_measurement_se": 0.10,
    },
}


def apply_stable_event_qc(events: pd.DataFrame, parameters: dict[str, Any]) -> pd.DataFrame:
    profile = str(parameters.get("stability_profile", "balanced"))
    if profile not in {*STABLE_PRESETS, "custom"}:
        raise ValueError(f"Unknown stability profile: {profile}")
    criteria = dict(STABLE_PRESETS.get(profile, STABLE_PRESETS["balanced"]))
    if profile == "custom":
        for key in criteria:
            if key in parameters:
                criteria[key] = float(parameters[key])
    frame = events.copy()
    min_reads = criteria["min_reads"]
    min_replicates = int(criteria["min_supported_replicates"])
    for arm in ("numerator", "denominator"):
        read_lists = frame.apply(
            lambda row: _arm_read_totals(row, arm), axis=1
        )
        frame[f"{arm}_supported_replicates"] = read_lists.map(
            lambda values: sum(value >= min_reads for value in values)
        )
    checks = {
        "replicate_count": frame["numerator_replicates"].ge(min_replicates)
        & frame["denominator_replicates"].ge(min_replicates),
        "read_support": frame["numerator_supported_replicates"].ge(min_replicates)
        & frame["denominator_supported_replicates"].ge(min_replicates),
        "direction_consistency": pd.to_numeric(
            frame["replicate_pairwise_direction_consistency"], errors="coerce"
        ).ge(criteria["min_direction_consistency"]),
        "leave_one_out": pd.to_numeric(
            frame["leave_one_out_direction_consistency"], errors="coerce"
        ).ge(criteria["min_direction_consistency"]),
        "measurement_uncertainty": pd.to_numeric(
            frame["max_measurement_se"], errors="coerce"
        ).le(criteria["max_measurement_se"]),
    }
    for name, values in checks.items():
        frame[f"passes_{name}"] = values.fillna(False)
    frame["stable_event"] = np.logical_and.reduce(
        [frame[f"passes_{name}"].to_numpy(dtype=bool) for name in checks]
    )
    frame["stability_profile"] = profile
    frame["stability_failure_reason"] = frame.apply(
        lambda row: ", ".join(
            name for name in checks if not bool(row[f"passes_{name}"])
        )
        or "passed",
        axis=1,
    )
    for key, value in criteria.items():
        frame[f"criterion_{key}"] = value
    return frame


def _arm_read_totals(row: pd.Series, arm: str) -> list[float]:
    reverse = bool(row.get("reverse_direction", False))
    if arm == "numerator":
        suffix = "2" if reverse else "1"
    else:
        suffix = "1" if reverse else "2"
    return _paired_sums(
        parse_number_list(row.get(f"IJC_SAMPLE_{suffix}")),
        parse_number_list(row.get(f"SJC_SAMPLE_{suffix}")),
    )


def significance_columns(events: pd.DataFrame, parameters: dict[str, Any]) -> pd.DataFrame:
    frame = events.copy()
    fdr_cutoff = float(parameters.get("fdr", 0.05))
    dpsi_cutoff = float(parameters.get("min_abs_dpsi", 0.10))
    frame["significant"] = frame["fdr"].le(fdr_cutoff) & frame["abs_dpsi"].ge(dpsi_cutoff)
    frame["stable_significant"] = frame["significant"] & frame.get(
        "stable_event", pd.Series(False, index=frame.index)
    ).fillna(False)
    frame["direction"] = np.select(
        [frame["dpsi"].gt(0), frame["dpsi"].lt(0)],
        ["higher inclusion in numerator", "lower inclusion in numerator"],
        default="no change",
    )
    frame["fdr_cutoff"] = fdr_cutoff
    frame["min_abs_dpsi"] = dpsi_cutoff
    return frame


def landscape_tables(
    data: BulkEventData, *, max_pca_events: int = 2000
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    events = data.events
    summary = (
        events.groupby(["comparison_id", "comparison_label", "event_type"], dropna=False)
        .agg(
            tested_events=("event_id", "size"),
            unique_events=("event_id", "nunique"),
            median_abs_dpsi=("abs_dpsi", "median"),
            median_fdr=("fdr", "median"),
        )
        .reset_index()
    )
    pca = sample_pca(data.psi_matrix, max_events=max_pca_events)
    similarity = sample_similarity(data.psi_matrix)
    return summary, pca, similarity


def sample_pca(matrix: pd.DataFrame, max_events: int = 2000) -> pd.DataFrame:
    if matrix.shape[1] < 3 or matrix.shape[0] < 2:
        return pd.DataFrame(columns=["row_id", "sample_id", "PC1", "PC2", "explained_PC1", "explained_PC2"])
    working = matrix.copy()
    variance = working.var(axis=1, skipna=True).sort_values(ascending=False)
    selected = variance[variance.gt(0)].head(max_events).index
    working = working.loc[selected]
    if len(working) < 2:
        return pd.DataFrame(columns=["row_id", "sample_id", "PC1", "PC2", "explained_PC1", "explained_PC2"])
    working = working.apply(lambda row: row.fillna(row.median()), axis=1).fillna(0.5)
    values = working.T.to_numpy(dtype=float)
    values = (values - np.mean(values, axis=0)) / np.where(np.std(values, axis=0) > 0, np.std(values, axis=0), 1.0)
    n_components = min(2, values.shape[0], values.shape[1])
    model = PCA(n_components=n_components, random_state=17)
    scores = model.fit_transform(values)
    second = scores[:, 1] if n_components > 1 else np.zeros(len(scores))
    explained = list(model.explained_variance_ratio_) + [0.0, 0.0]
    return pd.DataFrame(
        {
            "row_id": matrix.columns.astype(str),
            "sample_id": matrix.columns.astype(str),
            "PC1": scores[:, 0],
            "PC2": second,
            "explained_PC1": explained[0],
            "explained_PC2": explained[1],
            "events_used": len(working),
        }
    )


def sample_similarity(matrix: pd.DataFrame) -> pd.DataFrame:
    if matrix.empty:
        return pd.DataFrame(columns=["row_id", "sample_a", "sample_b", "spearman_r", "shared_events"])
    rows = []
    for index_a, sample_a in enumerate(matrix.columns):
        for index_b, sample_b in enumerate(matrix.columns):
            aligned = pd.concat(
                [matrix.iloc[:, index_a], matrix.iloc[:, index_b]], axis=1
            ).dropna()
            correlation, _ = safe_spearman(aligned.iloc[:, 0], aligned.iloc[:, 1])
            rows.append(
                {
                    "row_id": f"{sample_a}::{sample_b}",
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "spearman_r": correlation,
                    "shared_events": len(aligned),
                }
            )
    return pd.DataFrame(rows)


def load_deg_results(context: ProjectContext, contrast_ids: list[str]) -> pd.DataFrame:
    frames = []
    for asset in context.manifest.assets:
        if asset.kind != "deg_result" or not asset.path.is_file():
            continue
        contrast_id = str(asset.metadata.get("contrast_id", ""))
        if contrast_id not in contrast_ids:
            continue
        frame = pd.read_csv(asset.path, sep="\t", low_memory=False)
        gene_column = next((column for column in ("gene_id", "GeneID", "gene") if column in frame.columns), None)
        effect_column = next((column for column in ("log2FoldChange", "log2FC", "effect") if column in frame.columns), None)
        fdr_column = next((column for column in ("padj", "FDR", "qvalue") if column in frame.columns), None)
        if gene_column is None or effect_column is None or fdr_column is None:
            raise ValueError(f"DEG file has no recognized gene/effect/FDR columns: {asset.path}")
        output = pd.DataFrame(
            {
                "gene_id": frame[gene_column].map(normalize_gene_id),
                "gene_symbol": frame.get("gene_name", frame.get("gene_symbol", frame[gene_column])).map(clean_identifier),
                "log2fc": pd.to_numeric(frame[effect_column], errors="coerce"),
                "deg_fdr": pd.to_numeric(frame[fdr_column], errors="coerce"),
            }
        )
        if bool(asset.metadata.get("reverse_direction", False)):
            output["log2fc"] = -output["log2fc"]
        output["comparison_id"] = contrast_id
        output["row_id"] = contrast_id + "::" + output["gene_id"]
        output["source_path"] = str(asset.path)
        frames.append(output)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def load_dtu_results(context: ProjectContext, contrast_ids: list[str]) -> pd.DataFrame:
    frames = []
    for asset in context.manifest.assets:
        if asset.kind not in {"dtu_result", "dexseq_result"} or not asset.path.is_file():
            continue
        contrast_id = str(asset.metadata.get("contrast_id", ""))
        if contrast_id not in contrast_ids:
            continue
        frame = pd.read_csv(asset.path, sep="\t", low_memory=False)
        gene_column = next((column for column in ("groupID", "geneID", "gene_id", "GeneID") if column in frame.columns), None)
        transcript_column = next((column for column in ("featureID", "txID", "transcript", "transcript_id") if column in frame.columns), None)
        fdr_column = next((column for column in ("padj", "qvalue", "adj_pval", "FDR") if column in frame.columns), None)
        is_stager = {"gene", "transcript"}.issubset(frame.columns)
        if gene_column is None or (fdr_column is None and not is_stager):
            continue
        if fdr_column is not None:
            adjusted_p = pd.to_numeric(frame[fdr_column], errors="coerce")
            adjusted_p_scope = fdr_column
        else:
            # stageR exports gene- and transcript-confirmation adjusted p-values
            # in columns named `gene` and `transcript`. Prefer the transcript
            # value when available and retain the gene screen value otherwise.
            transcript_adjusted = pd.to_numeric(frame["transcript"], errors="coerce")
            gene_adjusted = pd.to_numeric(frame["gene"], errors="coerce")
            adjusted_p = transcript_adjusted.fillna(gene_adjusted)
            adjusted_p_scope = "stageR transcript confirmation, falling back to gene screen"
        output = pd.DataFrame(
            {
                "gene_id": frame[gene_column].map(normalize_gene_id),
                "transcript_id": frame[transcript_column].map(clean_identifier) if transcript_column else "",
                "dtu_fdr": adjusted_p,
            }
        )
        output["comparison_id"] = contrast_id
        output["source_kind"] = asset.kind
        output["adjusted_p_scope"] = adjusted_p_scope
        output["row_id"] = contrast_id + "::" + output["gene_id"] + "::" + output["transcript_id"]
        output["source_path"] = str(asset.path)
        frames.append(output)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def cross_comparison_summary(events: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    contrasts = list(dict.fromkeys(events["comparison_id"].astype(str)))
    pair_rows = []
    membership_rows = []
    tested = {
        contrast: set(events.loc[events["comparison_id"].eq(contrast), "event_id"].astype(str))
        for contrast in contrasts
    }
    significant = {
        contrast: set(
            events.loc[
                events["comparison_id"].eq(contrast) & events["stable_significant"],
                "event_id",
            ].astype(str)
        )
        for contrast in contrasts
    }
    for index, comparison_a in enumerate(contrasts):
        for comparison_b in contrasts[index + 1 :]:
            universe = tested[comparison_a] & tested[comparison_b]
            set_a = significant[comparison_a] & universe
            set_b = significant[comparison_b] & universe
            shared = set_a & set_b
            lookup_a = (
                events.loc[events["comparison_id"].eq(comparison_a)]
                .drop_duplicates("event_id")
                .set_index("event_id")
            )
            lookup_b = (
                events.loc[events["comparison_id"].eq(comparison_b)]
                .drop_duplicates("event_id")
                .set_index("event_id")
            )
            effects_a = [float(lookup_a.loc[event_id, "dpsi"]) for event_id in shared]
            effects_b = [float(lookup_b.loc[event_id, "dpsi"]) for event_id in shared]
            same_direction = sum(np.sign(left) == np.sign(right) for left, right in zip(effects_a, effects_b, strict=True))
            correlation, correlation_p = safe_spearman(effects_a, effects_b)
            pair_id = f"{comparison_a}__vs__{comparison_b}"
            pair_rows.append(
                {
                    "row_id": pair_id,
                    "comparison_a": comparison_a,
                    "comparison_b": comparison_b,
                    "tested_in_both": len(universe),
                    "significant_a": len(set_a),
                    "significant_b": len(set_b),
                    "shared_exact_events": len(shared),
                    "a_only_events": len(set_a - set_b),
                    "b_only_events": len(set_b - set_a),
                    "jaccard": len(shared) / len(set_a | set_b) if set_a | set_b else 0.0,
                    "containment_a_in_b": len(shared) / len(set_a) if set_a else math.nan,
                    "containment_b_in_a": len(shared) / len(set_b) if set_b else math.nan,
                    "overlap_pvalue": hypergeometric_overlap(len(universe), len(set_a), len(set_b), len(shared)),
                    "same_direction_events": same_direction,
                    "direction_concordance": same_direction / len(shared) if shared else math.nan,
                    "effect_spearman_r": correlation,
                    "effect_spearman_pvalue": correlation_p,
                    "identity_rule": "exact normalized rMATS event structure; effect direction is evaluated separately",
                }
            )
            for event_id in sorted(set_a | set_b):
                in_a = event_id in set_a
                in_b = event_id in set_b
                row_a = lookup_a.loc[event_id] if event_id in lookup_a.index else None
                row_b = lookup_b.loc[event_id] if event_id in lookup_b.index else None
                dpsi_a = float(row_a["dpsi"]) if row_a is not None else math.nan
                dpsi_b = float(row_b["dpsi"]) if row_b is not None else math.nan
                membership_rows.append(
                    {
                        "row_id": f"{pair_id}::{event_id}",
                        "pair_id": pair_id,
                        "event_id": event_id,
                        "gene_id": clean_identifier((row_a if row_a is not None else row_b).get("gene_id", "")),
                        "gene_symbol": clean_identifier((row_a if row_a is not None else row_b).get("gene_symbol", "")),
                        "event_type": clean_identifier((row_a if row_a is not None else row_b).get("event_type", "")),
                        "membership": "shared" if in_a and in_b else ("A specific" if in_a else "B specific"),
                        "dpsi_a": dpsi_a,
                        "dpsi_b": dpsi_b,
                        "direction_relation": (
                            "same direction"
                            if in_a and in_b and np.sign(dpsi_a) == np.sign(dpsi_b)
                            else "opposite direction"
                            if in_a and in_b
                            else "not shared"
                        ),
                    }
                )
    summary = pd.DataFrame(pair_rows)
    if not summary.empty:
        summary["overlap_fdr"] = bh_adjust(summary["overlap_pvalue"])
        summary["effect_spearman_fdr"] = bh_adjust(summary["effect_spearman_pvalue"])
    return summary, pd.DataFrame(membership_rows)


def build_candidate_evidence(
    events: pd.DataFrame,
    deg: pd.DataFrame,
    dtu: pd.DataFrame,
    parameters: dict[str, Any],
) -> pd.DataFrame:
    fdr = float(parameters.get("fdr", 0.05))
    deg_fc = float(parameters.get("deg_min_abs_log2fc", 1.0))
    event_summary = (
        events.groupby(["comparison_id", "gene_id", "gene_symbol"], dropna=False)
        .agg(
            tested_as_events=("event_id", "nunique"),
            significant_as_events=("significant", "sum"),
            stable_significant_as_events=("stable_significant", "sum"),
            strongest_abs_dpsi=("abs_dpsi", "max"),
            best_as_fdr=("fdr", "min"),
        )
        .reset_index()
    )
    if not deg.empty:
        deg_work = deg.copy()
        deg_work["deg_significant"] = deg_work["deg_fdr"].le(fdr) & deg_work["log2fc"].abs().ge(deg_fc)
        deg_summary = (
            deg_work.sort_values("deg_fdr")
            .drop_duplicates(["comparison_id", "gene_id"])
            [["comparison_id", "gene_id", "log2fc", "deg_fdr", "deg_significant"]]
        )
        event_summary = event_summary.merge(deg_summary, on=["comparison_id", "gene_id"], how="outer")
    if not dtu.empty:
        dtu_summary = (
            dtu.groupby(["comparison_id", "gene_id"], dropna=False)
            .agg(best_dtu_fdr=("dtu_fdr", "min"), significant_dtu_transcripts=("dtu_fdr", lambda values: int(pd.to_numeric(values, errors="coerce").le(fdr).sum())))
            .reset_index()
        )
        event_summary = event_summary.merge(dtu_summary, on=["comparison_id", "gene_id"], how="outer")
    frame = event_summary.copy()

    def column_or_default(name: str, default: object) -> pd.Series:
        if name in frame.columns:
            return frame[name]
        return pd.Series(default, index=frame.index)

    frame["gene_symbol"] = column_or_default("gene_symbol", "").fillna("")
    frame["as_supported"] = column_or_default("stable_significant_as_events", 0).fillna(0).gt(0)
    frame["deg_supported"] = column_or_default("deg_significant", False).map(
        lambda value: False if pd.isna(value) else bool(value)
    )
    frame["dtu_supported"] = column_or_default("significant_dtu_transcripts", 0).fillna(0).gt(0)
    frame["evidence_dimensions"] = frame[["as_supported", "deg_supported", "dtu_supported"]].sum(axis=1)
    frame["evidence_class"] = np.select(
        [
            frame["evidence_dimensions"].ge(2),
            frame["as_supported"],
            frame["deg_supported"] | frame["dtu_supported"],
        ],
        ["convergent multi-omic", "stable AS only", "expression/isoform only"],
        default="no threshold-passing evidence",
    )
    frame["row_id"] = frame["comparison_id"].fillna("").astype(str) + "::" + frame["gene_id"].fillna("").astype(str)
    return frame.sort_values(
        ["evidence_dimensions", "stable_significant_as_events", "strongest_abs_dpsi"],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)
