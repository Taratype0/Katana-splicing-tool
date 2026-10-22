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
from katana_next.science.gene_filter import apply_gene_filter


EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
READ_SUPPORT_POLICY_VERSION = "bulk-all-samples-counts-v1"
DEFAULT_MIN_READS = 10
_BULK_EVENT_LOADER_VERSION = "4.0-all-declared-samples"
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
    if isinstance(value, str) and not value.strip():
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
            values.append(math.nan)
            continue
        try:
            number = float(token)
        except ValueError:
            values.append(math.nan)
            continue
        values.append(number if math.isfinite(number) else math.nan)
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


def _declared_sample_names(asset: Any, contrast: Any, arm: str) -> list[str]:
    """An asset may reorder an arm's native vector, but never reduce its scope."""
    key = f"{arm}_sample_ids"
    declared = [str(value) for value in getattr(contrast, arm).sample_ids]
    if not declared or len(declared) != len(set(declared)):
        raise ValueError(f"{contrast.contrast_id}: {arm} requires nonempty unique declared samples")
    names = asset.metadata[key] if key in asset.metadata else declared
    if not isinstance(names, (list, tuple)):
        raise ValueError(f"{contrast.contrast_id}: {key} must map every declared sample")
    normalized = [str(value) for value in names]
    if len(normalized) != len(set(normalized)) or len(normalized) != len(declared) or set(normalized) != set(declared):
        raise ValueError(f"{contrast.contrast_id}: {key} must cover every declared {arm} sample exactly once; partial, duplicate or unmapped IDs are invalid")
    return normalized


def _sample_names(asset: Any, contrast: Any, arm: str, count: int) -> list[str]:
    normalized = _declared_sample_names(asset, contrast, arm)
    if len(normalized) < count:
        raise ValueError(f"{contrast.contrast_id}: {count} {arm} replicate values but only {len(normalized)} mapped samples. Correct the sample mapping before analysis.")
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
        source = apply_gene_filter(context, source)
        required = {"FDR", "IncLevelDifference", "IncLevel1", "IncLevel2"}
        missing = sorted(required - set(source.columns))
        if missing:
            raise ValueError(f"rMATS file {asset.path.name} is missing columns: {missing}")
        if source.empty:
            continue
        contrast = _contrast(context, contrast_id)
        frame = source.copy()
        for arm in ("numerator", "denominator"):
            frame[f"{arm}_expected_samples"] = len(_declared_sample_names(asset, contrast, arm))
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
        frame["source_row"] = frame.index.to_numpy(dtype=int)

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
        for arm_index, arm in ((1, arm1), (2, "denominator" if arm1 == "numerator" else "numerator")):
            expected_count = len(_declared_sample_names(asset, contrast, arm))
            for column in (f"IncLevel{arm_index}", f"IJC_SAMPLE_{arm_index}", f"SJC_SAMPLE_{arm_index}"):
                if not parsed[column].map(len).isin({0, expected_count}).all():
                    raise ValueError(f"{asset.path.name}: {column} replicate count does not match its {expected_count} declared samples. NA positions must be retained.")
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
                "median_psi": _finite_median(values),
                "max_absolute_disagreement": _finite_range(values),
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
        "loader_version": _BULK_EVENT_LOADER_VERSION,
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
        metadata = json.loads(paths["metadata"].read_text(encoding="utf-8"))
        if metadata.get("loader_version") != _BULK_EVENT_LOADER_VERSION:
            return None
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
                "schema_version": "3.0",
                "loader_version": _BULK_EVENT_LOADER_VERSION,
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
        sum((difference > 0) == expected for difference in differences) / len(differences)
        if differences and math.isfinite(dpsi) and dpsi != 0
        else math.nan
    )
    loo_effects = []
    for group, other, sign in (
        (numerator_psi, denominator_psi, 1.0),
        (denominator_psi, numerator_psi, -1.0),
    ):
        group = [value for value in group if math.isfinite(value)]
        other = [value for value in other if math.isfinite(value)]
        for index in range(len(group)):
            reduced = [value for position, value in enumerate(group) if position != index]
            if reduced and other:
                loo_effects.append(sign * (sum(reduced) / len(reduced) - sum(other) / len(other)))
    loo_effects = [effect for effect in loo_effects if math.isfinite(effect) and effect != 0]
    loo_consistency = (
        sum((effect > 0) == expected for effect in loo_effects) / len(loo_effects)
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
        "numerator_replicates": sum(math.isfinite(value) for value in numerator_psi),
        "denominator_replicates": sum(math.isfinite(value) for value in denominator_psi),
        "numerator_mean_reads": sum(numerator_reads) / len(numerator_reads) if numerator_reads else math.nan,
        "denominator_mean_reads": sum(denominator_reads) / len(denominator_reads) if denominator_reads else math.nan,
        "replicate_pairwise_direction_consistency": pairwise,
        "leave_one_out_direction_consistency": loo_consistency,
        "max_measurement_se": max(uncertainty) if uncertainty else math.nan,
    }


def _observation_matrix(values: dict[tuple[str, str], list[float]]) -> pd.DataFrame:
    if not values:
        return pd.DataFrame()
    matrix = pd.Series({key: _finite_median(observations) for key, observations in values.items()}, dtype=float).unstack()
    matrix = matrix.sort_index().sort_index(axis=1)
    matrix.index.name = "event_id"
    return matrix


def _finite_median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    n = len(finite)
    return (finite[n // 2] if n % 2 else (finite[n // 2 - 1] + finite[n // 2]) / 2) if n else math.nan


def _finite_range(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return max(finite) - min(finite) if finite else math.nan


def _validate_read_support_parameters(parameters: dict[str, Any]) -> str:
    """Reject changed legacy gates; their old default values are inert compatibility input."""
    legacy_defaults = {
        "stability_profile": "balanced",
        "min_supported_replicates": 2,
        "min_direction_consistency": 0.75,
        "max_measurement_se": 0.12,
        "apply_measurement_se_gate": False,
    }
    if "min_reads" in parameters and parameters["min_reads"] != DEFAULT_MIN_READS:
        raise ValueError(f"{READ_SUPPORT_POLICY_VERSION} requires fixed min_reads={DEFAULT_MIN_READS} in every declared sample")
    requested = parameters.get("read_support_policy", READ_SUPPORT_POLICY_VERSION)
    if requested != READ_SUPPORT_POLICY_VERSION:
        raise ValueError(f"Unknown read-support policy: {requested}")
    ignored = {}
    for name, previous_default in legacy_defaults.items():
        if name not in parameters:
            continue
        if parameters[name] != previous_default:
            raise ValueError(f"Legacy QC option {name}={parameters[name]!r} is incompatible with {READ_SUPPORT_POLICY_VERSION}; direction, replicate-minimum and SE gates are no longer applied")
        ignored[name] = parameters[name]
    return json.dumps(ignored, sort_keys=True)


def _policy_vector(value: object) -> list[float]:
    if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
        result = []
        for item in value:
            try:
                result.append(float(item) if not isinstance(item, (bool, np.bool_)) else math.nan)
            except (TypeError, ValueError):
                result.append(math.nan)
        return result
    return parse_number_list(value)


def _arm_read_support(row: pd.Series, arm: str, suffix: str) -> dict[str, Any]:
    psi = _policy_vector(row.get(f"IncLevel{suffix}"))
    ijc = _policy_vector(row.get(f"IJC_SAMPLE_{suffix}"))
    sjc = _policy_vector(row.get(f"SJC_SAMPLE_{suffix}"))
    key = f"{arm}_expected_samples"
    expected = len(psi)
    reasons = []
    basis = "PSI_vector_length" if key not in row else "declared_contrast_samples"
    if key in row:
        value = row[key]
        try:
            if isinstance(value, (bool, np.bool_)) or not math.isfinite(float(value)) or float(value) != int(value) or int(value) < 1:
                raise ValueError
            expected = int(value)
        except (TypeError, ValueError, OverflowError):
            expected = 0
            reasons.append(f"{arm}:invalid_declared_sample_count")
    if expected < 1:
        reasons.append(f"{arm}:empty_arm")
    for label, values in (("PSI", psi), ("IJC", ijc), ("SJC", sjc)):
        if len(values) != expected:
            reasons.append(f"{arm}:{label}_vector_length_{len(values)}_expected_{expected}")
    shape_valid = expected > 0 and len(psi) == len(ijc) == len(sjc) == expected
    count_shape_valid = expected > 0 and len(ijc) == len(sjc) == expected
    usable = read_supported = supported = invalid = low_reads = 0
    observed_positions = min(expected, max(len(psi), len(ijc), len(sjc)))
    for index in range(observed_positions):
        p = psi[index] if index < len(psi) else math.nan
        i = ijc[index] if index < len(ijc) else math.nan
        s = sjc[index] if index < len(sjc) else math.nan
        psi_valid = math.isfinite(p) and 0 <= p <= 1
        counts_valid = math.isfinite(i) and math.isfinite(s) and i >= 0 and s >= 0 and math.isfinite(i + s)
        read_pass = counts_valid and i + s >= DEFAULT_MIN_READS
        usable += int(psi_valid and counts_valid)
        read_supported += int(read_pass)
        supported += int(psi_valid and read_pass)
        invalid += int(not (psi_valid and counts_valid))
        low_reads += int(counts_valid and not read_pass)
        if not psi_valid:
            reasons.append(f"{arm}:sample_{index + 1}:invalid_PSI")
        if not counts_valid:
            reasons.append(f"{arm}:sample_{index + 1}:invalid_counts")
        elif not read_pass:
            reasons.append(f"{arm}:sample_{index + 1}:counts_below_{DEFAULT_MIN_READS}")
    if observed_positions < expected:
        invalid += expected - observed_positions
        reasons.append(f"{arm}:missing_sample_positions_{expected - observed_positions}")
    return {
        f"{arm}_expected_samples": expected,
        f"{arm}_sample_scope_basis": basis,
        f"{arm}_observed_psi_samples": len(psi),
        f"{arm}_observed_ijc_samples": len(ijc),
        f"{arm}_observed_sjc_samples": len(sjc),
        f"{arm}_read_count_supported_samples": read_supported,
        f"{arm}_supported_replicates": supported,
        f"{arm}_invalid_samples": invalid,
        f"{arm}_low_read_samples": low_reads,
        f"{arm}_vector_lengths_match": shape_valid,
        f"{arm}_data_usable": shape_valid and usable == expected,
        f"{arm}_all_read_counts_pass": count_shape_valid and read_supported == expected,
        "_reasons": reasons,
    }


def apply_stable_event_qc(events: pd.DataFrame, parameters: dict[str, Any]) -> pd.DataFrame:
    """Compatibility entry for fixed ten-count, all-declared-sample eligibility.

    Every sample in each selected contrast arm is checked independently. Native
    PSI, count, P and FDR columns are never rewritten. Existing direction/SE
    diagnostics are preserved as observations and do not select events.
    """
    ignored = _validate_read_support_parameters(parameters)
    retired = {
        "stability_profile", "differential_repeat_supported",
        "passes_direction_consistency", "passes_leave_one_out",
        "passes_measurement_uncertainty", "criterion_min_supported_replicates",
        "criterion_min_direction_consistency", "criterion_max_measurement_se",
    }
    frame = events.drop(columns=list(retired), errors="ignore").copy()
    audits = []
    for _, row in frame.iterrows():
        reverse = bool(row.get("reverse_direction", False))
        numerator = _arm_read_support(row, "numerator", "2" if reverse else "1")
        denominator = _arm_read_support(row, "denominator", "1" if reverse else "2")
        reasons = numerator.pop("_reasons") + denominator.pop("_reasons")
        audit = {**numerator, **denominator}
        audit["passes_replicate_count"] = numerator["numerator_vector_lengths_match"] and denominator["denominator_vector_lengths_match"]
        audit["passes_data_usability"] = numerator["numerator_data_usable"] and denominator["denominator_data_usable"]
        audit["passes_read_counts"] = numerator["numerator_all_read_counts_pass"] and denominator["denominator_all_read_counts_pass"]
        audit["passes_read_support"] = audit["passes_data_usability"] and audit["passes_read_counts"]
        audit["measurement_eligible"] = audit["passes_read_support"]
        audit["stable_event"] = audit["measurement_eligible"]
        audit["read_support_failure_reason"] = "; ".join(reasons) or "passed"
        audit["stability_failure_reason"] = audit["read_support_failure_reason"]
        audits.append(audit)
    audit_frame = pd.DataFrame(audits, index=frame.index)
    for column in audit_frame:
        frame[column] = audit_frame[column]
    if frame.empty:
        for column in ("passes_replicate_count", "passes_data_usability", "passes_read_counts", "passes_read_support", "measurement_eligible", "stable_event"):
            frame[column] = pd.Series(index=frame.index, dtype=bool)
        for column in ("read_support_failure_reason", "stability_failure_reason"):
            frame[column] = pd.Series(index=frame.index, dtype=object)
    frame["read_support_policy"] = READ_SUPPORT_POLICY_VERSION
    frame["criterion_min_reads"] = DEFAULT_MIN_READS
    frame["measurement_se_gate_enabled"] = False
    frame["ignored_legacy_qc_options"] = ignored
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
        frame = apply_gene_filter(context, frame)
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
        frame = apply_gene_filter(context, frame)
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
            # value only; a gene screen cannot confirm a specific transcript.
            transcript_adjusted = pd.to_numeric(frame["transcript"], errors="coerce")
            gene_adjusted = pd.to_numeric(frame["gene"], errors="coerce")
            adjusted_p = transcript_adjusted
            adjusted_p_scope = "stageR transcript confirmation only"
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
        output["gene_screen_fdr"] = pd.to_numeric(frame["gene"], errors="coerce") if is_stager else float("nan")
        output["row_id"] = contrast_id + "::" + output["gene_id"] + "::" + output["transcript_id"]
        output["source_path"] = str(asset.path)
        frames.append(output)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def cross_comparison_summary(events: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    contrasts = list(dict.fromkeys(events["comparison_id"].astype(str)))
    pair_rows = []
    membership_rows = []
    present = {
        contrast: set(events.loc[events["comparison_id"].eq(contrast), "event_id"].astype(str))
        for contrast in contrasts
    }
    testable = (np.isfinite(pd.to_numeric(events["fdr"], errors="coerce"))
                if "fdr" in events else pd.Series(True, index=events.index))
    testable &= np.isfinite(pd.to_numeric(events["dpsi"], errors="coerce"))
    tested = {contrast: set(events.loc[events["comparison_id"].eq(contrast) & testable, "event_id"].astype(str))
              for contrast in contrasts}
    significant = {
        contrast: set(
            events.loc[
                events["comparison_id"].eq(contrast) & events["stable_significant"] & testable,
                "event_id",
            ].astype(str)
        )
        for contrast in contrasts
    }
    for index, comparison_a in enumerate(contrasts):
        for comparison_b in contrasts[index + 1 :]:
            universe = tested[comparison_a] & tested[comparison_b]
            set_a = significant[comparison_a]
            set_b = significant[comparison_b]
            test_a, test_b = set_a & universe, set_b & universe
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
                    "significant_a_testable_in_both": len(test_a),
                    "significant_b_testable_in_both": len(test_b),
                    "a_significant_missing_in_b": len(set_a - present[comparison_b]),
                    "b_significant_missing_in_a": len(set_b - present[comparison_a]),
                    "a_significant_untestable_in_b": len(set_a & (present[comparison_b] - tested[comparison_b])),
                    "b_significant_untestable_in_a": len(set_b & (present[comparison_a] - tested[comparison_a])),
                    "shared_exact_events": len(shared),
                    "a_only_events": len(set_a - set_b),
                    "b_only_events": len(set_b - set_a),
                    "jaccard": len(shared) / len(set_a | set_b) if set_a | set_b else math.nan,
                    "containment_a_in_b": len(shared) / len(set_a) if set_a else math.nan,
                    "containment_b_in_a": len(shared) / len(set_b) if set_b else math.nan,
                    "testable_containment_a_in_b": len(shared) / len(test_a) if test_a else math.nan,
                    "testable_containment_b_in_a": len(shared) / len(test_b) if test_b else math.nan,
                    "overlap_pvalue": hypergeometric_overlap(len(universe), len(test_a), len(test_b), len(shared)) if universe else math.nan,
                    "containment_denominator": "complete significant set, including counterpart-missing events",
                    "overlap_null": "exchangeable event sets within jointly tested universe; shared-control dependence is not calibrated",
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
                        "counterpart_status": (
                            "both significant" if in_a and in_b else
                            "not tested in B" if row_b is None else
                            "not tested in A" if row_a is None else
                            "present but not testable in B" if event_id not in tested[comparison_b] else
                            "present but not testable in A" if event_id not in tested[comparison_a] else
                            "tested but not in strict significant set"
                        ),
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


def usage_gene_support(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep exon- and transcript-usage minima separate; neither is a gene test."""
    keys = ["comparison_id", "gene_id"]
    result = frame[keys].drop_duplicates().copy()
    kinds = frame.get("source_kind", pd.Series("dtu_result", index=frame.index))
    for kind, name in (("dtu_result", "best_dtu_fdr"), ("dexseq_result", "best_dexseq_fdr")):
        subset = frame.loc[kinds.eq(kind), [*keys, "dtu_fdr"]].copy()
        subset["dtu_fdr"] = pd.to_numeric(subset["dtu_fdr"], errors="coerce")
        values = subset.groupby(keys, as_index=False).agg(**{name: ("dtu_fdr", "min")})
        result = result.merge(values, on=keys, how="left", validate="one_to_one")
    return result


def candidate_display_labels(frame: pd.DataFrame, labels: dict[str, str]) -> pd.Series:
    names = frame["gene_symbol"].where(frame["gene_symbol"].ne(""), frame["gene_id"])
    display = names + " | " + frame["comparison_id"].map(labels).fillna(frame["comparison_id"])
    duplicate = display.duplicated(keep=False)
    return display.where(~duplicate, display + " [" + frame["row_id"] + "]")


def candidate_pair_plot_counts(pairs: pd.DataFrame, labels: dict[str, str]) -> pd.DataFrame:
    counts = pairs.groupby(["comparison_A", "comparison_B", "membership"], observed=True).size().reset_index(name="genes")
    label_counts = pd.Series(labels).value_counts()
    display = {key: value if label_counts[value] == 1 else f"{value} [{key}]" for key, value in labels.items()}
    for side in ("A", "B"):
        key = f"comparison_{side}"
        counts[key + "_label"] = counts[key].map(display).fillna(counts[key])
    counts["pair"] = counts.comparison_A_label + " compared with " + counts.comparison_B_label
    return counts


def build_candidate_evidence(
    events: pd.DataFrame,
    deg: pd.DataFrame,
    dtu: pd.DataFrame,
    parameters: dict[str, Any],
) -> pd.DataFrame:
    fdr = float(parameters.get("fdr", 0.05))
    deg_fc = float(parameters.get("deg_min_abs_log2fc", 1.0))
    if events.empty:
        events = pd.DataFrame(columns=["comparison_id", "gene_id", "gene_symbol", "event_id", "significant", "stable_significant", "abs_dpsi", "fdr"])
    events = events.copy()
    events["strict_abs_dpsi"] = events["abs_dpsi"].where(events["stable_significant"].fillna(False).astype(bool))
    event_summary = (
        events.groupby(["comparison_id", "gene_id"], dropna=False)
        .agg(
            gene_symbol=("gene_symbol", "first"),
            tested_as_events=("event_id", "nunique"),
            significant_as_events=("significant", "sum"),
            stable_significant_as_events=("stable_significant", "sum"),
            strongest_abs_dpsi=("strict_abs_dpsi", "max"),
            strongest_tested_abs_dpsi=("abs_dpsi", "max"),
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
        kinds = dtu.get("source_kind", pd.Series("dtu_result", index=dtu.index))
        for prefix, kind, feature in (("dtu", "dtu_result", "transcripts"), ("dexseq", "dexseq_result", "exons")):
            evidence = dtu.loc[kinds.eq(kind)].copy()
            if evidence.empty:
                continue
            evidence["dtu_fdr"] = pd.to_numeric(evidence["dtu_fdr"], errors="coerce")
            if "transcript_id" in evidence:
                evidence = evidence.sort_values("dtu_fdr").drop_duplicates(["comparison_id", "gene_id", "transcript_id"])
            summary = evidence.groupby(["comparison_id", "gene_id"], dropna=False).agg(
                **{f"best_{prefix}_fdr": ("dtu_fdr", "min"),
                   f"significant_{prefix}_{feature}": ("dtu_fdr", lambda values: int(values.le(fdr).sum()))}
            ).reset_index()
            event_summary = event_summary.merge(summary, on=["comparison_id", "gene_id"], how="outer")
    if not events.empty and "dpsi" in events:
        dominant = events.loc[events["stable_significant"].fillna(False)].sort_values(
            ["abs_dpsi", "event_id"], ascending=[False, True]).drop_duplicates(["comparison_id", "gene_id"])
        dominant = dominant[["comparison_id", "gene_id", "event_id", "dpsi"]].rename(
            columns={"event_id": "dominant_as_event_id", "dpsi": "dominant_as_dpsi"})
        event_summary = event_summary.merge(dominant, on=["comparison_id", "gene_id"], how="left")
    frame = event_summary.copy()

    def column_or_default(name: str, default: object) -> pd.Series:
        if name in frame.columns:
            return frame[name]
        return pd.Series(default, index=frame.index)

    frame["gene_symbol"] = column_or_default("gene_symbol", "").fillna("")
    for evidence in (deg, dtu):
        if not evidence.empty and "gene_symbol" in evidence:
            names = evidence.dropna(subset=["gene_symbol"]).drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"]
            frame["gene_symbol"] = frame["gene_symbol"].where(frame["gene_symbol"].ne(""), frame["gene_id"].map(names)).fillna("")
    frame["as_available"] = column_or_default("tested_as_events", np.nan).notna()
    frame["deg_available"] = column_or_default("deg_fdr", np.nan).notna()
    frame["dtu_available"] = column_or_default("best_dtu_fdr", np.nan).notna()
    frame["dexseq_available"] = column_or_default("best_dexseq_fdr", np.nan).notna()
    frame["as_supported"] = pd.to_numeric(column_or_default("stable_significant_as_events", 0), errors="coerce").fillna(0).gt(0)
    frame["deg_supported"] = column_or_default("deg_significant", False).map(
        lambda value: False if pd.isna(value) else bool(value)
    )
    frame["dtu_supported"] = column_or_default("significant_dtu_transcripts", 0).fillna(0).gt(0)
    frame["dexseq_supported"] = column_or_default("significant_dexseq_exons", 0).fillna(0).gt(0)
    frame["evidence_dimensions"] = frame[["as_supported", "deg_supported", "dtu_supported", "dexseq_supported"]].sum(axis=1)
    frame["evidence_class"] = np.select(
        [
            frame["evidence_dimensions"].ge(2),
            frame["as_supported"],
            frame["deg_supported"] | frame["dtu_supported"] | frame["dexseq_supported"],
        ],
        ["multiple supporting analyses", "stable AS only", "expression/usage only"],
        default="no threshold-passing evidence",
    )
    frame["row_id"] = frame["comparison_id"].fillna("").astype(str) + "::" + frame["gene_id"].fillna("").astype(str)
    frame = frame.sort_values(
        ["evidence_dimensions", "stable_significant_as_events", "strongest_abs_dpsi"],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    frame["rank_within_comparison"] = frame.groupby("comparison_id").cumcount() + 1
    frame["ranking_score"] = frame["evidence_dimensions"].astype(int)
    return frame


def pairwise_candidate_membership(candidates: pd.DataFrame, comparisons: list[str]) -> pd.DataFrame:
    """Gene membership and direction summaries; never an AS-event reversal test."""
    from itertools import combinations
    rows = []
    for a, b in combinations(comparisons, 2):
        sides = []
        for contrast, suffix in ((a, "A"), (b, "B")):
            side = candidates.loc[candidates.comparison_id.eq(contrast)].copy()
            if side.gene_id.duplicated().any():
                raise ValueError("Candidate table must have one row per gene and comparison")
            side = side.drop(columns=["row_id", "comparison_id"]).set_index("gene_id")
            sides.append(side.add_suffix(f"_{suffix}"))
        pair = sides[0].join(sides[1], how="outer").reset_index()
        if pair.empty:
            continue
        for suffix in ("A", "B"):
            pair[f"is_candidate_{suffix}"] = pair[f"evidence_dimensions_{suffix}"].fillna(0).gt(0)
        pair = pair.loc[pair.is_candidate_A | pair.is_candidate_B].copy()
        pair["comparison_A"], pair["comparison_B"] = a, b
        pair["membership"] = np.select([pair.is_candidate_A & pair.is_candidate_B, pair.is_candidate_A],
                                         ["shared", "A_specific"], default="B_specific")
        pair["n_comparisons_candidate"] = pair[["is_candidate_A", "is_candidate_B"]].sum(axis=1)
        for domain, effect in (("deg", "log2fc"), ("as", "dominant_as_dpsi")):
            def effect_column(suffix):
                return pd.to_numeric(pair.get(f"{effect}_{suffix}", pd.Series(np.nan, index=pair.index)), errors="coerce")
            left, right = effect_column("A"), effect_column("B")
            supported = pair[f"{domain}_supported_A"].eq(True) & pair[f"{domain}_supported_B"].eq(True)
            pair[f"{domain}_direction_comparable"] = supported & left.notna() & right.notna() & left.ne(0) & right.ne(0)
            pair[f"{domain}_direction_changed"] = pair[f"{domain}_direction_comparable"] & np.sign(left).ne(np.sign(right))
        pair["gene_level_direction_changed"] = pair.deg_direction_changed | pair.as_direction_changed
        pair["direction_pattern"] = np.select(
            [pair.gene_level_direction_changed, pair.deg_direction_comparable | pair.as_direction_comparable],
            ["gene_level_direction_changed", "comparable_directions_agree"], default="direction_not_comparable")
        pair["best_pair_rank"] = pair[["rank_within_comparison_A", "rank_within_comparison_B"]].min(axis=1)
        pair["row_id"] = a + "::" + b + "::" + pair.gene_id.astype(str)
        rows.append(pair)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["row_id", "comparison_A", "comparison_B", "gene_id", "membership", "direction_pattern"])
