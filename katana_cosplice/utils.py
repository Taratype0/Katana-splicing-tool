from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def normalize_gene_id(value: object) -> str:
    text = str(value or "").strip().strip('"')
    return re.sub(r"\.\d+$", "", text)


def normalize_sample_id(value: object) -> str:
    name = Path(str(value or "").strip()).name
    for suffix in (".markdup.sorted.bam", "_sorted.bam", ".sorted.bam", ".bam"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_vector(value: object) -> list[float]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    output: list[float] = []
    for token in str(value).split(","):
        token = token.strip()
        if token in {"", "NA", "nan", "NaN"}:
            output.append(float("nan"))
        else:
            try:
                output.append(float(token))
            except ValueError:
                output.append(float("nan"))
    return output


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=float)
    result = np.full(array.shape, np.nan, dtype=float)
    valid = np.isfinite(array)
    if not valid.any():
        return result
    p = array[valid]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def dataframe_digest(frame: pd.DataFrame) -> str:
    hashed = pd.util.hash_pandas_object(frame, index=True).values.tobytes()
    return hashlib.sha256(hashed).hexdigest()


def file_signature(path: Path | None) -> dict[str, object]:
    if path is None:
        return {"path": "", "exists": False}
    resolved = path.resolve()
    if not resolved.exists():
        return {"path": str(resolved), "exists": False}
    stat = resolved.stat()
    return {"path": str(resolved), "exists": True, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def fingerprint(payload: dict[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
