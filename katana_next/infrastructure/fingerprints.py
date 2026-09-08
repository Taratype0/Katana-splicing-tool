from __future__ import annotations

import hashlib
import json
from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel

from katana_next.domain.models import FileFingerprint


_CHUNK_BYTES = 1024 * 1024
_FULL_HASH_LIMIT = 64 * 1024 * 1024


def _json_default(value: object) -> object:
    if isinstance(value, BaseModel):
        return value.model_dump(mode="json")
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Cannot fingerprint value of type {type(value).__name__}")


def canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        default=_json_default,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )


def fingerprint_value(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def fingerprint_file(path: str | Path, *, force_full: bool = False) -> FileFingerprint:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        return FileFingerprint(path=resolved, exists=False)

    stat = resolved.stat()
    digest = hashlib.sha256()
    if force_full or stat.st_size <= _FULL_HASH_LIMIT:
        with resolved.open("rb") as handle:
            for chunk in iter(lambda: handle.read(4 * _CHUNK_BYTES), b""):
                digest.update(chunk)
        method = "sha256-full"
    else:
        offsets = (
            0,
            max(0, stat.st_size // 2 - _CHUNK_BYTES // 2),
            max(0, stat.st_size - _CHUNK_BYTES),
        )
        with resolved.open("rb") as handle:
            for offset in offsets:
                handle.seek(offset)
                digest.update(handle.read(_CHUNK_BYTES))
        digest.update(str(stat.st_size).encode("ascii"))
        method = "sha256-sampled-3x1MiB"

    return FileFingerprint(
        path=resolved,
        exists=True,
        size=stat.st_size,
        modified_ns=stat.st_mtime_ns,
        digest=digest.hexdigest(),
        hash_method=method,
    )


def workflow_fingerprint(
    *,
    workflow_id: str,
    workflow_version: str,
    project_fingerprint: str,
    inputs: list[FileFingerprint],
    parameters: dict[str, Any],
) -> str:
    return fingerprint_value(
        {
            "workflow_id": workflow_id,
            "workflow_version": workflow_version,
            "project_fingerprint": project_fingerprint,
            "inputs": [item.model_dump(mode="json") for item in inputs],
            "parameters": parameters,
        }
    )
