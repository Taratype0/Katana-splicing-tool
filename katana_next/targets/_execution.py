"""Small portable execution primitives; no project database or scheduler access."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import threading
import time
import uuid


PLOT_LOCK = threading.RLock()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def portable(value):
    """JSON with explicit nulls, retaining finite native floats without rounding."""
    if isinstance(value, dict):
        return {str(k): portable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, set)):
        return [portable(v) for v in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "tolist"):
        return portable(value.tolist())
    if hasattr(value, "isoformat"):
        return value.isoformat()
    try:
        import pandas as pd
        if pd.isna(value):
            return None
    except (ValueError, TypeError):
        pass
    return str(value)


def canonical(value):
    return json.dumps(portable(value), ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def fingerprint(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(portable(value), ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        write_json(temporary, value)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def bounded(root, relative):
    """Reject Windows/Unix escapes and symlink ancestors on either platform."""
    relative = str(relative)
    pure = PurePosixPath(relative)
    if (not relative or "\\" in relative or ":" in relative or pure.is_absolute()
            or any(part in {"..", ".", ""} for part in pure.parts)):
        raise ValueError(f"Unsafe portable member path: {relative!r}")
    root = Path(root).resolve()
    path = root.joinpath(*pure.parts)
    if not path.resolve().is_relative_to(root):
        raise ValueError(f"Member escapes its result root: {relative!r}")
    current = root
    for part in pure.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"Symlink result member is not portable: {relative!r}")
    return path


def ledger(root, *, exclude=()):
    root = Path(root)
    excluded = set(exclude)
    result = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        if path.is_symlink():
            raise ValueError(f"Symlink product is not accepted: {relative}")
        if path.is_file() and relative not in excluded:
            result[relative] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    return result


def verify_ledger(root, members):
    if not isinstance(members, dict) or not members:
        raise ValueError("Result has no verified product ledger")
    for relative, proof in members.items():
        path = bounded(root, relative)
        if not path.is_file():
            raise ValueError(f"Missing result member: {relative}")
        if path.stat().st_size != proof["bytes"] or sha256(path) != proof["sha256"]:
            raise ValueError(f"Result member hash mismatch: {relative}")


@contextmanager
def owned_claim(path, identity, *, wait_seconds=0):
    """Never infer that another process's claim is stale or remove its claim."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    claim = {"identity": identity, "token": token, "pid": os.getpid(), "created_at": utc_now()}
    deadline = time.monotonic() + wait_seconds
    while True:
        try:
            fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            break
        except FileExistsError:
            if time.monotonic() >= deadline:
                raise RuntimeError(f"Target output is claimed by another execution: {path.name}") from None
            time.sleep(.05)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(claim, stream)
        yield claim
    finally:
        try:
            current = json.loads(path.read_text(encoding="utf-8"))
            if current.get("token") == token:
                path.unlink()
        except (OSError, ValueError):
            pass
