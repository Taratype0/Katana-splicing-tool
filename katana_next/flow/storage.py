"""Bounded cleanup of work owned by one managed, successfully checkpointed run.

This module never cleans arbitrary Nextflow directories, shared caches, inputs
or results.  Callers keep the exclusive run lock while running Nextflow and
while validating the independent checkpoint.  Deleting a stage's work disables
Nextflow ``-resume`` for that stage; the independently verified outputs remain.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
from uuid import uuid4


SCOPE_SCHEMA = "katana.flow.managed_storage/v1"
CLEANUP_SCHEMA = "katana.flow.storage_cleanup/v1"
_STAGES = frozenset({"upstream", "downstream"})
_HELD_LOCKS: dict[str, tuple[str, tuple[int, int]]] = {}
_MARKER = ".katana-managed-run.json"
_LOCK = ".katana-managed-run.lock"


def _now():
    return datetime.now(timezone.utc).isoformat()


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     default=str, allow_nan=False).encode()).hexdigest()


def _signed(value):
    return {**value, "record_sha256": _hash(value)}


def _read(path):
    before = path.lstat()
    if _link(before):
        raise ValueError("Storage control records cannot be links")
    value = json.loads(path.read_text(encoding="utf-8"))
    after = path.lstat()
    if _identity(before) != _identity(after):
        raise ValueError("Storage control record changed while being read")
    if value.get("record_sha256") != _hash({k: v for k, v in value.items() if k != "record_sha256"}):
        raise ValueError("Storage control record identity is invalid")
    return value


def _write_new(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, allow_nan=False, default=str)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def _link(info):
    return stat.S_ISLNK(info.st_mode) or bool(getattr(info, "st_file_attributes", 0) &
        getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400))


def _identity(info):
    return (info.st_dev, info.st_ino, info.st_mode, info.st_size, info.st_mtime_ns,
            getattr(info, "st_file_attributes", 0))


def _node(info):
    return info.st_dev, info.st_ino, stat.S_IFMT(info.st_mode)


def _bounded_path(value):
    raw = Path(value).expanduser()
    if ".." in raw.parts:
        raise ValueError("Storage paths cannot contain parent traversal")
    location = raw.absolute()
    for parent in reversed([location, *location.parents]):
        if parent.exists() or parent.is_symlink():
            info = parent.lstat()
            if _link(info):
                raise ValueError(f"Storage roots/control paths cannot traverse links: {parent}")
    return location


def _overlaps(left, right):
    return left == right or left.is_relative_to(right) or right.is_relative_to(left)


def _scope(root):
    root = _bounded_path(root)
    value = _read(root / _MARKER)
    if value.get("schema") != SCOPE_SCHEMA or value.get("scope_root") != str(root):
        raise ValueError("This directory is not the exact registered managed scope")
    for item in value["protected_roots"]:
        if _overlaps(root, Path(item)):
            raise ValueError("Managed storage overlaps protected source/input/result roots")
    return root, value


def create_managed_scope(scope_root, *, protected_roots=(), run_id=None):
    """Register one new, empty scope outside protected input/result/code roots."""
    root = _bounded_path(scope_root)
    protected = {Path(__file__).resolve().parents[2], Path.home().resolve()}
    # A user home itself is forbidden; descendants outside the code tree are
    # allowed, so ordinary ~/runs/<name> locations continue to work.
    if root == Path(root.anchor) or root == Path.home().resolve() or len(root.parts) < 3:
        raise ValueError("Choose a dedicated run directory, not a filesystem/home root")
    protected.discard(Path.home().resolve())
    for value in protected_roots:
        path = Path(value).expanduser().absolute()
        protected.update((path, path.resolve()))
    if any(_overlaps(root, item) for item in protected):
        raise ValueError("Managed storage must be separate from protected roots")
    if root.exists() and any(root.iterdir()):
        raise FileExistsError("Managed scope must be new or empty; unknown directories are never adopted")
    root.mkdir(parents=True, exist_ok=True)
    marker = _signed({"schema": SCOPE_SCHEMA, "scope_root": str(root),
                      "run_id": str(run_id or uuid4()), "created_at": _now(),
                      "protected_roots": sorted(map(str, protected))})
    _write_new(root / _MARKER, marker)
    return marker


@contextmanager
def managed_run_lock(scope_root):
    """Exclusive lifetime lock; an existing/stale claim is never auto-removed."""
    root, owner = _scope(scope_root)
    token = str(uuid4())
    lock = root / _LOCK
    record = _signed({"scope_sha256": owner["record_sha256"], "token": token,
                      "pid": os.getpid(), "created_at": _now()})
    _write_new(lock, record)
    identity = _node(lock.lstat())
    _HELD_LOCKS[str(root)] = token, identity
    try:
        yield record
    finally:
        held = _HELD_LOCKS.pop(str(root), None)
        if held != (token, identity) or not lock.exists() or _node(lock.lstat()) != identity:
            raise RuntimeError("Run lock ownership changed; no lock was removed")
        current = _read(lock)
        if current.get("token") != token:
            raise RuntimeError("Run lock content changed; no lock was removed")
        lock.unlink()


def _assert_lock(root, owner):
    held = _HELD_LOCKS.get(str(root))
    if held is None:
        raise RuntimeError("Cleanup/registration requires the held exclusive managed-run lock")
    lock = root / _LOCK
    value = _read(lock)
    if (_node(lock.lstat()) != held[1] or value.get("token") != held[0] or
            value.get("scope_sha256") != owner["record_sha256"]):
        raise RuntimeError("Managed-run lock identity changed")


def _stage_name(stage):
    if stage not in _STAGES:
        raise ValueError("Managed storage stage must be upstream or downstream")
    return stage


def register_stage(scope_root, stage, work_root=None):
    """Own exactly scope/work/{stage}; never adopt a nonempty unowned work root."""
    root, owner = _scope(scope_root)
    _assert_lock(root, owner)
    stage = _stage_name(stage)
    expected = root / "work" / stage
    work = _bounded_path(work_root if work_root is not None else expected)
    if work != expected:
        raise ValueError("Stage work must be the exact scope/work/stage directory")
    control = root / f".katana-stage-{stage}.json"
    if control.exists():
        value = _read(control)
        if (value.get("scope_sha256") != owner["record_sha256"] or
                value.get("stage") != stage or value.get("work_root") != str(work)):
            raise ValueError("Registered stage identity differs from this run")
        if (root / f"cleanup-{stage}.pending.json").exists():
            raise RuntimeError("An earlier cleanup has an unresolved outcome; inspect its pending record, do not retry blindly")
        cleanup_record = root / f"cleanup-{stage}.receipt.json"
        if cleanup_record.exists():
            completed = _read(cleanup_record)
            if (work.exists() or completed.get("status") != "completed" or
                    completed.get("scope_sha256") != owner["record_sha256"] or
                    completed.get("stage") != stage):
                raise ValueError("Previously cleaned stage has unexpected new data or changed ownership")
        elif not work.is_dir() or list(_node(work.lstat())) != value["work_identity"]:
            raise ValueError("Registered work directory is missing or changed identity")
        return value
    if work.exists() and (not work.is_dir() or any(work.iterdir())):
        raise FileExistsError("Unowned/nonempty stage work is never adopted")
    work.mkdir(parents=True, exist_ok=True)
    value = _signed({"schema": SCOPE_SCHEMA + "/stage", "scope_sha256": owner["record_sha256"],
                     "stage": stage, "work_root": str(work), "created_at": _now(),
                     "work_identity": list(_node(work.lstat()))})
    _write_new(control, value)
    return value


def _inventory(work):
    rows = []
    stack = [work]
    while stack:
        path = stack.pop()
        info = path.lstat()
        link = _link(info)
        if not (link or stat.S_ISDIR(info.st_mode) or stat.S_ISREG(info.st_mode)):
            raise ValueError(f"Special/unrecognized file in stage work: {path}")
        rows.append({"relative_path": path.relative_to(work).as_posix(),
                     "identity": list(_identity(info)), "node": list(_node(info)),
                     "link": link, "directory": stat.S_ISDIR(info.st_mode) and not link,
                     "size": info.st_size, "links": info.st_nlink,
                     "allocated": getattr(info, "st_blocks", None) * 512 if hasattr(info, "st_blocks") else None})
        if stat.S_ISDIR(info.st_mode) and not link:
            stack.extend(path / name for name in sorted(os.listdir(path), reverse=True))
    return rows


def _sizes(rows):
    files = [row for row in rows if not row["directory"] and not row["link"]]
    counts = Counter(tuple(row["node"]) for row in files)
    unique = {tuple(row["node"]): row for row in files}
    known = all(row["allocated"] is not None for row in unique.values())
    return {"regular_files": len(files), "symlinks_or_reparse_points": sum(row["link"] for row in rows),
            "directories": sum(row["directory"] for row in rows),
            "logical_file_bytes": sum(row["size"] for row in files),
            "unique_allocated_bytes_estimate": sum(row["allocated"] for row in unique.values()) if known else None,
            "reclaimable_allocated_bytes_estimate": sum(row["allocated"] for key, row in unique.items()
                if row["links"] <= counts[key]) if known else None,
            "quota_reclaimed_bytes": None,
            "accounting_note": "Allocation is a stat-based estimate; external hardlinks are not counted as reclaimed. Filesystem/project quota return is not measured."}


def _verify_inventory(work, rows):
    current = _inventory(work)
    if current != rows:
        raise RuntimeError("Stage work changed after cleanup planning; no deletion is authorized")


def _remove_verified_tree(work, rows):
    """Unlink leaves without following staged links; check identities per entry."""
    nodes = {row["relative_path"]: row for row in rows}
    for row in sorted(rows, key=lambda item: len(Path(item["relative_path"]).parts), reverse=True):
        relative = Path(row["relative_path"])
        path = work if row["relative_path"] == "." else work / relative
        # Parent identities ignore directory mtimes, which change as their
        # already-verified child entries are removed.
        for parent in [path.parent, *path.parent.parents]:
            if parent == work.parent:
                break
            key = parent.relative_to(work).as_posix()
            expected = nodes[key]
            observed = parent.lstat()
            if _link(observed) or list(_node(observed)) != expected["node"]:
                raise RuntimeError("Cleanup parent changed identity; deletion stopped")
        observed = path.lstat()
        match = list(_node(observed)) == row["node"] if row["directory"] else list(_identity(observed)) == row["identity"]
        if not match or _link(observed) != row["link"]:
            raise RuntimeError("Cleanup member changed identity; deletion stopped")
        if row["directory"]:
            path.rmdir()
        elif row["link"] and stat.S_ISDIR(observed.st_mode) and not stat.S_ISLNK(observed.st_mode):
            path.rmdir()  # Windows directory junction itself, never its target.
        else:
            path.unlink()


def cleanup_stage(scope_root, stage, *, native_status, checkpoint_validator, dry_run=False):
    """Delete only registered work after native success and verified handoff.

    ``native_status`` must declare state/status COMPLETED and exit_code 0 or
    '0:0'. ``checkpoint_validator(work_root)`` must actually reload the external
    published outputs and return status completed/success/passed. A dry run uses
    the same gates. Failures/running stages are never deleted.
    """
    root, owner = _scope(scope_root)
    _assert_lock(root, owner)
    stage = _stage_name(stage)
    if (not isinstance(native_status, dict) or
            str(native_status.get("state", native_status.get("status", ""))).upper() != "COMPLETED" or
            isinstance(native_status.get("exit_code"), bool) or native_status.get("exit_code") not in (0, "0:0")):
        raise ValueError("Cleanup requires explicit completed native execution with exit code zero")
    registration = _read(root / f".katana-stage-{stage}.json")
    work = _bounded_path(root / "work" / stage)
    if (registration.get("scope_sha256") != owner["record_sha256"] or
            registration.get("work_root") != str(work) or registration.get("stage") != stage):
        raise ValueError("Stage is not owned by this managed run")
    receipt_path = root / f"cleanup-{stage}.receipt.json"
    if receipt_path.exists():
        prior = _read(receipt_path)
        if work.exists() or prior.get("scope_sha256") != owner["record_sha256"] or prior.get("status") != "completed":
            raise ValueError("Previously cleaned stage has unexpected new data or changed ownership")
        return prior
    pending = root / f"cleanup-{stage}.pending.json"
    if pending.exists():
        raise RuntimeError("An earlier cleanup has an unresolved outcome; inspect its pending record, do not retry blindly")
    if not work.is_dir() or list(_node(work.lstat())) != registration["work_identity"]:
        raise ValueError("Registered work directory is missing or changed identity")
    if not callable(checkpoint_validator):
        raise TypeError("A real independent-checkpoint validation callback is required")
    checkpoint = checkpoint_validator(work)
    if not isinstance(checkpoint, dict) or checkpoint.get("status") not in {"completed", "success", "passed"}:
        raise ValueError("Independent retained outputs have not passed checkpoint validation")
    root2, owner2 = _scope(root)
    _assert_lock(root2, owner2)
    if owner2 != owner or list(_node(work.lstat())) != registration["work_identity"]:
        raise RuntimeError("Managed scope or stage changed during checkpoint validation")
    rows = _inventory(work)
    accounting = _sizes(rows)
    plan = {"schema": CLEANUP_SCHEMA, "status": "planned", "stage": stage,
            "scope_sha256": owner["record_sha256"], "run_id": owner["run_id"],
            "work_root": str(work), "native_status": native_status, "checkpoint": checkpoint,
            "accounting": accounting, "inventory_sha256": _hash(rows), "created_at": _now(),
            "resume_lost": f"Nextflow -resume from deleted {stage} work is unavailable; use the retained verified checkpoint.",
            "scope": "Only the registered stage work; input/results/shared environments are excluded"}
    if dry_run:
        return _signed({**plan, "dry_run": True, "status": "planned"})
    _verify_inventory(work, rows)
    quarantine_parent = _bounded_path(root / ".cleanup-quarantine")
    quarantine_parent.mkdir(exist_ok=True)
    quarantine = quarantine_parent / f"{stage}-{uuid4()}"
    plan["quarantine"] = str(quarantine)
    _assert_lock(root, owner)
    _write_new(pending, _signed(plan))
    # Same filesystem/scope rename makes this phase's deletion independent of
    # old path references. The run lock prevents a second managed writer.
    try:
        _assert_lock(root, owner)
        work.rename(quarantine)
        if list(_node(quarantine.lstat())) != registration["work_identity"]:
            raise RuntimeError("Work identity changed at quarantine handoff")
        _verify_inventory(quarantine, rows)
        _assert_lock(root, owner)
        _remove_verified_tree(quarantine, rows)
        receipt = _signed({**plan, "status": "completed", "finished_at": _now(),
                           "deleted_entries": len(rows), "dry_run": False})
        _write_new(receipt_path, receipt)
        pending.unlink()
        if not any(quarantine_parent.iterdir()):
            quarantine_parent.rmdir()
        return receipt
    except BaseException as error:
        # Preserve pending ownership/evidence for explicit inspection. Never
        # erase a claim or repeat an uncertain deletion automatically.
        failure = root / f"cleanup-{stage}.failure-{uuid4()}.json"
        _write_new(failure, _signed({**plan, "status": "failed", "finished_at": _now(),
                                    "error": f"{type(error).__name__}: {error}"}))
        raise


__all__ = ["create_managed_scope", "managed_run_lock", "register_stage", "cleanup_stage"]
