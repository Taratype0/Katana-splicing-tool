"""Verify lossless native-table projections without changing existing projects."""
from __future__ import annotations

import csv
import hashlib
import io
from pathlib import Path, PurePosixPath


def dexseq_projection_identity(asset, input_root):
    """Return native path/hash and TSV hash for a recorded CSV projection.

    Older imported projects deliberately retain the upstream CSV hash as
    provenance. It cannot be used as the expected hash of the generated TSV.
    Reproduce only the lossless serialization, never a statistical calculation.
    """
    metadata = asset.metadata
    projection = metadata.get("native_table_projection")
    if not projection or asset.kind != "dexseq_result":
        return None
    name = projection.get("native_path", "")
    relative = PurePosixPath(name)
    root = Path(input_root).resolve()
    if (not name or relative.is_absolute() or ".." in relative.parts or
            "\\" in name or ":" in name or relative.as_posix() != name):
        raise ValueError("Invalid native projection source path")
    native = root.joinpath(*relative.parts).resolve(strict=True)
    if not native.is_relative_to(root) or not native.is_file():
        raise ValueError("Native projection source escapes project inputs")
    if native == asset.path.resolve():
        return None
    if (projection.get("native_delimiter") != "," or projection.get("unit") != "exon_bin"
            or projection.get("native_statistics_recomputed") is not False):
        raise ValueError("Unsupported native table projection")
    content = native.read_bytes()
    native_sha = hashlib.sha256(content).hexdigest()
    if native_sha != projection.get("native_sha256") or native_sha != metadata.get("upstream_sha256"):
        raise ValueError("Native projection source SHA256 mismatch")
    # The original importer uses Path.read_text(), including universal newline
    # translation inside quoted fields. Reproduce that exact serialization.
    rows = list(csv.reader(io.StringIO(content.decode("utf-8-sig"), newline=None), delimiter=","))
    if (not rows or not {"groupID", "featureID", "padj"}.issubset(rows[0]) or
            len(rows[0]) != len(set(rows[0])) or any(len(row) != len(rows[0]) for row in rows[1:])):
        raise ValueError("Invalid native DEXSeq projection source")
    buffer = io.StringIO(newline="")
    csv.writer(buffer, delimiter="\t", lineterminator="\n").writerows(rows)
    projected_sha = hashlib.sha256(buffer.getvalue().encode("utf-8")).hexdigest()
    if metadata.get("sha256") and metadata["sha256"] != projected_sha:
        raise ValueError("Declared projection SHA256 differs from native serialization")
    return native, native_sha, projected_sha
