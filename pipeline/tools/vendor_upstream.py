#!/usr/bin/env python3
"""Copy the executable upstream closure without data, work files or Git metadata.

This maintenance command is deliberately separate from scientific execution.
The copied bytes and any uncommitted contrast-splitting changes are identified
in the resulting source ledger. It never modifies the source checkout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


def copy_sources(source: Path, destination: Path) -> dict:
    source = source.resolve(strict=True)
    destination = destination.resolve()
    if (destination / "UPSTREAM_SOURCE.json").exists():
        raise ValueError("A source snapshot already exists; preserve it and create a new destination")
    pending = [Path("workflows/katana_core.nf")]
    closure: set[Path] = set()
    while pending:
        relative = pending.pop()
        if relative in closure:
            continue
        path = source / relative
        closure.add(relative)
        text = path.read_text(encoding="utf-8-sig")
        for included in re.findall(r"from\s+['\"]([^'\"]+)['\"]", text):
            if not included.startswith("."):
                raise ValueError(f"Unpinned external include in {relative}: {included}")
            candidate = (path.parent / included).resolve()
            candidate = candidate / "main.nf" if candidate.is_dir() else candidate.with_suffix(".nf")
            pending.append(candidate.relative_to(source))
        # Preserve module environment specifications and attribution.
        for name in ("environment.yml", "meta.yml"):
            sibling = path.parent / name
            if sibling.is_file():
                closure.add(sibling.relative_to(source))
        for template in re.findall(r"template\s+['\"]([^'\"]+)['\"]", text):
            closure.add((path.parent / "templates" / template).relative_to(source))
    closure.update(path.relative_to(source) for path in (source / "bin").iterdir() if path.is_file())
    closure.update(map(Path, (
        "assets/software-provenance.lock.json", "assets/migration-source-manifest.json",
        "conf/base.config", "conf/modules.config", "modules.json", "LICENSE", "CITATIONS.md",
    )))
    files = []
    for relative in sorted(closure):
        original = source / relative
        if not original.is_file():
            raise ValueError(f"Missing upstream source: {relative}")
        # Keep root config/license names distinct from the unified entry point.
        output = relative
        if relative.parts[0] == "conf":
            output = Path("conf") / ("upstream_" + relative.name)
        elif relative == Path("LICENSE"):
            output = Path("UPSTREAM_LICENSE")
        elif relative == Path("CITATIONS.md"):
            output = Path("UPSTREAM_CITATIONS.md")
        target = destination / output
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, target)
        if original.parent == source / "bin":
            target.chmod(0o755)
        files.append({"source": relative.as_posix(), "destination": output.as_posix(),
                      "sha256": hashlib.sha256(original.read_bytes()).hexdigest()})
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=source, text=True)
    manifest = {
        "schema": "katana.upstream_source_snapshot/v1", "base_commit": revision,
        "working_tree_modifications_included": bool(status.strip()),
        "working_tree_status_sha256": hashlib.sha256(status.encode()).hexdigest(),
        "source_kind": "existing_bulk_upstream_with_per_contrast_differential_tasks",
        "official_nfcore_release": False,
        "data_or_work_directories_copied": False, "files": files,
    }
    (destination / "UPSTREAM_SOURCE.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    result = copy_sources(args.source, args.destination)
    print(json.dumps({"copied_files": len(result["files"]), "base_commit": result["base_commit"]}))


if __name__ == "__main__":
    main()
