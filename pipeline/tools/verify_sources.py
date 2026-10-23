#!/usr/bin/env python3
"""Verify every vendored upstream source against its snapshot or explicit patch."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    source = json.loads((root / "UPSTREAM_SOURCE.json").read_text(encoding="utf-8"))
    integration = json.loads((root / "INTEGRATION_CHANGES.json").read_text(encoding="utf-8"))
    patches = {item["path"]: item for item in integration["changes"]}
    originals = {item["destination"]: item for item in source["files"]}
    if len(originals) != len(source["files"]) or len(patches) != len(integration["changes"]):
        raise ValueError("Duplicate source or patch identities")
    if set(patches) - set(originals):
        raise ValueError("An integration patch has no copied source identity")
    for relative, original in originals.items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Source path escapes the pipeline")
        expected = original["sha256"]
        if relative in patches:
            if patches[relative]["snapshot_sha256"] != expected:
                raise ValueError("Patch references a different source snapshot")
            expected = patches[relative]["integrated_sha256"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source changed without a recorded integration revision: {relative}")
    print(json.dumps({"status": "passed", "verified_source_files": len(originals),
                      "explicit_integration_changes": len(patches)}))


if __name__ == "__main__":
    main()
