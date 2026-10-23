"""Relocatable, content-bound author source supplied by the installed package."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path


MAPP_COMMIT = "6ad6e12275e6d7aa1254d1a8933cfbce6766b16c"
MAPP_PROVIDER = "katana.bundled_mapp/v1"
COUNTER = "modules/CREATE_SITECOUNT_MATRICES/scripts/create-kmers-sitecount-matrix.py"
MAEI = (
    "modules/MAEI/scripts/optimize-MARA-model-parameters.R",
    "modules/MAEI/scripts/model.cpp",
    "modules/MAEI/scripts/optim.c",
    "modules/MAEI/scripts/zeroin.c",
    "modules/MAEI/scripts/calculate-Zscores.py",
    "modules/MAEI/scripts/extend-tables.py",
    "modules/MAEI/scripts/calculate_Zscores_statistical_significancy.R",
)


def _source_root():
    return Path(__file__).resolve().parents[1] / "vendor" / "MAPP"


def verified_mapp_sources():
    root = _source_root()
    manifest = json.loads((root / "SOURCE.json").read_text(encoding="utf-8"))
    if (manifest.get("commit") != MAPP_COMMIT or manifest.get("modified") is not False
            or manifest.get("repository") != "https://github.com/gruber-sciencelab/MAPP"):
        raise ValueError("Bundled MAPP author source identity is invalid")
    records = {row["path"]: row for row in manifest["files"]}
    expected = {"LICENSE", COUNTER, *MAEI}
    if set(records) != expected or len(records) != len(manifest["files"]):
        raise ValueError("Bundled MAPP source inventory is incomplete or ambiguous")
    for relative, row in records.items():
        data = (root / relative).read_bytes()
        if len(data) != row["size_bytes"] or hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Bundled MAPP source changed: {relative}")
    descriptor = {"native_source_provider": MAPP_PROVIDER, "native_commit": MAPP_COMMIT,
                  "counter_sha256": records[COUNTER]["sha256"],
                  "source_sha256": {name: records[name]["sha256"] for name in MAEI}}
    return root, descriptor


def mapp_input_defaults(inputs):
    """New plans bind content rather than a machine-specific package directory."""
    values = deepcopy(inputs)
    if values.get("native_root") and not values.get("native_source_provider"):
        return values  # Explicit existing source remains an explicit contract.
    _root, descriptor = verified_mapp_sources()
    for key, value in descriptor.items():
        if key in values and values[key] != value:
            raise ValueError(f"Bundled MAPP {key} differs from the requested author source")
        values[key] = value
    return values


def resolve_bundled_stage(spec):
    """Resolve source files inside this worker without rewriting the frozen plan."""
    inputs = spec.get("inputs", {})
    if not inputs.get("native_source_provider"):
        return spec
    if inputs["native_source_provider"] != MAPP_PROVIDER:
        raise ValueError("Unknown bundled native source provider")
    adapter = spec.get("adapter") or spec.get("workflow_id", "").removeprefix("flow.")
    if adapter not in {"mapp.prepare", "mapp.maei"}:
        raise ValueError("Bundled MAPP source is only valid for MAPP stages")
    root, descriptor = verified_mapp_sources()
    if any(inputs.get(key) != value for key, value in descriptor.items()):
        raise ValueError("Worker MAPP source differs from the plan's content binding")
    resolved = deepcopy(spec)
    resolved["inputs"].update(native_root=str(root),
        native_source_files=[str(root / name) for name in (COUNTER, *MAEI, "SOURCE.json", "LICENSE")])
    return resolved
