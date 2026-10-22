"""A separate methods/current-view page; frozen engine reports remain references."""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
import re


NETWORK_LABEL_METHODS = (
    "Module membership correlation (kME) and intramodular connectivity follow WGCNA methods "
    "(Langfelder and Horvath, 2008; doi:10.1186/1471-2105-9-559). "
    "Katana's optional co-assignment, rank-fraction, preset-selection-frequency and estimability-coverage cutoffs are operational rules "
    "without independent calibration; they are not universal significance cutoffs. "
    "Continuous measurements and complete within-module rankings remain available whether labels are on or off. "
    "A score, correlation or preset-selection frequency is not a calibrated probability, P value or FDR. "
    "Native test P values and Benjamini-Hochberg FDR retain their own test families; "
    "FDR <=0.05 is the declared significance convention (Benjamini and Hochberg, 1995; "
    "doi:10.1111/j.2517-6161.1995.tb02031.x). Association, bootstrap consistency and external preservation "
    "answer separate questions and do not establish causality."
)


def _save_immutable(path: Path, payload: bytes) -> None:
    if path.is_file():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Refusing to replace an existing network-label report output: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(payload)


def create_network_label_report(raw_report: Path, out: Path, *, policy: dict, table_views: dict[str, Path]) -> Path:
    """Create an explicit companion without relabelling or altering raw HTML."""
    raw_report, out = Path(raw_report), Path(out)
    if raw_report.resolve() == out.resolve():
        raise ValueError("The current-view report must be separate from the frozen raw report")
    raw_bytes = raw_report.read_bytes()
    raw_hash = hashlib.sha256(raw_bytes).hexdigest()
    inventory, rows = {}, []
    used_names = set()
    for key, source in sorted(table_views.items()):
        source = Path(source)
        if not source.is_file():
            raise FileNotFoundError(source)
        stem = re.sub(r"[^A-Za-z0-9_.-]", "_", str(key)) or "view"
        name = stem + (source.suffix or ".tsv")
        if name in used_names:
            raise ValueError("Current network table-view names must be unique")
        used_names.add(name)
        data = source.read_bytes()
        relative = Path("tables") / name
        _save_immutable(out.parent / relative, data)
        inventory[str(key)] = {"file": relative.as_posix(), "sha256": hashlib.sha256(data).hexdigest()}
        rows.append(f"<li><a href='{html.escape(relative.as_posix(), quote=True)}'>{html.escape(str(key))}</a></li>")
    enabled = policy.get("apply_threshold_labels")
    state = "on" if enabled is True else "off" if enabled is False else "not recorded"
    provenance = {"schema": "katana.network_label_report/v1", "network_label_policy": policy,
                  "legacy_report_reference": {"file_name": raw_report.name, "sha256": raw_hash,
                                              "role": "historical/internal reference; not a current label view"},
                  "current_table_views": inventory}
    serialized = json.dumps(provenance, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
    _save_immutable(out.parent / "NETWORK_LABEL_REPORT_PROVENANCE.json", serialized.encode("utf-8"))
    content = ("<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Network descriptive-label views</title>"
               "<style>body{font:16px Arial,sans-serif;max-width:980px;margin:40px auto;padding:20px;line-height:1.5}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style></head><body>"
               "<main><h1>Network descriptive-label views</h1>"
               f"<section><h2>Current label policy</h2><p>Threshold labels: {state}. "
               "Off means labels were not applied, not that members failed a threshold. "
               "Module membership and complete continuous measurements/rankings are retained.</p>"
               f"<pre>{html.escape(json.dumps(policy, indent=2, sort_keys=True, allow_nan=False))}</pre></section>"
               "<section><h2>Current table views</h2><p>These files retain every assessed row and explicit unavailable states.</p><ul>"
               + "".join(rows) + "</ul></section>"
               f"<section><h2>Methods and threshold basis</h2><p>{html.escape(NETWORK_LABEL_METHODS)}</p></section>"
               f"<section><h2>Frozen raw report reference</h2><p>{html.escape(raw_report.name)} is a historical/internal reference. "
               "Its original threshold labels, figures and text have not been converted into current-policy results, "
               "and its body is not presented as a new report here.</p>"
               f"<p>Original SHA256: <code>{raw_hash}</code></p></section>"
               "<p><a href='NETWORK_LABEL_REPORT_PROVENANCE.json'>Current view and source provenance</a></p></main></body></html>")
    _save_immutable(out, content.encode("utf-8"))
    if raw_report.read_bytes() != raw_bytes:
        raise RuntimeError("Frozen raw report changed while building its companion")
    return out
