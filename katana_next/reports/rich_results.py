"""Read-only network and MotifTrace report assets from completed workflows.

This adapter does not discover modules, scan RNA, rank motifs, or fit a model.
It preserves every registered figure/table/browser, including the original
measurement population, contrast and statistical family.  A result from another
project is never borrowed to fill a missing section.  Optional search roots are
accepted only through completed ``parameters.json`` workflow manifests.
"""
from __future__ import annotations

from collections import deque
from hashlib import sha256
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Iterable
from urllib.parse import unquote, urlsplit


ADAPTER_VERSION = "rich-results-v1"
NETWORK_WORKFLOWS = frozenset({
    "programs.cosplice", "programs.psi_network", "programs.djexpress",
    "programs.psi_projection", "programs.module_functions", "programs.hub_events",
    "programs.module_explanation", "programs.module_diagnostics",
    "programs.explain_expression", "programs.diagnostics",
    "programs.external_preservation", "programs.external_validation",
    "single_cell.programs",
})
MOTIF_WORKFLOWS = frozenset({
    "prioritize.rbp_support", "prioritize.rbp_mechanism",
    "prioritize.event_motifs", "prioritize.binding_enrichment",
    "prioritize.se_rbp_impact", "prioritize.rbp_sequence_results",
    "prioritize.motiftrace_consequences", "prioritize.motiftrace_usage",
    "prioritize.motiftrace_parallel", "prioritize.mapp_native",
})
_RS_REQUIRED = frozenset({
    "prioritize.rbp_mechanism", "prioritize.event_motifs",
    "prioritize.binding_enrichment", "prioritize.se_rbp_impact",
    "prioritize.motiftrace_consequences",
})
_FIXED_RESULTS = frozenset({
    "prioritize.rbp_sequence_results", "prioritize.motiftrace_usage",
    "prioritize.motiftrace_parallel",
})
_STATUSES = {"complete", "completed", "available"}
_TABLE_SUFFIXES = (".tsv", ".csv", ".parquet", ".tsv.gz", ".csv.gz")
_NON_RESULTS = frozenset({"__pycache__", ".git", ".katana", ".nextflow"})
_API_ARTIFACT = re.compile(r"/api/projects/[^/\s\"'<>]+/artifacts/([^/\s?\"'<>]+)/file(?:\?[^\s\"'<>]*)?")
_CSS_URL = re.compile(r"url\(\s*['\"]?([^)'\"\s]+)", re.I)
_SCRIPT_FILE = re.compile(r"\b(?:fetch|import)\(\s*['\"]([^'\"]+)['\"]\s*\)")


def _chapter(workflow: str) -> int | None:
    if workflow in NETWORK_WORKFLOWS:
        return 6
    if workflow in MOTIF_WORKFLOWS:
        return 7
    return None


def _kind(path: Path) -> str | None:
    name = path.name.lower()
    if path.suffix.lower() in {".png", ".svg", ".jpg", ".jpeg", ".webp"}:
        return "image"
    if path.suffix.lower() in {".html", ".htm"}:
        return "html"
    if name.endswith((".vl.json", ".vega.json", ".vegalite.json")):
        return "vega"
    if name.endswith(_TABLE_SUFFIXES):
        return "table"
    return None


def _enum(value: Any) -> str:
    return str(getattr(value, "value", value) or "").lower()


def _inside(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _path(value: Any, project_root: Path | None = None) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (project_root or Path.cwd()) / path
    return path.resolve()


def _root_for(path: Path) -> Path:
    """The closest completed workflow root is a portability boundary."""
    for parent in (path.parent, *path.parents):
        marker = parent / "parameters.json"
        if marker.is_file():
            try:
                item = json.loads(marker.read_text(encoding="utf-8-sig"))
            except (OSError, ValueError):
                continue
            if item.get("completed") is True and _chapter(str(item.get("workflow_id", ""))):
                return parent
    return path.parent


def _completion(artifact: dict[str, Any], path: Path | None) -> tuple[bool, str]:
    """A registered file does not override an unfinished producer marker."""
    meta = artifact.get("metadata") or {}
    if meta.get("completion_verified") is True:
        return True, "completed_workflow_manifest"
    if path is None:
        return False, "missing_source_path"
    for parent in (path.parent, *path.parents):
        marker = parent / "parameters.json"
        if not marker.is_file():
            continue
        try:
            value = json.loads(marker.read_text(encoding="utf-8-sig"))
        except (ValueError, OSError):
            continue
        if value.get("workflow_id") != artifact.get("analysis_id"):
            continue
        if value.get("completed") is not True:
            return False, "producer_not_completed"
        identifier = artifact.get("artifact_id")
        generated = {row.get("artifact_id") for row in value.get("generated_artifacts", []) if isinstance(row, dict)}
        if identifier and meta.get("artifact_contract_version") == "3.0-correctness" and identifier not in generated:
            return False, "artifact_not_in_completed_producer_manifest"
        return True, "completed_workflow_manifest"
    # Historical/externally registered immutable attachments can live outside
    # their writer run. The caller-selected completed catalogue remains their
    # evidence; do not fabricate a local producer marker or claim it was checked.
    return True, "caller_selected_completed_catalogue"


def _scope(artifact: dict[str, Any]) -> dict[str, Any]:
    workflow = str(artifact.get("analysis_id", ""))
    meta = dict(artifact.get("metadata") or {})
    parameters = dict(artifact.get("parameters") or {})
    recorded = dict(meta.get("das_selection") or {})
    event_set = recorded.get("event_set", parameters.get("event_set"))
    if workflow in NETWORK_WORKFLOWS:
        mode, note = "measurement_population", (
            "Uses the recorded measured-PSI population. The DAS / RS-DAS display "
            "switch does not filter module membership, refit this network or change its P/FDR."
        )
    elif workflow == "prioritize.mapp_native":
        mode, note = "measurement_population", (
            "Official native MAEI fits use their declared measured SE/transcript TPM and motif-count population. "
            "The DAS display switch does not refit them, and native Bonferroni P values are not empirical "
            "MotifTrace/MAPP parallel calibration."
        )
    elif workflow in _RS_REQUIRED:
        mode, note = "rs_required", (
            "RS-DAS is required by this selected-event workflow. This saved collection "
            "retains its recorded read-support scope; an unrecorded historical policy "
            "is not relabelled as the current every-sample ten-read policy."
        )
    elif workflow in _FIXED_RESULTS:
        mode, note = "fixed_completed_result", (
            "Imported completed motif / RNA-usage collection. Its own sample scope, "
            "models, motif sites and statistical families remain unchanged by the DAS switch."
        )
    else:
        mode, note = "switchable_recorded_result", (
            "The saved foreground and statistical family apply only to this collection. "
            "A result for another or unrecorded DAS selection cannot substitute for it."
        )
    return {"mode": mode, "recorded_event_set": event_set,
            "read_support_policy": meta.get("read_support_policy", recorded.get("read_support_policy")),
            "scope_recorded": bool(event_set or mode in {"measurement_population", "fixed_completed_result"}),
            "description": note}


class _References(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.urls: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for key, value in attrs:
            if value and key in {"src", "href", "data", "poster"}:
                self.urls.append(value)
            elif value and key == "srcset":
                self.urls.extend(part.strip().split(" ")[0] for part in value.split(","))


def _references(path: Path) -> list[str]:
    if path.suffix.lower() not in {".html", ".htm", ".css", ".js", ".svg", ".json"}:
        return []
    try:
        text = path.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError):
        return []
    found = _CSS_URL.findall(text) + _SCRIPT_FILE.findall(text)
    if path.suffix.lower() in {".html", ".htm", ".svg"}:
        parser = _References()
        parser.feed(text)
        found.extend(parser.urls)
    if path.suffix.lower() == ".json":
        try:
            value = json.loads(text)
        except ValueError:
            return found
        def urls(item: Any) -> None:
            if isinstance(item, dict):
                for key, child in item.items():
                    if key == "url" and isinstance(child, str):
                        found.append(child)
                    else:
                        urls(child)
            elif isinstance(item, list):
                for child in item:
                    urls(child)
        urls(value)
    return list(dict.fromkeys(found))


def _dependencies(path: Path, root: Path, artifacts_by_id: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Read local references, including app artifact links, without fetching URLs."""
    queue = deque([(path, root, "")])
    files = {path}
    local: dict[Path, str] = {}
    rewrites: list[dict[str, str]] = []
    problems: list[str] = []
    external: set[str] = set()
    while queue:
        current, boundary, destination_prefix = queue.popleft()
        refs = _references(current)
        if current.suffix.lower() in {".html", ".htm", ".js"}:
            try:
                text = current.read_text(encoding="utf-8-sig")
                refs.extend(match.group(0) for match in _API_ARTIFACT.finditer(text))
            except (OSError, UnicodeError):
                pass
        for raw in dict.fromkeys(refs):
            uri = urlsplit(raw)
            if uri.scheme in {"http", "https"} or uri.netloc:
                external.add(raw)
                continue
            if uri.scheme in {"data", "javascript", "mailto", "tel"} or not uri.path:
                continue
            api = _API_ARTIFACT.fullmatch(raw)
            if api:
                artifact = artifacts_by_id.get(api.group(1))
                target = _path(artifact.get("path")) if artifact else None
                if (target is None or not target.is_file() or
                        _enum(artifact.get("status", "complete")) not in _STATUSES or
                        not _completion(artifact, target)[0]):
                    problems.append("Unresolved app artifact link: " + raw)
                    continue
                relative = "_related/" + api.group(1) + "/" + target.name
                local.setdefault(target, relative)
                relative = local[target]
                rewrites.append({"source": str(current), "url": raw, "target": relative})
                next_boundary = target.parent
                next_prefix = str(Path(relative).parent).replace("\\", "/")
            else:
                if uri.scheme or uri.path.startswith(("/", "\\")):
                    problems.append("Non-portable absolute link: " + raw)
                    continue
                target = (current.parent / unquote(uri.path)).resolve()
                if not _inside(target, boundary):
                    problems.append("Local link outside completed run boundary: " + raw)
                    continue
                if not target.is_file():
                    problems.append("Missing local attachment: " + raw)
                    continue
                local.setdefault(target, (Path(destination_prefix) / target.relative_to(boundary)).as_posix())
                next_boundary, next_prefix = boundary, destination_prefix
            if target not in files:
                files.add(target)
                queue.append((target, next_boundary, next_prefix))
    return {"local_dependencies": [{"path": str(p), "relative_path": relative}
                                    for p, relative in sorted(local.items(), key=lambda pair: pair[1])],
            "local_link_rewrites": rewrites, "portability_warnings": sorted(set(problems)),
            "external_urls": sorted(external)}


def _from_completed_manifests(search_roots: Iterable[Path]) -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    for root in dict.fromkeys(Path(p).resolve() for p in search_roots):
        if not root.is_dir():
            continue
        for marker in sorted(root.rglob("parameters.json")):
            if any(part in _NON_RESULTS for part in marker.relative_to(root).parts):
                continue
            try:
                value = json.loads(marker.read_text(encoding="utf-8-sig"))
            except (ValueError, OSError):
                continue
            workflow = str(value.get("workflow_id", ""))
            if value.get("completed") is not True or _chapter(workflow) is None:
                continue
            rows = value.get("generated_artifacts", [])
            if not isinstance(rows, list):
                continue
            for item in rows:
                if not isinstance(item, dict):
                    continue
                raw = _path(item.get("path"), marker.parent)
                # Relocated exports keep their immutable basename/run directory.
                if raw is None or not raw.is_file():
                    relocated = marker.parent / Path(str(item.get("path", ""))).name
                    raw = relocated if relocated.is_file() else raw
                found.append({**item, "path": str(raw) if raw else "", "analysis_id": workflow,
                              "run_fingerprint": value.get("run_fingerprint", ""), "status": "complete",
                              "parameters": value.get("parameters", {}),
                              "metadata": {"workflow_version": value.get("workflow_version"),
                                           "completion_manifest": str(marker),
                                           "completion_verified": True}})
            for exported in value.get("generated_plot_files", []) + value.get("generated_table_files", []):
                export = _path(exported, marker.parent)
                if export is None or not export.is_file():
                    name = Path(str(exported)).name
                    export = marker.parent / ("plots" if name.endswith((".png", ".svg")) else "") / name
                if export.is_file():
                    found.append({"path": str(export), "name": export.stem.replace("_", " "),
                                  "analysis_id": workflow, "status": "complete",
                                  "run_fingerprint": value.get("run_fingerprint", ""),
                                  "parameters": value.get("parameters", {}),
                                  "metadata": {"workflow_version": value.get("workflow_version"),
                                               "completion_manifest": str(marker),
                                               "completion_verified": True}})
    return found


def collect_rich_results(
    artifacts: list[dict[str, Any]], project_root: Path | None = None,
    search_roots: list[Path] | None = None,
) -> list[dict[str, Any]]:
    """Return all existing chapter 6/7 assets and explicit missing-section states.

    ``artifacts`` must already belong to the caller's selected project/snapshot.
    This function accepts dictionaries rather than opening or mutating a project
    catalogue.  ``search_roots`` are explicit completed-result roots, not a broad
    search for similarly named experiments.  Nothing is recomputed or refiltered.
    """
    project_root = Path(project_root).resolve() if project_root else None
    records = [dict(row) for row in artifacts]
    if search_roots:
        records.extend(_from_completed_manifests(search_roots))
    project_ids = {str(row["project_id"]) for row in records if row.get("project_id")}
    if len(project_ids) > 1:
        raise ValueError("Rich report inputs must belong to one explicitly selected project; cohorts cannot be pooled")
    by_id = {str(row["artifact_id"]): row for row in records if row.get("artifact_id")}
    results: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    observed: set[int] = set()
    for artifact in records:
        workflow = str(artifact.get("analysis_id", ""))
        chapter = _chapter(workflow)
        if chapter is None:
            continue
        status = _enum(artifact.get("status", "complete"))
        if status not in _STATUSES:
            continue
        metadata = dict(artifact.get("metadata") or {})
        params = dict(artifact.get("parameters") or {})
        contrasts = artifact.get("contrast_ids") or params.get("contrast_ids") or metadata.get("displayed_comparisons")
        if not contrasts:
            single = metadata.get("comparison_id", metadata.get("contrast_id"))
            contrasts = [single] if single else []
        if isinstance(contrasts, str):
            contrasts = [contrasts]
        contrast_id = str(contrasts[0]) if len(contrasts) == 1 else "all"
        main = _path(artifact.get("path"), project_root)
        completed, completion_evidence = _completion(artifact, main)
        if not completed:
            continue
        paths = []
        if main:
            paths.append(main)
        exports = metadata.get("exports") or {}
        if isinstance(exports, dict):
            paths.extend(p for p in (_path(v, project_root) for v in exports.values()) if p)
        for path in dict.fromkeys(paths):
            kind = _kind(path)
            if kind is None or not path.is_file():
                continue
            scope = _scope(artifact)
            key = (str(path), workflow, str(artifact.get("run_fingerprint", "")), tuple(contrasts),
                   scope["recorded_event_set"])
            if key in seen:
                continue
            seen.add(key)
            observed.add(chapter)
            root = _root_for(path)
            declared_root = _path(metadata.get("portable_root"), project_root)
            if declared_root is not None:
                if not declared_root.is_dir() or not _inside(path, declared_root):
                    raise ValueError("Declared portable root must contain its completed source asset")
                root = declared_root
            deps = _dependencies(path, root, by_id) if kind in {"html", "vega"} else {
                "local_dependencies": [], "local_link_rewrites": [],
                "portability_warnings": [], "external_urls": [],
            }
            stem = path.name
            for ending in (".vl.json", ".tsv.gz", ".csv.gz", path.suffix):
                if stem.endswith(ending):
                    stem = stem[:-len(ending)]
                    break
            title = str(artifact.get("name") or stem.replace("_", " "))
            if path != main and main and _kind(main) != kind:
                title += " — " + path.name
            caption = str(metadata.get("caption") or metadata.get("limitations") or "")
            if not caption:
                caption = ("Completed network output; module/event identities, continuous rankings and "
                           "original association P/FDR are retained." if chapter == 6 else
                           "Completed RBP / MotifTrace output. Identity enrichment, motif matches, CLIP "
                           "and RNA availability retain their separate evidence and correction families.")
            if metadata.get("completion_scope") == "sequence_consequences_only":
                caption += (" Sequence-consequence substep completed; regional motif enrichment and "
                            "whole-comparison inference are not established by this asset.")
            if metadata.get("completion_scope") == "observed_network_baseline_only":
                caption += (" Observed network baseline completed; bootstrap module confidence "
                            "requires its own completed phase and is not established by this asset.")
            stat = path.stat()
            result_meta = {**metadata, "adapter_version": ADAPTER_VERSION,
                           "analysis_id": workflow, "workflow_version": metadata.get("workflow_version"),
                           "source_artifact_id": artifact.get("artifact_id"),
                           "project_id": artifact.get("project_id"),
                           "run_fingerprint": artifact.get("run_fingerprint"),
                           "contrast_ids": list(contrasts), "das_scope": scope,
                           "source_status": status, "new_computation": False,
                           "completion_evidence": completion_evidence,
                           "statistical_values_preserved": True, "all_rows_retained": True,
                           "portable_root": str(root), "relative_path": path.relative_to(root).as_posix(),
                           "source_bytes": stat.st_size, "source_mtime_ns": stat.st_mtime_ns,
                           "display_group_id": str(artifact.get("artifact_id") or path.parent.name) + ":" + stem,
                           "display_in_main": bool(metadata.get("display_in_main", True)) and not (path.suffix.lower() == ".svg" and path.with_suffix(".png").is_file()),
                           **deps}
            if main and main.is_file() and main.name.endswith('.vl.json') and kind in {'vega', 'image'}:
                # This exact export-to-spec relationship comes from one producer
                # artifact, never from scanning similar neighbouring filenames.
                result_meta['saved_plot_specification'] = {
                    'path': str(main), 'source_sha256': sha256(main.read_bytes()).hexdigest(),
                    'export_sha256': sha256(path.read_bytes()).hexdigest(),
                }
                # The specification and its PNG encode the same figure. Keep
                # both in outputs, but show only the PNG in the main report.
                if kind == 'vega' and any(p.suffix.lower() == '.png' and p.is_file() for p in paths):
                    result_meta['display_in_main'] = False
            if metadata.get("completion_scope") == "sequence_consequences_only":
                result_meta["partial_status"] = "sequence_consequences_available_inference_not_completed"
            identifier = sha256(json.dumps(key, ensure_ascii=True).encode()).hexdigest()[:20]
            results.append({"id": "rich-" + identifier, "chapter": chapter,
                            "title": title, "path": str(path), "kind": kind,
                            "contrast_id": contrast_id,
                            "cell_type": str(metadata.get("cell_type", "all")),
                            "module_id": metadata.get("module_id", metadata.get("module")),
                            "caption": caption, "status": "available", "metadata": result_meta,
                            "source_paths": [str(path), *[p["path"] for p in deps["local_dependencies"]]]})
    missing = {
        6: ("Katana-Splice-Network", "No completed network output was supplied for this project. "
            "Module discovery, activity, hubs, enrichment and diagnostics are not computed; "
            "a missing result is not a zero-module network."),
        7: ("RBP and MotifTrace", "No completed RBP / MotifTrace output was supplied for this project. "
            "RBP identity, motif enrichment, sequence logos/windows and RNA consequences "
            "are not computed; results from another cohort are not substituted."),
    }
    for chapter, (title, caption) in missing.items():
        if chapter not in observed:
            results.append({"id": f"rich-chapter-{chapter}-not-computed", "chapter": chapter,
                            "title": title, "path": "", "kind": "table", "contrast_id": "all",
                            "cell_type": "all", "module_id": None, "caption": caption,
                            "status": "not_computed", "source_paths": [],
                            "metadata": {"adapter_version": ADAPTER_VERSION, "new_computation": False,
                                         "reason": "no_completed_project_matched_assets"}})
    # An RBP identity table does not establish that MotifTrace or CLIP ran.
    # Missing optional branches remain explicit even when another branch in the
    # same report chapter has real outputs.
    families = (
        (7, "rbp_identity", "RBP identity / declared gene sets", {"prioritize.rbp_support"},
         "No completed RBP identity or declared-set enrichment result was supplied."),
        (7, "motif_enrichment", "Motif enrichment and matched background",
         {"prioritize.event_motifs", "prioritize.se_rbp_impact", "prioritize.rbp_sequence_results",
          "prioritize.rbp_mechanism", "prioritize.motiftrace_consequences"},
         "No completed motif enrichment result was supplied. Group enrichment FDR, model scores "
         "and per-site statistics cannot substitute for one another."),
        (7, "motiftrace_sequence", "MotifTrace RNA sequences, models and loci",
         {"prioritize.rbp_sequence_results", "prioritize.motiftrace_consequences", "prioritize.motiftrace_usage"},
         "No completed RNA-form / model / locus result was supplied. Sequence logos, windows, "
         "gained/lost/retained loci and RNA-usage availability are not inferred from DAS alone."),
        (7, "clip", "Reference CLIP and binding evidence",
         {"prioritize.binding_enrichment", "prioritize.rbp_mechanism"},
         "No completed matched CLIP result was supplied. Motif matches do not establish binding."),
        (7, "mapp_parallel", "Optional calibrated MAPP / MotifTrace parallel collection", {"prioritize.motiftrace_parallel"},
         "No completed calibrated parallel collection was supplied. Native MAEI fits and MotifTrace computations, "
         "when present, are reported separately; native P values are not relabelled as empirical parallel calibration."),
    )
    for chapter, family, title, workflows, caption in families:
        supporting = [row for row in results if row["status"] == "available"
                      and row["metadata"].get("analysis_id") in workflows]
        if family == "motif_enrichment":
            supporting = [row for row in supporting
                          if row["metadata"].get("completion_scope") != "sequence_consequences_only"]
        if not supporting:
            results.append({"id": "rich-" + family + "-not-computed", "chapter": chapter,
                            "title": title, "path": "", "kind": "table", "contrast_id": "all",
                            "cell_type": "all", "module_id": None, "caption": caption,
                            "status": "not_computed", "source_paths": [],
                            "metadata": {"adapter_version": ADAPTER_VERSION, "component": family,
                                         "new_computation": False,
                                         **({"requirement": "conditional_import", "optional_source_supplied": False}
                                            if family == "mapp_parallel" else {}),
                                         "reason": "no_completed_project_matched_component"}})
    baseline_only = any(row["status"] == "available" and
                        row["metadata"].get("completion_scope") == "observed_network_baseline_only"
                        for row in results)
    bootstrap_complete = any(row["status"] == "available" and
                             row["metadata"].get("completion_scope") in {"bootstrap_confidence", "full_network_with_bootstrap"}
                             for row in results)
    if baseline_only and not bootstrap_complete:
        results.append({"id": "rich-network-bootstrap-not-completed", "chapter": 6,
                        "title": "Bootstrap module confidence", "path": "", "kind": "table",
                        "contrast_id": "all", "cell_type": "all", "module_id": None,
                        "caption": "The observed network baseline is complete. Its independent bootstrap phase is not yet accepted; no module-confidence result is inferred from the baseline.",
                        "status": "not_computed", "source_paths": [],
                        "metadata": {"adapter_version": ADAPTER_VERSION, "component": "network_bootstrap",
                                     "reason": "completed_baseline_is_not_completed_bootstrap"}})
    return sorted(results, key=lambda row: (row["chapter"], row["contrast_id"],
                                            str(row["metadata"].get("analysis_id", "")),
                                            str(row["metadata"].get("run_fingerprint", "")),
                                            row["title"], row["path"]))


def _sha(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def copy_rich_asset(asset: dict[str, Any], destination: Path) -> dict[str, Any]:
    """Copy one asset and local browser attachments to a new portable directory.

    Only destination HTML/CSS/JS links to registered application artifacts may be
    rewritten.  Original files, P/FDR, matrices and images are not altered.  The
    receipt records before/after hashes, so a copied link rewrite is transparent.
    The caller may instead consume ``local_dependencies`` / ``local_link_rewrites``
    with its own report asset writer.
    """
    if asset.get("status") != "available":
        raise ValueError("A not-computed result cannot be copied as a completed asset")
    source = Path(asset["path"]).resolve()
    meta = asset["metadata"]
    portable_root = Path(meta["portable_root"]).resolve()
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError("The portable asset destination already exists; previous output is preserved")
    if _inside(destination, portable_root) or _inside(portable_root, destination):
        raise ValueError("Portable output must be independent of the source run directory")
    files = [{"path": str(source), "relative_path": meta["relative_path"]},
             *meta.get("local_dependencies", [])]
    bounded: dict[Path, Path] = {}
    for row in files:
        path = Path(row["path"]).resolve()
        target = (destination / row["relative_path"]).resolve()
        if not _inside(target, destination) or not path.is_file():
            raise ValueError("Missing or unbounded portable attachment: " + str(row))
        bounded[path] = target
    before = {path: (path.stat().st_size, path.stat().st_mtime_ns, _sha(path)) for path in bounded}
    destination.mkdir(parents=True)
    rewritten: list[dict[str, str]] = []
    for path, target in bounded.items():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        changes = [row for row in meta.get("local_link_rewrites", []) if Path(row["source"]).resolve() == path]
        if changes:
            text = target.read_text(encoding="utf-8-sig")
            for row in changes:
                rel = Path(os.path.relpath(destination / row["target"], target.parent)).as_posix()
                text = text.replace(row["url"], rel)
                rewritten.append({"file": str(target.relative_to(destination)), "from": row["url"], "to": rel})
            target.write_text(text, encoding="utf-8")
    rows = []
    for path, target in bounded.items():
        size, stamp, digest = before[path]
        if (path.stat().st_size, path.stat().st_mtime_ns, _sha(path)) != before[path]:
            raise ValueError("Original result changed during portable copy: " + str(path))
        rows.append({"source": str(path), "destination": target.relative_to(destination).as_posix(),
                     "source_bytes": size, "source_mtime_ns": stamp, "source_sha256": digest,
                     "portable_sha256": _sha(target), "original_preserved": True})
    return {"schema": "katana.rich_asset_copy/v1", "asset_id": asset["id"],
            "path": str(bounded[source]), "files": rows, "link_rewrites": rewritten,
            "portability_warnings": meta.get("portability_warnings", []),
            "originals_preserved": True, "new_scientific_computation": False}
