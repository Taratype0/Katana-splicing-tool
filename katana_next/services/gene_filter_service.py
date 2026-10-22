from __future__ import annotations

import json
from pathlib import Path
from uuid import NAMESPACE_URL, uuid5

import pandas as pd

from katana_next.domain.models import InputAsset
from katana_next.infrastructure.fingerprints import fingerprint_value, fingerprint_file
from katana_next.science.gene_filter import (filter_adapter, filter_config, gene_filter_details,
                                            gene_filter_footnote)
from src.analysis.gene_filter import project_gene_filter
from src.domain.models.gene_filter import GeneFilterConfig


def configure_gene_filter(projects, context, request):
    """Create/reopen an isolated project scope; never rewrite the source manifest."""
    values = request.model_dump() if hasattr(request, "model_dump") else dict(request)
    annotation = values.pop("annotation_path", "") or ""
    expression = values.pop("expression_matrix_path", "") or ""
    # Direct Python callers and typed API requests use the same new-scope
    # defaults. This does not change legacy manifest decoding.
    config = GeneFilterConfig.from_mapping({**GeneFilterConfig.for_new_analysis().to_dict(), **values})
    original_path = Path(context.manifest.metadata.get("gene_filter_source_manifest", context.manifest_path))
    _, original = projects.manifests.load(original_path)
    if not config.active:
        return projects.open_project(original_path)
    if not any(asset.kind in {"rmats_result", "deg_result", "dtu_result", "dexseq_result"} for asset in original.assets):
        raise ValueError("Project gene exclusion currently applies to bulk DEG/DAS/DTU inputs.")
    settings = {"config": config.to_dict(),
                "application_stage": "post_import_downstream_scope",
                "upstream_model_exclusions": "not_verified",
                "upstream_refit_performed": False}
    if annotation:
        path = Path(annotation).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f"Annotation file is missing: {path}")
        settings["gtf" if path.name.lower().endswith((".gtf", ".gtf.gz")) else "annotation_tsv"] = str(path)
    else:
        previous = context.manifest.metadata.get("gene_filter", {})
        settings.update({key: previous[key] for key in ("annotation_tsv", "gtf") if previous.get(key)})
    if expression:
        path = Path(expression).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f"Expression matrix is missing: {path}")
        settings["expression_matrix"] = str(path)
    elif context.manifest.metadata.get("gene_filter", {}).get("expression_matrix"):
        settings["expression_matrix"] = context.manifest.metadata["gene_filter"]["expression_matrix"]
    candidate = original.model_copy(update={"metadata": {**original.metadata, "gene_filter": settings}})
    details = gene_filter_details(candidate)  # Validate all required annotation/expression evidence before writing.
    source_identity = fingerprint_file(original_path, force_full=True).digest
    scope = fingerprint_value({"source": str(original_path), "source_sha256": source_identity,
                               "filter": details["identity"]})[:20]
    root = original.output_root / "gene_filtered" / scope
    projects._validate_output_root(root)
    manifest_path = root / "katana.project.json"
    assets = list(original.assets)
    adapter = filter_adapter(candidate)
    for key, path in (("annotation", adapter.tool_paths.get("annotation_tsv")),
                      ("gtf", adapter.tool_paths.get("gtf")),
                      ("expression", adapter.counts_path if config.min_expression > 0 else None)):
        if path:
            assets.append(InputAsset(asset_id=f"gene-filter-{key}", kind=f"gene_filter_{key}", path=Path(path)))
    candidate = candidate.model_copy(update={
        "project_id": str(uuid5(NAMESPACE_URL, original.project_id + ":gene-filter:" + scope)),
        "output_root": root / "results", "assets": tuple(assets),
        "metadata": {**candidate.metadata, "gene_filter_source_manifest": str(original_path),
                     "gene_filter_identity": details["identity"]},
    })
    encoded = json.dumps(candidate.model_dump(mode="json"), indent=2, ensure_ascii=False)
    root.mkdir(parents=True, exist_ok=True)
    try:
        with manifest_path.open("x", encoding="utf-8") as handle:
            handle.write(encoded)
    except FileExistsError:
        if manifest_path.read_text(encoding="utf-8") != encoded:
            raise ValueError("Existing filtered scope has different contents; source files or manifest changed.")
    filtered = projects.open_project(manifest_path)
    audit_root = filtered.manifest.output_root / "gene_filter"
    audit_root.mkdir(parents=True, exist_ok=True)
    engine = project_gene_filter(filter_adapter(candidate))
    engine.annotation_audit().to_csv(audit_root / "annotation_decisions.tsv", sep="\t", index=False)
    observations = []
    for asset in candidate.assets:
        if asset.kind not in {"deg_result", "rmats_result", "dtu_result", "dexseq_result"}:
            continue
        frame = pd.read_csv(asset.path, sep="," if asset.path.suffix.lower() == ".csv" else "\t", low_memory=False)
        decisions = engine.decisions(frame)
        observations.append({"asset_id": asset.asset_id, "kind": asset.kind,
                             "contrast_id": asset.metadata.get("contrast_id", ""),
                             "input_rows": len(frame), "excluded_rows": int(decisions.excluded.sum()),
                             "retained_rows": int((~decisions.excluded).sum())})
    pd.DataFrame(observations).to_csv(audit_root / "input_filter_counts.tsv", sep="\t", index=False)
    (audit_root / "gene_filter.json").write_text(json.dumps(details, indent=2, ensure_ascii=False), encoding="utf-8")
    (audit_root / "FOOTNOTE.txt").write_text(gene_filter_footnote(candidate), encoding="utf-8")
    return filtered
