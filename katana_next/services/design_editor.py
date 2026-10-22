"""Explicit, revision-checked metadata edits; imported contrasts remain immutable."""
import json
import os
from pathlib import Path
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field

from katana_next.domain.models import ProjectManifest
from katana_next.infrastructure.fingerprints import fingerprint_file
from katana_next.services.readiness import experimental_design


class DesignRow(BaseModel):
    model_config = ConfigDict(extra="forbid")
    sample_id: str
    biological_sample_id: str
    condition: str
    batch: str | None = None
    attributes: dict[str, str] = Field(default_factory=dict)


class DesignEdit(BaseModel):
    model_config = ConfigDict(extra="forbid")
    revision: str
    samples: list[DesignRow]
    reason: str


def editable_design(context):
    return {"revision": fingerprint_file(context.manifest_path, force_full=True).digest,
            "samples": [{k: s.model_dump(mode="json")[k] for k in DesignRow.model_fields}
                        for s in context.manifest.samples],
            "review": experimental_design(context.manifest),
            "scope": "Edits update declared metadata and invalidate cached analyses. Imported contrast arms and directions are never reversed. Cell-level metadata stays a separately declared input."}


def revise_design(projects, project_id, edit, *, apply=False):
    context = projects.get(project_id)
    edit = DesignEdit.model_validate(edit)
    if edit.revision != fingerprint_file(context.manifest_path, force_full=True).digest:
        raise ValueError("Experimental design changed since opening; reload before saving")
    if not edit.reason.strip():
        raise ValueError("Provide a reason for the design revision")
    old = {s.sample_id: s for s in context.manifest.samples}
    if len(edit.samples) != len(old) or {s.sample_id for s in edit.samples} != set(old):
        raise ValueError("All existing sample IDs must appear exactly once; input identity cannot be edited here")
    payload = context.manifest.model_dump(mode="json")
    for row in edit.samples:
        if row.batch is not None:
            if any(c in row.batch for c in "\t\r\n"):
                raise ValueError("Batch values cannot contain tabs/newlines")
            row.batch = row.batch.strip() or None
        for value in (row.sample_id, row.biological_sample_id, row.condition, *row.attributes.keys(), *row.attributes.values()):
            if not value.strip() or any(c in value for c in "\t\r\n"):
                raise ValueError("Design values must be nonempty and contain no tabs/newlines")
        if set(row.attributes) & {"sample_id", "biological_sample_id", "condition", "batch", "genome_bam", "bam", "quant_sf", "donor_id"}:
            raise ValueError("Factor names cannot override sample identity, condition or batch")
    payload["samples"] = [{**old[s.sample_id].model_dump(mode="json"), **s.model_dump()} for s in edit.samples]
    payload["metadata"]["design_revision"] = {"id": str(uuid4()), "reason": edit.reason.strip(),
                                                "parent_manifest_sha256": edit.revision}
    revised = ProjectManifest.model_validate(payload)
    result = {"applied": False, "review": experimental_design(revised),
              "scope": "Existing imported statistical contrasts retain their sample arms and effect direction; no analysis is launched."}
    if not apply:
        return result
    folder = context.manifest.output_root / "design_history"
    folder.mkdir(parents=True, exist_ok=True)
    backup = folder / f"{edit.revision}.manifest.json"
    if not backup.exists():
        backup.write_bytes(context.manifest_path.read_bytes())
    temporary = context.manifest_path.with_name(f".{context.manifest_path.name}.{uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(revised.model_dump(mode="json"), indent=2), encoding="utf-8")
        os.replace(temporary, context.manifest_path)
    finally:
        temporary.unlink(missing_ok=True)
    new = projects.open_project(context.manifest_path)
    return {**result, "applied": True, "backup": str(backup), "revision": editable_design(new)["revision"],
            "snapshot": new.snapshot.model_dump(mode="json")}


def sync_cosplice_design(config, manifest, output_root):
    """Apply edited metadata by exact sample ID, preserving all engine input paths."""
    if "design_revision" not in manifest.metadata:
        return
    import pandas as pd
    frame = pd.read_csv(config.sample_metadata, sep="\t", dtype=str, keep_default_na=False)
    if "sample_id" not in frame or frame.sample_id.duplicated().any():
        raise ValueError("CoSplice metadata must have unique sample_id to apply an edited design")
    declared = {s.sample_id: s for s in manifest.samples}
    if set(frame.sample_id) != set(declared):
        raise ValueError("Edited project sample IDs do not exactly match CoSplice metadata; explicitly remap before running")
    for index, row in frame.iterrows():
        s = declared[row.sample_id]
        values = {"condition": s.condition, "biological_sample_id": s.biological_sample_id,
                  "batch": s.batch or "", **s.attributes}
        for key, value in values.items():
            frame.loc[index, key] = value
    output_root.mkdir(parents=True, exist_ok=True)
    config.sample_metadata = output_root / "edited_sample_metadata.tsv"
    payload = frame.to_csv(sep="\t", index=False).encode("utf-8")
    # Stage fingerprints include metadata mtime. Reopening an unchanged design
    # or changing only a label view must not invalidate the scientific cache.
    if not config.sample_metadata.is_file() or config.sample_metadata.read_bytes() != payload:
        config.sample_metadata.write_bytes(payload)
