from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from katana_cosplice.architecture import ArtifactKind
from katana_cosplice.stage_registry import stage_spec
from katana_cosplice.utils import write_json


def _infer_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".png", ".jpg", ".jpeg", ".svg", ".pdf"}:
        return ArtifactKind.PLOT.value
    if suffix in {".tsv", ".csv", ".xlsx"}:
        return ArtifactKind.TABLE.value
    if suffix == ".json":
        return ArtifactKind.PARAMETERS.value
    if suffix in {".html", ".htm"}:
        return ArtifactKind.REPORT.value
    if suffix in {".log", ".txt"}:
        return ArtifactKind.LOG.value
    return ArtifactKind.OTHER.value


class ArtifactManifestStore:
    """Persistent inventory of stage contracts and generated artifacts."""

    schema_version = 1

    def __init__(self, output_root: Path, project_name: str) -> None:
        self.output_root = output_root.resolve()
        self.project_name = project_name
        self.path = self.output_root / "artifact_manifest.json"

    def read(self) -> dict[str, object]:
        if not self.path.exists():
            return {
                "schema_version": self.schema_version,
                "project_name": self.project_name,
                "stages": {},
            }
        return json.loads(self.path.read_text(encoding="utf-8"))

    def update_stage(
        self,
        stage: str,
        status: str,
        message: str,
        stage_fingerprint: str,
        outputs: list[Path] | None,
    ) -> None:
        manifest = self.read()
        stages = manifest.setdefault("stages", {})
        previous = stages.get(stage, {})
        generated = previous.get("generated_artifacts", [])
        if outputs is not None:
            generated = [self._record(path, stage) for path in outputs]
        spec = stage_spec(stage)
        stages[stage] = {
            "title": spec.title,
            "description": spec.description,
            "status": status,
            "message": message,
            "fingerprint": stage_fingerprint,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "dependencies": list(spec.dependencies),
            "optional_dependencies": list(spec.optional_dependencies),
            "explicit_only": spec.explicit_only,
            "expected_artifacts": [
                {
                    "key": artifact.key,
                    "path": artifact.relative_path,
                    "kind": artifact.kind.value,
                    "required": artifact.required,
                }
                for artifact in spec.artifacts
            ],
            "generated_artifacts": generated,
        }
        write_json(self.path, manifest)

    def _record(self, path: Path, stage: str) -> dict[str, object]:
        resolved = path.resolve()
        try:
            display_path = str(resolved.relative_to(self.output_root))
        except ValueError:
            display_path = str(resolved)
        normalized = display_path.replace("\\", "/")
        declared = {
            artifact.relative_path.replace("\\", "/"): artifact.kind.value
            for artifact in stage_spec(stage).artifacts
        }
        exists = resolved.exists()
        stat = resolved.stat() if exists else None
        return {
            "path": display_path,
            "kind": declared.get(normalized, _infer_kind(resolved)),
            "exists": exists,
            "size": stat.st_size if stat else 0,
            "modified_at": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat() if stat else "",
        }
