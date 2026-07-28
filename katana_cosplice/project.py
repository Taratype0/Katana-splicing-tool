from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from katana_cosplice.annotations import AnnotationRegistry
from katana_cosplice.config import CoSpliceConfig, load_config
from katana_cosplice.pipeline import CoSplicePipeline, Progress
from katana_cosplice.stage_registry import stage_specs


@dataclass
class CoSpliceProject:
    """Public facade shared by standalone Python and host-application adapters."""

    config: CoSpliceConfig
    annotations: AnnotationRegistry = field(default_factory=AnnotationRegistry)

    @classmethod
    def open(cls, config_path: str | Path) -> "CoSpliceProject":
        return cls(load_config(Path(config_path)))

    def pipeline(self, progress: Progress = None) -> CoSplicePipeline:
        return CoSplicePipeline(self.config, progress=progress)

    def stages(self) -> list[dict[str, object]]:
        return [
            {
                "key": spec.key,
                "title": spec.title,
                "description": spec.description,
                "dependencies": list(spec.dependencies),
                "optional_dependencies": list(spec.optional_dependencies),
                "explicit_only": spec.explicit_only,
            }
            for spec in stage_specs(visible_only=True)
        ]

    def state(self) -> dict[str, object]:
        path = self.config.output_dir / "state.json"
        if not path.exists():
            return {"project_name": self.config.project_name, "stages": {}}
        return json.loads(path.read_text(encoding="utf-8"))

    def artifact_manifest(self) -> dict[str, object]:
        path = self.config.output_dir / "artifact_manifest.json"
        if not path.exists():
            return {"schema_version": 1, "project_name": self.config.project_name, "stages": {}}
        return json.loads(path.read_text(encoding="utf-8"))
