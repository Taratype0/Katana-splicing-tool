from __future__ import annotations

import json
from pathlib import Path

from katana_next.domain.models import ProjectManifest


class ManifestStore:
    def load(self, path: str | Path) -> tuple[Path, ProjectManifest]:
        manifest_path = Path(path).expanduser().resolve()
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Project manifest not found: {manifest_path}")
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest = ProjectManifest.model_validate(payload)
        return manifest_path, self._resolve_paths(manifest, manifest_path.parent)

    def save(self, manifest: ProjectManifest, path: str | Path) -> Path:
        target = Path(path).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(manifest.model_dump(mode="json"), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        return target

    @staticmethod
    def _resolve_paths(manifest: ProjectManifest, base: Path) -> ProjectManifest:
        payload = manifest.model_dump()
        input_root = Path(payload["input_root"])
        output_root = Path(payload["output_root"])
        payload["input_root"] = (
            (base / input_root).resolve()
            if not input_root.is_absolute()
            else input_root.resolve()
        )
        payload["output_root"] = (
            (base / output_root).resolve()
            if not output_root.is_absolute()
            else output_root.resolve()
        )

        annotation = dict(payload["annotation"])
        for key in ("gtf", "genome_fasta", "transcript_fasta"):
            value = annotation.get(key)
            if value is not None:
                value_path = Path(value)
                annotation[key] = (
                    payload["input_root"] / value_path
                    if not value_path.is_absolute()
                    else value_path
                ).resolve()
        payload["annotation"] = annotation

        assets = []
        for asset in payload["assets"]:
            item = dict(asset)
            value_path = Path(item["path"])
            item["path"] = (
                payload["input_root"] / value_path
                if not value_path.is_absolute()
                else value_path
            ).resolve()
            assets.append(item)
        payload["assets"] = assets

        samples = []
        for sample in payload["samples"]:
            item = dict(sample)
            item["assets"] = {
                key: (
                    payload["input_root"] / Path(value)
                    if not Path(value).is_absolute()
                    else Path(value)
                ).resolve()
                for key, value in item["assets"].items()
            }
            samples.append(item)
        payload["samples"] = samples
        return ProjectManifest.model_validate(payload)
