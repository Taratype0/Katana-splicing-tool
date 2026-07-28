from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


ARTIFACT_REGISTRY_SCHEMA_VERSION = 1
_SAFE_TOKEN = re.compile(r"[^A-Za-z0-9_.-]+")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_token(value: str) -> str:
    token = _SAFE_TOKEN.sub("_", str(value or "").strip()).strip("._")
    return token or "project"


def _stable_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"), default=str)


def artifact_fingerprint(value: object) -> str:
    return hashlib.sha256(_stable_json(value).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class AnalysisArtifact:
    artifact_id: str
    module_id: str
    artifact_key: str
    kind: str
    entity_type: str
    scope: str = "project"
    scope_id: str = "project"
    comparison_ids: tuple[str, ...] = ()
    path: str = ""
    status: str = "finished"
    fingerprint: str = ""
    source_artifact_ids: tuple[str, ...] = ()
    consumer_modules: tuple[str, ...] = ()
    row_count: int | None = None
    description: str = ""
    selection_rule: str = ""
    parameters: dict[str, object] = field(default_factory=dict)
    created_at: str = ""
    updated_at: str = ""

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["comparison_ids"] = list(self.comparison_ids)
        payload["source_artifact_ids"] = list(self.source_artifact_ids)
        payload["consumer_modules"] = list(self.consumer_modules)
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "AnalysisArtifact":
        return cls(
            artifact_id=str(payload.get("artifact_id", "")),
            module_id=str(payload.get("module_id", "")),
            artifact_key=str(payload.get("artifact_key", "")),
            kind=str(payload.get("kind", "other")),
            entity_type=str(payload.get("entity_type", "other")),
            scope=str(payload.get("scope", "project")),
            scope_id=str(payload.get("scope_id", "project")),
            comparison_ids=tuple(str(item) for item in payload.get("comparison_ids", []) or []),
            path=str(payload.get("path", "")),
            status=str(payload.get("status", "finished")),
            fingerprint=str(payload.get("fingerprint", "")),
            source_artifact_ids=tuple(str(item) for item in payload.get("source_artifact_ids", []) or []),
            consumer_modules=tuple(str(item) for item in payload.get("consumer_modules", []) or []),
            row_count=int(payload["row_count"]) if payload.get("row_count") is not None else None,
            description=str(payload.get("description", "")),
            selection_rule=str(payload.get("selection_rule", "")),
            parameters=dict(payload.get("parameters", {}) or {}),
            created_at=str(payload.get("created_at", "")),
            updated_at=str(payload.get("updated_at", "")),
        )


class AnalysisArtifactRegistry:
    """Project-wide artifact and selection handoff registry.

    The registry stores only lightweight metadata plus explicit selection TSVs.
    Large analysis outputs stay in their module-owned output directories.
    """

    def __init__(self, output_root: str | Path, project_name: str) -> None:
        self.output_root = Path(output_root).resolve()
        self.project_name = str(project_name or "Katana project")
        self.registry_dir = self.output_root / "_artifacts"
        self.path = self.registry_dir / "analysis_artifact_registry.json"

    def read(self) -> dict[str, object]:
        if not self.path.exists():
            return {
                "schema_version": ARTIFACT_REGISTRY_SCHEMA_VERSION,
                "project_name": self.project_name,
                "updated_at": "",
                "artifacts": {},
            }
        try:
            payload = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {
                "schema_version": ARTIFACT_REGISTRY_SCHEMA_VERSION,
                "project_name": self.project_name,
                "updated_at": "",
                "artifacts": {},
            }
        payload.setdefault("schema_version", ARTIFACT_REGISTRY_SCHEMA_VERSION)
        payload.setdefault("project_name", self.project_name)
        payload.setdefault("artifacts", {})
        return payload

    def register(
        self,
        *,
        module_id: str,
        artifact_key: str,
        kind: str,
        entity_type: str,
        path: str | Path | None = None,
        scope: str = "project",
        scope_id: str = "project",
        comparison_ids: Iterable[str] = (),
        status: str = "finished",
        fingerprint: str = "",
        source_artifact_ids: Iterable[str] = (),
        consumer_modules: Iterable[str] = (),
        row_count: int | None = None,
        description: str = "",
        selection_rule: str = "",
        parameters: dict[str, object] | None = None,
    ) -> AnalysisArtifact:
        module_id = str(module_id).strip()
        artifact_key = str(artifact_key).strip()
        scope_id = str(scope_id or "project").strip()
        if not module_id or not artifact_key:
            raise ValueError("module_id and artifact_key are required.")
        artifact_id = f"{module_id}:{scope}:{scope_id}:{artifact_key}"
        now = _utc_now()
        payload = self.read()
        previous = dict((payload.get("artifacts", {}) or {}).get(artifact_id, {}) or {})
        record = AnalysisArtifact(
            artifact_id=artifact_id,
            module_id=module_id,
            artifact_key=artifact_key,
            kind=str(kind or "other"),
            entity_type=str(entity_type or "other"),
            scope=str(scope or "project"),
            scope_id=scope_id,
            comparison_ids=tuple(dict.fromkeys(str(item) for item in comparison_ids if str(item).strip())),
            path=self._display_path(path),
            status=str(status or "finished"),
            fingerprint=str(fingerprint or ""),
            source_artifact_ids=tuple(dict.fromkeys(str(item) for item in source_artifact_ids if str(item).strip())),
            consumer_modules=tuple(dict.fromkeys(str(item) for item in consumer_modules if str(item).strip())),
            row_count=row_count,
            description=str(description or ""),
            selection_rule=str(selection_rule or ""),
            parameters=dict(parameters or {}),
            created_at=str(previous.get("created_at") or now),
            updated_at=now,
        )
        artifacts = payload.setdefault("artifacts", {})
        artifacts[artifact_id] = record.to_dict()
        payload["updated_at"] = now
        self._write(payload)
        return record

    def publish_selection(
        self,
        *,
        module_id: str,
        artifact_key: str,
        entity_type: str,
        frame: pd.DataFrame,
        scope: str,
        scope_id: str,
        comparison_ids: Iterable[str] = (),
        id_columns: Iterable[str] = (),
        source_artifact_ids: Iterable[str] = (),
        consumer_modules: Iterable[str] = (),
        selection_rule: str = "",
        parameters: dict[str, object] | None = None,
    ) -> AnalysisArtifact:
        selection = frame.copy() if frame is not None else pd.DataFrame()
        dedupe_columns = [column for column in id_columns if column in selection.columns]
        if dedupe_columns and not selection.empty:
            selection = selection.drop_duplicates(dedupe_columns, keep="first").reset_index(drop=True)
        selection_dir = (
            self.registry_dir
            / "selections"
            / _safe_token(module_id)
            / _safe_token(scope)
            / _safe_token(scope_id)
        )
        selection_dir.mkdir(parents=True, exist_ok=True)
        path = selection_dir / f"{_safe_token(artifact_key)}.tsv"
        selection.to_csv(path, sep="\t", index=False)
        fingerprint_payload = {
            "module_id": module_id,
            "artifact_key": artifact_key,
            "scope": scope,
            "scope_id": scope_id,
            "columns": list(selection.columns),
            "rows": selection.fillna("").astype(str).to_dict("records"),
            "parameters": parameters or {},
        }
        return self.register(
            module_id=module_id,
            artifact_key=artifact_key,
            kind="selection",
            entity_type=entity_type,
            path=path,
            scope=scope,
            scope_id=scope_id,
            comparison_ids=comparison_ids,
            fingerprint=artifact_fingerprint(fingerprint_payload),
            source_artifact_ids=source_artifact_ids,
            consumer_modules=consumer_modules,
            row_count=int(len(selection)),
            description=f"Explicit {entity_type} selection published for downstream modules.",
            selection_rule=selection_rule,
            parameters=parameters,
        )

    def get(self, artifact_id: str) -> AnalysisArtifact | None:
        raw = (self.read().get("artifacts", {}) or {}).get(str(artifact_id))
        return AnalysisArtifact.from_dict(raw) if isinstance(raw, dict) else None

    def find(
        self,
        *,
        module_id: str | None = None,
        artifact_key: str | None = None,
        kind: str | None = None,
        entity_type: str | None = None,
        scope: str | None = None,
        scope_id: str | None = None,
        comparison_id: str | None = None,
        status: str | None = None,
    ) -> list[AnalysisArtifact]:
        records = [
            AnalysisArtifact.from_dict(raw)
            for raw in (self.read().get("artifacts", {}) or {}).values()
            if isinstance(raw, dict)
        ]
        if module_id is not None:
            records = [record for record in records if record.module_id == module_id]
        if artifact_key is not None:
            records = [record for record in records if record.artifact_key == artifact_key]
        if kind is not None:
            records = [record for record in records if record.kind == kind]
        if entity_type is not None:
            records = [record for record in records if record.entity_type == entity_type]
        if scope is not None:
            records = [record for record in records if record.scope == scope]
        if scope_id is not None:
            records = [record for record in records if record.scope_id == scope_id]
        if comparison_id is not None:
            records = [record for record in records if comparison_id in record.comparison_ids]
        if status is not None:
            records = [record for record in records if record.status == status]
        return sorted(records, key=lambda record: (record.updated_at, record.artifact_id), reverse=True)

    def latest_selection(
        self,
        *,
        entity_type: str,
        scope_id: str,
        module_id: str | None = None,
        artifact_key: str | None = None,
    ) -> tuple[pd.DataFrame, AnalysisArtifact | None]:
        records = self.find(
            module_id=module_id,
            artifact_key=artifact_key,
            kind="selection",
            entity_type=entity_type,
            scope_id=scope_id,
            status="finished",
        )
        if not records:
            return pd.DataFrame(), None
        record = records[0]
        path = self.resolve_path(record)
        if not path.exists():
            return pd.DataFrame(), record
        try:
            return pd.read_csv(path, sep="\t"), record
        except (OSError, pd.errors.ParserError, UnicodeDecodeError):
            return pd.DataFrame(), record

    def resolve_path(self, artifact: AnalysisArtifact) -> Path:
        path = Path(artifact.path)
        return path if path.is_absolute() else self.output_root / path

    def _display_path(self, path: str | Path | None) -> str:
        if path is None or str(path).strip() == "":
            return ""
        resolved = Path(path).resolve()
        try:
            return resolved.relative_to(self.output_root).as_posix()
        except ValueError:
            return str(resolved)

    def _write(self, payload: dict[str, object]) -> None:
        self.registry_dir.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True, default=str),
            encoding="utf-8",
        )
        os.replace(temporary, self.path)
