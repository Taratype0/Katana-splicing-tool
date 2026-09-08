from __future__ import annotations

import json
import threading
from pathlib import Path

import duckdb

from katana_next.domain.models import AnalysisArtifact, ProjectSnapshot, SelectionSet


class ArtifactCatalog:
    """Project-local metadata catalog backed by DuckDB."""

    def __init__(self, output_root: str | Path) -> None:
        self.output_root = Path(output_root).resolve()
        self.catalog_dir = self.output_root / ".katana"
        self.catalog_dir.mkdir(parents=True, exist_ok=True)
        self.path = self.catalog_dir / "catalog.duckdb"
        self._lock = threading.RLock()
        self._initialize()

    def _connect(self) -> duckdb.DuckDBPyConnection:
        return duckdb.connect(str(self.path))

    def _initialize(self) -> None:
        with self._lock, self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS projects (
                    project_id VARCHAR PRIMARY KEY,
                    manifest_path VARCHAR NOT NULL,
                    manifest_fingerprint VARCHAR NOT NULL,
                    snapshot_json VARCHAR NOT NULL,
                    opened_at TIMESTAMP NOT NULL
                )
                """
            )
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS artifacts (
                    artifact_id VARCHAR PRIMARY KEY,
                    project_id VARCHAR NOT NULL,
                    analysis_id VARCHAR NOT NULL,
                    kind VARCHAR NOT NULL,
                    entity_type VARCHAR NOT NULL,
                    status VARCHAR NOT NULL,
                    fingerprint VARCHAR NOT NULL,
                    run_fingerprint VARCHAR NOT NULL,
                    path VARCHAR NOT NULL,
                    created_at TIMESTAMP NOT NULL,
                    record_json VARCHAR NOT NULL
                )
                """
            )
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS selections (
                    selection_id VARCHAR PRIMARY KEY,
                    project_id VARCHAR NOT NULL,
                    source_artifact_id VARCHAR NOT NULL,
                    entity_type VARCHAR NOT NULL,
                    name VARCHAR NOT NULL,
                    created_at TIMESTAMP NOT NULL,
                    record_json VARCHAR NOT NULL
                )
                """
            )
            connection.execute(
                "CREATE INDEX IF NOT EXISTS artifact_run_idx "
                "ON artifacts(project_id, analysis_id, run_fingerprint, status)"
            )

    def upsert_project(self, snapshot: ProjectSnapshot) -> None:
        payload = snapshot.model_dump_json()
        with self._lock, self._connect() as connection:
            connection.execute("BEGIN TRANSACTION")
            try:
                connection.execute(
                    "DELETE FROM projects WHERE project_id = ?", [snapshot.project_id]
                )
                connection.execute(
                    "INSERT INTO projects VALUES (?, ?, ?, ?, ?)",
                    [
                        snapshot.project_id,
                        str(snapshot.manifest_path),
                        snapshot.manifest_fingerprint,
                        payload,
                        snapshot.opened_at,
                    ],
                )
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise

    def get_project(self, project_id: str) -> ProjectSnapshot | None:
        with self._lock, self._connect() as connection:
            row = connection.execute(
                "SELECT snapshot_json FROM projects WHERE project_id = ?", [project_id]
            ).fetchone()
        return ProjectSnapshot.model_validate_json(row[0]) if row else None

    def register_artifact(self, artifact: AnalysisArtifact) -> AnalysisArtifact:
        payload = artifact.model_dump_json()
        with self._lock, self._connect() as connection:
            connection.execute("BEGIN TRANSACTION")
            try:
                connection.execute(
                    """
                    DELETE FROM artifacts
                    WHERE project_id = ? AND analysis_id = ?
                      AND run_fingerprint = ? AND fingerprint = ?
                    """,
                    [
                        artifact.project_id,
                        artifact.analysis_id,
                        artifact.run_fingerprint,
                        artifact.fingerprint,
                    ],
                )
                connection.execute(
                    "DELETE FROM artifacts WHERE artifact_id = ?",
                    [artifact.artifact_id],
                )
                connection.execute(
                    "INSERT INTO artifacts VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [
                        artifact.artifact_id,
                        artifact.project_id,
                        artifact.analysis_id,
                        artifact.kind.value,
                        artifact.entity_type,
                        artifact.status.value,
                        artifact.fingerprint,
                        artifact.run_fingerprint,
                        str(artifact.path),
                        artifact.created_at,
                        payload,
                    ],
                )
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        return artifact

    def get_artifact(self, artifact_id: str) -> AnalysisArtifact | None:
        with self._lock, self._connect() as connection:
            row = connection.execute(
                "SELECT record_json FROM artifacts WHERE artifact_id = ?", [artifact_id]
            ).fetchone()
        return AnalysisArtifact.model_validate_json(row[0]) if row else None

    def list_artifacts(
        self,
        project_id: str,
        *,
        analysis_id: str | None = None,
        kind: str | None = None,
    ) -> list[AnalysisArtifact]:
        query = "SELECT record_json FROM artifacts WHERE project_id = ?"
        parameters: list[str] = [project_id]
        if analysis_id:
            query += " AND analysis_id = ?"
            parameters.append(analysis_id)
        if kind:
            query += " AND kind = ?"
            parameters.append(kind)
        query += " ORDER BY created_at DESC, artifact_id"
        with self._lock, self._connect() as connection:
            rows = connection.execute(query, parameters).fetchall()
        return [AnalysisArtifact.model_validate_json(row[0]) for row in rows]

    def cached_artifacts(
        self, project_id: str, analysis_id: str, run_fingerprint: str
    ) -> list[AnalysisArtifact]:
        with self._lock, self._connect() as connection:
            rows = connection.execute(
                """
                SELECT record_json
                FROM artifacts
                WHERE project_id = ? AND analysis_id = ?
                  AND run_fingerprint = ? AND status = 'complete'
                ORDER BY artifact_id
                """,
                [project_id, analysis_id, run_fingerprint],
            ).fetchall()
        artifacts = [AnalysisArtifact.model_validate_json(row[0]) for row in rows]
        return [artifact for artifact in artifacts if artifact.path.exists()]

    def save_selection(self, selection: SelectionSet) -> SelectionSet:
        with self._lock, self._connect() as connection:
            connection.execute("BEGIN TRANSACTION")
            try:
                connection.execute(
                    "DELETE FROM selections WHERE selection_id = ?",
                    [selection.selection_id],
                )
                connection.execute(
                    "INSERT INTO selections VALUES (?, ?, ?, ?, ?, ?, ?)",
                    [
                        selection.selection_id,
                        selection.project_id,
                        selection.source_artifact_id,
                        selection.entity_type,
                        selection.name,
                        selection.created_at,
                        selection.model_dump_json(),
                    ],
                )
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        return selection

    def list_selections(self, project_id: str) -> list[SelectionSet]:
        with self._lock, self._connect() as connection:
            rows = connection.execute(
                "SELECT record_json FROM selections WHERE project_id = ? "
                "ORDER BY created_at DESC",
                [project_id],
            ).fetchall()
        return [SelectionSet.model_validate_json(row[0]) for row in rows]

    def export_manifest(self, project_id: str) -> Path:
        snapshot = self.get_project(project_id)
        if snapshot is None:
            raise KeyError(project_id)
        artifacts = self.list_artifacts(project_id)
        selections = self.list_selections(project_id)
        target = self.catalog_dir / "artifact_manifest.json"
        target.write_text(
            json.dumps(
                {
                    "schema_version": "2.0",
                    "project": snapshot.model_dump(mode="json"),
                    "artifacts": [item.model_dump(mode="json") for item in artifacts],
                    "selections": [item.model_dump(mode="json") for item in selections],
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        return target
