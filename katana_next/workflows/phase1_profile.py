from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any

from katana_next.domain.models import AnalysisArtifact, ArtifactKind
from katana_next.infrastructure.fingerprints import fingerprint_value
from katana_next.infrastructure.jobs import JobContext
from katana_next.services.project_service import ProjectContext


class Phase1ProfileWorkflow:
    workflow_id = "platform.phase1_profile"
    version = "1.0.0"

    @staticmethod
    def source_asset(context: ProjectContext) -> Path:
        matches = [
            asset.path
            for asset in context.manifest.assets
            if asset.kind == "phase1_event_table"
        ]
        if not matches:
            raise ValueError(
                "The platform fixture workflow requires an input asset with "
                "kind='phase1_event_table'."
            )
        return matches[0]

    def run(
        self,
        context: ProjectContext,
        job: JobContext,
        run_fingerprint: str,
        parameters: dict[str, Any],
    ) -> list[str]:
        source = self.source_asset(context)
        job.report(0.12, "Reading fixture input")
        rows = self._read_rows(source)
        job.check_cancelled()

        job.report(0.38, "Normalizing table schema")
        normalized = []
        for index, row in enumerate(rows, start=1):
            record = dict(row)
            record.setdefault("event_id", f"event_{index:04d}")
            record["row_id"] = str(record["event_id"])
            record["abs_effect"] = abs(self._as_float(record.get("effect")))
            record["evidence_state"] = self._evidence_state(record)
            normalized.append(record)
        job.check_cancelled()

        directory = context.results.artifact_directory(
            self.workflow_id, run_fingerprint
        )
        table_path = directory / "event_profile.parquet"
        row_count = context.results.write_table(normalized, table_path)
        job.report(0.68, "Writing plot specification")
        plot_spec = self._plot_spec(normalized)
        plot_path = context.results.write_plot_spec(
            plot_spec, directory / "event_type_overview.vl.json"
        )

        metadata: dict[str, Any] = {
            "workflow_version": self.version,
            "project_snapshot_fingerprint": context.snapshot.manifest_fingerprint,
            "source": str(source),
            "support_prioritization_only": True,
        }
        table = AnalysisArtifact(
            project_id=context.manifest.project_id,
            analysis_id=self.workflow_id,
            run_fingerprint=run_fingerprint,
            parameters=parameters,
            metadata=metadata,
            name="Phase 1 event profile",
            kind=ArtifactKind.TABLE,
            entity_type="splice_event",
            fingerprint=fingerprint_value({"run": run_fingerprint, "kind": "table"}),
            path=table_path,
            row_count=row_count,
        )
        plot = AnalysisArtifact(
            project_id=context.manifest.project_id,
            analysis_id=self.workflow_id,
            run_fingerprint=run_fingerprint,
            parameters=parameters,
            metadata=metadata,
            name="Event type overview",
            kind=ArtifactKind.PLOT_SPEC,
            entity_type="splice_event_summary",
            fingerprint=fingerprint_value({"run": run_fingerprint, "kind": "plot"}),
            path=plot_path,
            parent_artifact_ids=(table.artifact_id,),
        )
        job.report(0.88, "Registering artifacts")
        context.catalog.register_artifact(table)
        context.catalog.register_artifact(plot)
        context.catalog.export_manifest(context.manifest.project_id)
        return [table.artifact_id, plot.artifact_id]

    @staticmethod
    def _read_rows(path: Path) -> list[dict[str, Any]]:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    @staticmethod
    def _as_float(value: object) -> float:
        if not isinstance(value, (str, int, float)):
            return 0.0
        try:
            return float(value)
        except ValueError:
            return 0.0

    def _evidence_state(self, row: dict[str, Any]) -> str:
        significant = str(row.get("significant", "")).strip().lower() in {
            "1",
            "true",
            "yes",
        }
        supported = self._as_float(row.get("read_support")) >= 20
        if significant and supported:
            return "significant + supported"
        if significant:
            return "significant, low support"
        if supported:
            return "supported, not significant"
        return "insufficient evidence"

    @staticmethod
    def _plot_spec(rows: list[dict[str, Any]]) -> dict[str, Any]:
        counts = Counter(str(row.get("event_type") or "Unknown") for row in rows)
        values = [
            {"event_type": event_type, "events": count}
            for event_type, count in sorted(counts.items())
        ]
        return {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "description": "Non-analytical Phase 1 fixture used to verify the platform contract.",
            "autosize": {"type": "fit", "contains": "padding"},
            "width": "container",
            "height": 245,
            "data": {"values": values},
            "mark": {"type": "bar", "cornerRadiusEnd": 5, "color": "#0f766e"},
            "encoding": {
                "x": {
                    "field": "event_type",
                    "type": "nominal",
                    "sort": "-y",
                    "axis": {"title": "Event type", "labelAngle": 0},
                },
                "y": {
                    "field": "events",
                    "type": "quantitative",
                    "axis": {"title": "Events in fixture"},
                },
                "tooltip": [
                    {"field": "event_type", "type": "nominal"},
                    {"field": "events", "type": "quantitative"},
                ],
            },
            "config": {
                "background": "transparent",
                "view": {"stroke": None},
                "axis": {
                    "domainColor": "#64748b",
                    "gridColor": "#dfe5e3",
                    "labelColor": "#24313a",
                    "titleColor": "#24313a",
                    "labelFont": "IBM Plex Sans",
                    "titleFont": "IBM Plex Sans",
                },
            },
        }
