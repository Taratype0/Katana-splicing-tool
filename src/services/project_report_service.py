from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from src.services.artifact_registry import AnalysisArtifactRegistry


class ProjectReportService:
    """Build a lightweight, reproducible report from cached module artifacts."""

    def __init__(self, output_root: str | Path, project_name: str) -> None:
        self.output_root = Path(output_root).resolve()
        self.project_name = str(project_name or "Katana project")
        self.report_dir = self.output_root / "00_project_report"

    def build(
        self,
        *,
        registry: AnalysisArtifactRegistry,
        project_summary: dict[str, object],
        module_states: dict[str, dict[str, object]],
    ) -> dict[str, Path]:
        self.report_dir.mkdir(parents=True, exist_ok=True)
        records = registry.find()
        inventory = pd.DataFrame([record.to_dict() for record in records])
        inventory_path = self.report_dir / "analysis_artifact_inventory.tsv"
        inventory.to_csv(inventory_path, sep="\t", index=False)

        state_rows = [
            {
                "module_id": module_id,
                **dict(details or {}),
            }
            for module_id, details in module_states.items()
        ]
        states = pd.DataFrame(state_rows)
        states_path = self.report_dir / "module_status.tsv"
        states.to_csv(states_path, sep="\t", index=False)

        generated_at = datetime.now(timezone.utc).isoformat()
        provenance = {
            "schema_version": 1,
            "generated_at": generated_at,
            "project": project_summary,
            "artifact_registry": str(registry.path),
            "artifact_count": len(records),
            "module_states": module_states,
        }
        provenance_path = self.report_dir / "analysis_provenance.json"
        provenance_path.write_text(
            json.dumps(provenance, indent=2, ensure_ascii=False, default=str),
            encoding="utf-8",
        )

        report_path = self.report_dir / "KATANA_PROJECT_REPORT.html"
        report_path.write_text(
            self._render_html(
                generated_at=generated_at,
                project_summary=project_summary,
                inventory=inventory,
                states=states,
            ),
            encoding="utf-8",
        )
        return {
            "report": report_path,
            "inventory": inventory_path,
            "module_status": states_path,
            "provenance": provenance_path,
        }

    def _render_html(
        self,
        *,
        generated_at: str,
        project_summary: dict[str, object],
        inventory: pd.DataFrame,
        states: pd.DataFrame,
    ) -> str:
        project_rows = "".join(
            f"<div><dt>{html.escape(str(key).replace('_', ' ').title())}</dt>"
            f"<dd>{html.escape(str(value))}</dd></div>"
            for key, value in project_summary.items()
        )
        status_counts = (
            states["status"].fillna("not_run").astype(str).value_counts().to_dict()
            if "status" in states.columns
            else {}
        )
        artifact_counts = (
            inventory["entity_type"].fillna("other").astype(str).value_counts().to_dict()
            if "entity_type" in inventory.columns
            else {}
        )
        status_cards = "".join(
            f"<article><strong>{html.escape(str(value))}</strong><span>{html.escape(str(key))}</span></article>"
            for key, value in sorted(status_counts.items())
        ) or "<article><strong>0</strong><span>cached modules</span></article>"
        artifact_cards = "".join(
            f"<article><strong>{html.escape(str(value))}</strong><span>{html.escape(str(key))} artifacts</span></article>"
            for key, value in sorted(artifact_counts.items())
        ) or "<article><strong>0</strong><span>registered artifacts</span></article>"
        state_table = self._frame_html(states, ["module_id", "status", "message", "output_folder"])
        artifact_table = self._frame_html(
            inventory,
            [
                "module_id",
                "artifact_key",
                "entity_type",
                "scope_id",
                "row_count",
                "status",
                "path",
                "source_artifact_ids",
            ],
        )
        return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(self.project_name)} | Katana report</title>
<style>
:root {{ --ink:#17211c; --muted:#647067; --paper:#f4f2ea; --card:#fffef9; --line:#d9d6ca; --accent:#0b6b53; --gold:#d99a2b; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; color:var(--ink); background:linear-gradient(135deg,#e9efe8 0,#f4f2ea 38%,#f7eddd 100%); font:15px/1.5 "Segoe UI",sans-serif; }}
main {{ width:min(1320px,94vw); margin:36px auto 72px; }}
header {{ padding:32px; border:1px solid var(--line); border-radius:22px; background:rgba(255,254,249,.94); box-shadow:0 20px 55px rgba(43,54,46,.09); }}
.eyebrow {{ color:var(--accent); font-size:12px; font-weight:800; letter-spacing:.14em; text-transform:uppercase; }}
h1 {{ margin:8px 0 4px; font:700 clamp(30px,4vw,54px)/1.05 Georgia,serif; }}
.muted {{ color:var(--muted); }}
dl {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:10px 26px; margin:24px 0 0; }}
dl div {{ padding-top:10px; border-top:1px solid var(--line); }}
dt {{ color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.08em; }}
dd {{ margin:3px 0 0; font-weight:650; overflow-wrap:anywhere; }}
section {{ margin-top:22px; padding:26px; border:1px solid var(--line); border-radius:18px; background:var(--card); }}
h2 {{ margin:0 0 16px; font:700 26px/1.2 Georgia,serif; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; }}
.cards article {{ padding:17px; border-radius:14px; color:white; background:var(--accent); }}
.cards article:nth-child(even) {{ background:#254d43; }}
.cards strong,.cards span {{ display:block; }}
.cards strong {{ font-size:27px; }}
.cards span {{ opacity:.8; }}
.table-wrap {{ overflow:auto; max-height:520px; border:1px solid var(--line); border-radius:12px; }}
table {{ width:100%; border-collapse:collapse; font-size:13px; background:white; }}
th {{ position:sticky; top:0; z-index:1; padding:10px; color:white; background:#26362f; text-align:left; }}
td {{ padding:9px 10px; border-top:1px solid #ece9df; vertical-align:top; max-width:430px; overflow-wrap:anywhere; }}
tr:nth-child(even) td {{ background:#faf8f1; }}
</style>
</head>
<body><main>
<header>
  <div class="eyebrow">Katana reproducible analysis report</div>
  <h1>{html.escape(self.project_name)}</h1>
  <div class="muted">Generated {html.escape(generated_at)} from cached outputs and AnalysisArtifactRegistry. Building this report does not run analyses.</div>
  <dl>{project_rows}</dl>
</header>
<section><h2>Module readiness</h2><div class="cards">{status_cards}</div></section>
<section><h2>Artifact handoff</h2><div class="cards">{artifact_cards}</div></section>
<section><h2>Module status</h2><div class="table-wrap">{state_table}</div></section>
<section><h2>Provenance and lineage</h2><p class="muted">Each row records the producing module, scope, file, fingerprint and upstream artifact IDs.</p><div class="table-wrap">{artifact_table}</div></section>
</main></body></html>"""

    @staticmethod
    def _frame_html(frame: pd.DataFrame, preferred_columns: list[str]) -> str:
        if frame.empty:
            return "<p>No cached records are available.</p>"
        columns = [column for column in preferred_columns if column in frame.columns]
        if not columns:
            columns = list(frame.columns)
        display = frame[columns].copy()
        for column in display.columns:
            display[column] = display[column].map(
                lambda value: ", ".join(str(item) for item in value)
                if isinstance(value, (list, tuple))
                else value
            )
        return display.to_html(index=False, border=0, escape=True)
