from __future__ import annotations

import json
from html import escape
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from katana_cosplice.branding import MODULE_NAME


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def read_cached_table(path: Path) -> pd.DataFrame:
    return _read_table(path)


def _records(frame: pd.DataFrame, limit: int = 50_000) -> list[dict[str, object]]:
    if frame.empty:
        return []
    cleaned = frame.head(limit).where(pd.notna(frame), None)
    return cleaned.to_dict(orient="records")


def build_module_browser(output_dir: Path, project_name: str, target: Path) -> Path:
    """Build a self-contained, read-only browser from cached CoSplice artifacts."""
    assignments = _read_table(output_dir / "04c_consensus_modules" / "consensus_event_module_assignments.tsv")
    if assignments.empty:
        assignments = _read_table(output_dir / "04_network" / "observed" / "event_module_assignments.tsv")
    hubs = _read_table(output_dir / "06e_module_interpretation" / "module_hub_events.tsv")
    regulators = _read_table(output_dir / "06e_module_interpretation" / "module_regulator_prioritization.tsv")
    enrichment = _read_table(output_dir / "06e_module_interpretation" / "module_go_kegg_enrichment.tsv")
    isoforms = _read_table(output_dir / "06e_module_interpretation" / "module_isoform_switch_enrichment.tsv")
    preservation = _read_table(output_dir / "06d_external_preservation" / "external_module_preservation.tsv")
    module_column = next(
        (name for name in ("consensus_module", "module") if name in assignments),
        "module",
    )
    payload = {
        "assignments": _records(assignments),
        "hubs": _records(hubs),
        "regulators": _records(regulators),
        "enrichment": _records(enrichment),
        "isoforms": _records(isoforms),
        "preservation": _records(preservation),
        "moduleColumn": module_column,
    }
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{escape(project_name)} module browser</title>
<style>
body{{font-family:Arial,sans-serif;margin:24px;color:#15242b;background:#f3f7f6}}
h1{{color:#087e8b}}.controls{{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}}
select,input{{font-size:15px;padding:8px 10px;min-width:220px;border:1px solid #9fb4b7;border-radius:6px}}
.tabs button{{padding:9px 14px;margin-right:6px;border:0;border-radius:6px;background:#dce9e8;cursor:pointer}}
.tabs button.active{{background:#087e8b;color:white}}.note{{color:#53656d}}
.table-wrap{{overflow:auto;max-height:68vh;background:white;border:1px solid #c8d7d7}}
table{{border-collapse:collapse;width:100%;font-size:13px}}th{{position:sticky;top:0;background:#dce9e8}}
th,td{{padding:7px 9px;border-bottom:1px solid #e1e9e8;text-align:left;white-space:nowrap}}
</style></head><body>
<h1>{MODULE_NAME} module browser: {escape(project_name)}</h1>
<p class="note">Read-only browser generated from cached outputs. Module, hub, RBP, enrichment, isoform-switch and preservation evidence are prioritization layers, not causal claims.</p>
<div class="controls"><label>Module <select id="module"></select></label><label>Search <input id="search" placeholder="event, gene, RBP, term..."></label></div>
<div class="tabs" id="tabs"></div><p id="status" class="note"></p><div class="table-wrap"><table id="table"></table></div>
<script>
const DATA={json.dumps(payload, ensure_ascii=True)};
const views=['assignments','hubs','regulators','enrichment','isoforms','preservation'];
let active='assignments';
const moduleKeys=['consensus_module','module','reference_module'];
function valueModule(row){{for(const key of moduleKeys) if(row[key]!==undefined&&row[key]!==null) return String(row[key]); return '';}}
function modules(){{const out=new Set(); for(const view of views) for(const row of DATA[view]){{const m=valueModule(row);if(m)out.add(m)}} return [...out].sort();}}
function init(){{const s=document.getElementById('module');s.innerHTML='<option value="">All modules</option>'+modules().map(x=>`<option>${{x}}</option>`).join('');
 const tabs=document.getElementById('tabs');tabs.innerHTML=views.map(v=>`<button data-view="${{v}}">${{v.replaceAll('_',' ')}}</button>`).join('');
 tabs.onclick=e=>{{if(e.target.dataset.view){{active=e.target.dataset.view;render()}}}};s.onchange=render;document.getElementById('search').oninput=render;render();}}
function render(){{const module=document.getElementById('module').value;const q=document.getElementById('search').value.toLowerCase();
 const rows=DATA[active].filter(r=>(!module||valueModule(r)===module)&&(!q||JSON.stringify(r).toLowerCase().includes(q)));
 document.querySelectorAll('#tabs button').forEach(b=>b.classList.toggle('active',b.dataset.view===active));
 const table=document.getElementById('table');if(!rows.length){{table.innerHTML='';document.getElementById('status').textContent='No cached rows found for this selection.';return;}}
 const cols=[...new Set(rows.flatMap(r=>Object.keys(r)))];table.innerHTML='<thead><tr>'+cols.map(c=>`<th>${{c}}</th>`).join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+cols.map(c=>`<td>${{r[c]??''}}</td>`).join('')+'</tr>').join('')+'</tbody>';
 document.getElementById('status').textContent=`${{rows.length}} cached rows shown.`;}}
init();
</script></body></html>""",
        encoding="utf-8",
    )
    return target


def build_validation_pdf(output_dir: Path, project_name: str, state: dict[str, object], target: Path) -> Path:
    """Render a compact PDF validation summary without recomputing any analysis."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(target) as pdf:
        fig = plt.figure(figsize=(11.7, 8.3))
        fig.text(0.07, 0.90, f"{MODULE_NAME} validation report", fontsize=24, weight="bold", color="#087e8b")
        fig.text(0.07, 0.84, project_name, fontsize=17, weight="bold")
        fig.text(0.07, 0.77, "Cached-output scientific audit", fontsize=13)
        fig.text(0.07, 0.70, "Modules are descriptive replicate-supported splicing programs.\nFunctional layers prioritize mechanisms; they do not establish causal regulation.", fontsize=11)
        stages = state.get("stages", {})
        y = 0.58
        for name, details in list(stages.items())[:18]:
            fig.text(0.08, y, str(name), fontsize=9, weight="bold")
            fig.text(0.30, y, str(details.get("status", "")), fontsize=9)
            fig.text(0.43, y, str(details.get("message", ""))[:95], fontsize=8)
            y -= 0.027
        plt.axis("off")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        summaries = [
            ("Consensus module sizes", output_dir / "04c_consensus_modules" / "consensus_module_summary.tsv", "consensus_module", "n_events"),
            ("Internal module stability", output_dir / "06_module_validation" / "module_bootstrap_stability.tsv", "module", "median_event_membership"),
            ("External preservation", output_dir / "06d_external_preservation" / "external_module_preservation.tsv", "consensus_module", "Zsummary_wgcna_like"),
            ("RBP regulator prioritization", output_dir / "06e_module_interpretation" / "module_regulator_prioritization.tsv", "rbp", "regulator_priority_score"),
            ("Isoform-switch support", output_dir / "06e_module_interpretation" / "module_isoform_switch_enrichment.tsv", "module", "n_switch_genes_in_module"),
        ]
        for title, path, label_column, value_column in summaries:
            frame = _read_table(path)
            if frame.empty or label_column not in frame or value_column not in frame:
                continue
            view = frame[[label_column, value_column]].copy()
            view[value_column] = pd.to_numeric(view[value_column], errors="coerce")
            view = view.dropna().sort_values(value_column, ascending=False).head(20)
            if view.empty:
                continue
            fig, ax = plt.subplots(figsize=(11.7, 8.3))
            ax.barh(view[label_column].astype(str)[::-1], view[value_column][::-1], color="#159a9c")
            ax.set_title(title, fontsize=18, weight="bold", pad=18)
            ax.set_xlabel(value_column)
            ax.spines[["top", "right"]].set_visible(False)
            fig.text(0.5, 0.02, f"Source: {path.name} (cached output)", ha="center", fontsize=8, color="#53656d")
            fig.tight_layout(rect=(0.06, 0.06, 0.98, 0.95))
            pdf.savefig(fig)
            plt.close(fig)
    return target
