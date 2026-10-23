"""Independent Gene / Event reports; never alter the seven-section renderer."""
from __future__ import annotations

from html import escape
import json
from pathlib import Path

import pandas as pd

from ._execution import PLOT_LOCK, portable, write_json


STYLE = """body{font:15px system-ui,sans-serif;color:#152d45;background:#f4f6f8;margin:0}
header{background:#102c46;color:white;padding:30px max(5%,24px)}h1{margin:0 0 10px;font-size:30px}
main{max-width:1280px;margin:25px auto;padding:0 24px}section,article{background:white;padding:24px;
margin:20px 0;border:1px solid #dde3e9;border-radius:3px}h2{font-size:21px;color:#102c46}
a{color:#2863A5}p,li{line-height:1.6}table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:9px;border-bottom:1px solid #e2e6eb;vertical-align:top;max-width:420px;
overflow-wrap:anywhere}th{background:#edf2f6}img{max-width:100%;height:auto}.scroll{overflow-x:auto}
.muted{color:#586b7d}.pill{display:inline-block;padding:4px 9px;background:#edf2f6;border-radius:3px}
.complete{color:#2863A5}.failed{color:#C63F4A}input,button,select{font:inherit;padding:7px;border:1px solid #c6cfd8;
background:white;border-radius:2px}button{cursor:pointer}nav{display:flex;gap:15px;flex-wrap:wrap;margin:14px 0}
pre{white-space:pre-wrap;overflow-wrap:anywhere;font:12px ui-monospace,monospace}.grid{display:grid;
grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px}.toolbar{display:flex;gap:12px;align-items:center;
margin:12px 0;flex-wrap:wrap}.toolbar input{min-width:260px}.summary td:first-child{width:190px;font-weight:600}
"""


def html_page(title, body, subtitle=""):
    return ("<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{escape(title)}</title><style>{STYLE}</style></head><body>"
            f"<header><h1>{escape(title)}</h1><p>{escape(subtitle)}</p></header>"
            f"<main>{body}</main></body></html>")


def _json_script(value):
    return json.dumps(portable(value), ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")


def render_table(path, title, frame, note=""):
    """All rows searchable and paginated offline, no remote JS/fetch required."""
    columns = [str(c) for c in frame.columns]
    rows = portable(frame.to_dict("records"))
    body = f"<p>{escape(note)}</p><p>{len(frame):,} complete rows. Search examines every row.</p>"
    body += """<div class='toolbar'><input id='search' type='search' placeholder='Search all rows' aria-label='Search all rows'>
    <select id='size' aria-label='Rows per page'><option>25</option><option>50</option><option>100</option><option>250</option></select>
    <button id='previous'>Previous</button><span id='page'></span><button id='next'>Next</button></div>
    <div class='scroll'><table><thead id='head'></thead><tbody id='rows'></tbody></table></div>"""
    body += f"<script id='table-data' type='application/json'>{_json_script({'columns': columns, 'rows': rows})}</script>"
    body += """<script>(()=>{'use strict';const data=JSON.parse(document.getElementById('table-data').textContent);
    let page=0;const search=document.getElementById('search'),size=document.getElementById('size');
    const cell=v=>v===null||v===undefined?'':typeof v==='object'?JSON.stringify(v):String(v);
    const tr=document.createElement('tr');for(const c of data.columns){const th=document.createElement('th');th.textContent=c;tr.append(th)}
    document.getElementById('head').append(tr);function draw(){const q=search.value.toLowerCase();
    const selected=data.rows.filter(r=>data.columns.some(c=>cell(r[c]).toLowerCase().includes(q)));
    const n=Number(size.value),pages=Math.max(1,Math.ceil(selected.length/n));page=Math.min(page,pages-1);
    const tbody=document.getElementById('rows');tbody.replaceChildren();for(const r of selected.slice(page*n,(page+1)*n)){
    const row=document.createElement('tr');for(const c of data.columns){const td=document.createElement('td');td.textContent=cell(r[c]);row.append(td)}tbody.append(row)}
    document.getElementById('page').textContent=`${selected.length} rows · Page ${page+1} / ${pages}`;
    document.getElementById('previous').disabled=page===0;document.getElementById('next').disabled=page===pages-1}
    search.addEventListener('input',()=>{page=0;draw()});size.addEventListener('change',()=>{page=0;draw()});
    document.getElementById('previous').addEventListener('click',()=>{page=Math.max(0,page-1);draw()});
    document.getElementById('next').addEventListener('click',()=>{page++;draw()});draw()})();</script>"""
    Path(path).write_text(html_page(title, body), encoding="utf-8")


def sample_measurements(request, target, tables):
    """Display native saved per-sample values; no model/statistical refit."""
    samples = {str(s["sample_id"]): s for s in request["project"].get("samples", [])}
    contrasts = {str(c["contrast_id"]): c for c in request["project"].get("contrasts", [])}
    rows = []
    seen = set()
    for item, frame in tables:
        for _, row in frame.iterrows():
            cid = str(row.get("comparison_id", row.get("contrast_id", item.get("metadata", {}).get("comparison_id", ""))))
            eid = str(row.get("event_id", row.get("event_uid", row.get("event_key", row.get("ID", "")))))
            # Wide sample matrices retain exact sample IDs and original units.
            for sid, sample in samples.items():
                if sid in frame.columns:
                    value = pd.to_numeric(pd.Series([row[sid]]), errors="coerce").iloc[0]
                    if pd.notna(value):
                        key = (item["evidence_id"], eid, cid, sid, "saved_matrix")
                        if key not in seen:
                            rows.append(dict(evidence_id=item["evidence_id"], event_id=eid, comparison_id=cid,
                                             sample_id=sid, condition=sample.get("condition"), measurement="saved_matrix",
                                             unit=item.get("metadata", {}).get("unit", "source units"), value=float(value)))
                            seen.add(key)
            # Original rMATS ordering is the declared contrast arm ordering.
            if cid in contrasts:
                contrast = contrasts[cid]
                for arm_number, arm_name in ((1, "numerator"), (2, "denominator")):
                    column = f"IncLevel{arm_number}"
                    if column not in frame.columns:
                        continue
                    raw = row[column]
                    if isinstance(raw, str):
                        values = raw.split(",")
                    elif isinstance(raw, (tuple, list)) or hasattr(raw, "tolist"):
                        values = list(raw)
                    else:
                        continue
                    ids = contrast.get(arm_name, {}).get("sample_ids", [])
                    flipped = any(str(row.get(flag, "")).lower() in {"1", "true", "yes"}
                                  for flag in ("splicing_direction_flipped", "reverse_direction"))
                    if flipped:
                        other_column = f"IncLevel{3 - arm_number}"
                        other = row.get(other_column)
                        if isinstance(other, str):
                            values = other.split(",")
                        elif isinstance(other, (tuple, list)) or hasattr(other, "tolist"):
                            values = list(other)
                    if len(ids) != len(values):
                        continue  # never guess which biological replicate is missing
                    for sid, value in zip(ids, values):
                        value = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
                        key = (item["evidence_id"], eid, cid, str(sid), "PSI")
                        if pd.notna(value) and key not in seen:
                            rows.append(dict(evidence_id=item["evidence_id"], event_id=eid, comparison_id=cid,
                                             sample_id=str(sid), condition=samples.get(str(sid), {}).get("condition"),
                                             measurement="PSI", unit="fraction", value=float(value)))
                            seen.add(key)
            # Already long-form sample data.
            if "sample_id" in frame.columns:
                sid = str(row["sample_id"])
                for column in ("psi", "PSI", "expression", "normalized_count", "TPM", "count", "value"):
                    if column not in frame.columns:
                        continue
                    value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
                    key = (item["evidence_id"], eid, cid, sid, column)
                    if sid in samples and pd.notna(value) and key not in seen:
                        rows.append(dict(evidence_id=item["evidence_id"], event_id=eid, comparison_id=cid,
                                         sample_id=sid, condition=samples[sid].get("condition"),
                                         measurement=column, unit=item.get("metadata", {}).get("unit", "source units"),
                                         value=float(value)))
                        seen.add(key)
    items = {item["evidence_id"]: item for item, _ in tables}
    for row in rows:
        sample = samples.get(row["sample_id"], {})
        item = items[row["evidence_id"]]
        attributes = sample.get("attributes", {})
        row.update(biological_sample_id=sample.get("biological_sample_id"),
                   donor_id=attributes.get("donor_id"), cell_type=attributes.get("cell_type"),
                   count_mode=item.get("metadata", {}).get("count_mode") or item.get("parameters", {}).get("count_mode") or "unknown",
                   source_statistical_family=item.get("source_fingerprint"))
    return pd.DataFrame(rows, columns=["evidence_id", "event_id", "comparison_id", "sample_id",
                                       "biological_sample_id", "donor_id", "cell_type", "condition",
                                       "measurement", "unit", "value", "count_mode", "source_statistical_family"])


def render_sample_plots(root, frame):
    if frame.empty:
        return []
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    root = Path(root)
    plots = []
    # Every event/measurement group is plotted; no arbitrary Top N selection.
    for i, ((evidence_id, mode, eid, cid, measurement, unit), group) in enumerate(
            frame.groupby(["evidence_id", "count_mode", "event_id", "comparison_id", "measurement", "unit"], sort=False, dropna=False)):
        relative = f"plots/sample-{i:05d}"
        (root / "plots").mkdir(exist_ok=True)
        with PLOT_LOCK:
            fig, ax = plt.subplots(figsize=(max(7, len(group) * .55), 4.2))
            conditions = list(dict.fromkeys(group["condition"].fillna("unknown").astype(str)))
            colors = {c: ("#2863A5", "#C63F4A", "#6F8498", "#9A5870")[j % 4] for j, c in enumerate(conditions)}
            ax.scatter(range(len(group)), group["value"], c=[colors[str(c)] for c in group["condition"].fillna("unknown")], s=42)
            ax.set_xticks(range(len(group)), group["sample_id"], rotation=45, ha="right")
            ax.set_ylabel(f"{measurement} ({unit})")
            ax.set_title(f"{eid or 'Selected gene'} | {cid or 'Saved sample values'} | {mode}")
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            for suffix in ("png", "svg"):
                fig.savefig(root / f"{relative}.{suffix}", dpi=170, facecolor="white")
            plt.close(fig)
        plots.append({"path": relative + ".png", "svg": relative + ".svg", "event_id": eid,
                      "comparison_id": cid, "measurement": measurement, "unit": unit,
                      "evidence_id": evidence_id, "count_mode": mode})
    return plots


def render_target(root, request, target, evidence, cards, samples, plots, sashimi, status):
    root = Path(root)
    label = target.get("gene_symbol") or target.get("gene_id") or target.get("event_id") or target["target_id"]
    title = f"Katana {target['target_kind'].title()} report | {label}"
    body = "<nav><a href='../../../INDEX.html'>Batch index</a><a href='gene-card.json'>Gene card JSON</a>"
    body += "<a href='event-cards.json'>All event cards JSON</a><a href='interactive-data.json'>Interactive data</a></nav>"
    body += f"<section><h2>Selected target</h2><p class='pill {escape(status)}'>{escape(status)}</p>"
    identity = [("Study", request["project"].get("name", request["project"]["project_id"])),
                ("Gene", f"{target.get('gene_symbol', '')} / {target.get('gene_id', '')}"),
                ("Target kind", target["target_kind"]), ("Exact event", target.get("event_id", "All exported events of the selected gene")),
                ("Comparisons", ", ".join(target.get("comparison_ids", []))),
                ("Reference", request["project"].get("references", {}).get("genome_build", "Not declared"))]
    body += "<table class='summary'>" + "".join(f"<tr><td>{escape(k)}</td><td>{escape(str(v))}</td></tr>" for k, v in identity) + "</table>"
    body += "<p>Effects, original P/FDR and test-family identities are inherited from completed source results. "
    body += "Target selection does not refit differential models, network discovery, enrichment or MotifTrace statistics.</p></section>"
    body += f"<section><h2>Gene / Event cards</h2><p>{len(cards):,} exact event/comparison cards. "
    body += "The gene card retains the full exported event set; no implicit three-event or Top N limit.</p>"
    body += "<a href='event-cards.html'>Browse every event card</a> · <a href='gene-card.json'>Complete gene card</a></section>"
    body += "<section><h2>Saved statistical, network and MotifTrace evidence</h2><div class='scroll'><table>"
    body += "<thead><tr><th>Analysis / evidence</th><th>Rows</th><th>Complete data</th><th>Native source</th></tr></thead><tbody>"
    for table in evidence:
        links = " · ".join(f"<a href='{escape(table[key], quote=True)}'>{text}</a>" for key, text in
                           (("html", "Browse all"), ("parquet", "Parquet"), ("tsv", "TSV"), ("json", "JSON")))
        body += (f"<tr><td>{escape(table['name'])}<br><span class='muted'>{escape(table['analysis_id'])}</span></td>"
                 f"<td>{table['rows']:,}</td><td>{links}</td><td><code>{escape(table['source_sha256'])}</code><br>"
                 f"{escape(str(table.get('metadata', {}).get('count_mode', 'unknown')))} · "
                 f"{escape(str(table.get('metadata', {}).get('completion_verification', 'exported completed source')))}</td></tr>")
    body += "</tbody></table></div><p>Network membership/activity and MotifTrace rows, when available, remain tied to their original populations, "
    body += "motif models and full statistical families. An absent analysis is unavailable, not a zero result.</p></section>"
    body += "<section><h2>Biological sample measurements</h2>"
    if plots:
        body += "<a href='sample-measurements.html'>Browse all sample values and biological identities</a>"
        for plot in plots:
            body += f"<figure><img src='{escape(plot['path'], quote=True)}' loading='lazy' alt='Saved per-sample {escape(plot['measurement'])}'><figcaption>"
            body += f"{escape(str(plot['event_id']))} | {escape(str(plot['comparison_id']))} · <a href='{escape(plot['svg'], quote=True)}'>SVG</a></figcaption></figure>"
    else:
        body += "<p>Saved per-sample values are unavailable. No sample values were inferred from summary statistics.</p>"
    body += "</section><section><h2>Exact-event native Sashimi</h2>"
    for record in sashimi:
        body += f"<article><h3>{escape(str(record.get('event_id', 'Requested event')))} | {escape(str(record.get('comparison_id', '')))}</h3>"
        body += f"<p class='{escape(record['status'])}'>{escape(record['status'])} · {escape(str(record.get('reader') or 'No BAM reader'))}</p>"
        if record.get("reason"):
            body += f"<p>{escape(record['reason'])}</p>"
        for path in record.get("files", []):
            if path.endswith(".png"):
                body += f"<img src='sashimi/{escape(path, quote=True)}' loading='lazy' alt='Native Sashimi event and biological replicates'>"
        body += "<p>" + " · ".join(f"<a href='sashimi/{escape(p, quote=True)}'>{escape(Path(p).name)}</a>" for p in record.get("files", [])) + "</p></article>"
    body += "</section><section><h2>Methods, source parameters and browsing settings</h2>"
    body += "<p>The portable request and unchanged selected evidence retain original source parameters separately from browsing thresholds. "
    body += "BAMs and large references are not included in the report download.</p><a href='target-provenance.json'>Complete provenance</a>"
    body += f"<pre>{escape(json.dumps(portable(target.get('source_selections', [])), indent=2, ensure_ascii=False))}</pre></section>"
    (root / "REPORT.html").write_text(html_page(title, body, "Independent selected-target dossier · separate from the seven-section study report"), encoding="utf-8")


def render_index(root, manifest):
    rows = []
    for target in manifest["targets"]:
        title = target.get("gene_symbol") or target.get("gene_id") or target["target_id"]
        link = f"<a href='{escape(target['report'], quote=True)}'>{escape(title)}</a>" if target.get("report") else escape(title)
        rows.append(f"<tr><td>{link}</td><td>{escape(target.get('target_kind', ''))}</td><td>{escape(str(target.get('event_id', '')))}</td>"
                    f"<td class='{escape(target['status'])}'>{escape(target['status'])}</td><td>{escape(target.get('reason', ''))}</td></tr>")
    body = (f"<section><h2>Study: {escape(manifest.get('project_name', manifest['project_id']))}</h2>"
            f"<p>{len(rows)} requested targets · {escape(manifest['status'])}</p>"
            "<nav><a href='result-manifest.json'>Verified result manifest</a><a href='request.json'>Portable input request</a>"
            "<a href='REPORT.zip'>Download complete target reports</a></nav>"
            "<p>Each target is independent. Completed source models are reused; missing read evidence is explicit. "
            "The seven-section study report and original analyses are unchanged.</p></section>"
            "<section><table><thead><tr><th>Gene / target</th><th>Kind</th><th>Exact event</th><th>Status</th><th>Diagnostic</th></tr></thead><tbody>"
            + "".join(rows) + "</tbody></table></section>")
    (Path(root) / "INDEX.html").write_text(html_page("Katana Gene / Event reports", body), encoding="utf-8")
