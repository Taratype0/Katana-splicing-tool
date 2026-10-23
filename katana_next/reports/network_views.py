"""Complete views of an immutable observed network; no network discovery here.

Optional matched-score tests and module ORA call the existing frozen kernels.
Neither operation changes event membership, PSI values or discovery parameters.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re

VERSION = "observed-network-report-v1"
BLUE, RED = "#2863A5", "#C63F4A"


def matched_module_associations(scores, samples, contrasts):
    """Run the frozen binary-trait adapter within each exact declared contrast.

    All three production rows per module and their original BH family survive.
    The display subset is numerator versus denominator, never an eight-group
    condition-versus-rest replacement. Unsupported designs remain unavailable.
    """
    import pandas as pd
    from katana_cosplice.trait_adapter import safe_trait_association

    if not scores.index.is_unique or not samples.sample_id.is_unique:
        raise ValueError("Unique frozen-score and metadata sample IDs are required")
    ids = [str(row["contrast_id"]) for row in contrasts]
    if len(ids) != len(set(ids)):
        raise ValueError("Declared contrast IDs must be unique")
    full, selected, audits = [], [], []
    for item in contrasts:
        cid = str(item["contrast_id"])
        arms = {key: list(item[key]["sample_ids"]) for key in ("numerator", "denominator")}
        audit = {"contrast_id": cid, "sample_ids": arms,
                 "effect_direction": "numerator minus denominator",
                 "score_orientation": "Frozen module-PC sign; not inclusion direction",
                 "discovery_repeated": False}
        try:
            a, b = arms["numerator"], arms["denominator"]
            if not a or not b or len(set(a)) != len(a) or len(set(b)) != len(b) or set(a) & set(b):
                raise ValueError("Contrast arms must be nonempty, disjoint and unique")
            wanted = a + b
            if not set(wanted) <= set(scores.index) or not set(wanted) <= set(samples.sample_id):
                raise ValueError("A declared comparison sample is absent from frozen scores or metadata")
            metadata = samples.set_index("sample_id").loc[wanted].reset_index().copy()
            metadata["declared_contrast_arm"] = ["numerator"] * len(a) + ["denominator"] * len(b)
            units = next((key for key in ("donor_id", "biological_sample_id") if key in metadata), None)
            # The frozen adapter resets its aligned metadata index. An unnamed
            # TSV index otherwise becomes an unrelated ``index`` column.
            # Name-only normalization preserves every sample ID and score.
            frozen_scores = scores.loc[wanted].rename_axis("sample_id")
            result = safe_trait_association(frozen_scores, metadata,
                factor="declared_contrast_arm", unit_column=units,
                batch_column=item.get("batch_column"))
            frame = result.condition_tests.copy()
            frame.insert(0, "contrast_id", cid)
            frame["correction_family_id"] = cid + "::frozen_binary_trait_family"
            full.append(frame)
            # A blocked/nuisance adapter returns omnibus terms, not a signed
            # Welch contrast. Keep its rows but never invent a signed effect.
            if {"contrast", "level"} <= set(frame):
                chosen = frame.loc[frame["contrast"].eq("numerator_vs_rest") & frame.level.eq("numerator")].copy()
                if len(chosen) != scores.shape[1]:
                    raise ValueError("Frozen adapter did not return one numerator Welch row per module")
                selected.append(chosen)
                audit["signed_effect_status"] = "available"
            else:
                audit["signed_effect_status"] = "not_computed: frozen blocked model emits omnibus terms only"
            for kind, extra in (("nuisance", result.nuisance_tests), ("coefficient", result.coefficients)):
                if len(extra):
                    extra = extra.copy()
                    extra.insert(0, "contrast_id", cid)
                    extra["row_scope"] = kind
                    full.append(extra)
            audit.update({"status": "completed", "adapter_audit": result.audit})
        except (ValueError, KeyError) as exc:
            audit.update({"status": "not_computed", "reason": str(exc)})
        audits.append(audit)
    return (pd.concat(full, ignore_index=True) if full else pd.DataFrame(),
            pd.concat(selected, ignore_index=True) if selected else pd.DataFrame(), audits)


def _fingerprint(path):
    path = Path(path)
    stat = path.stat()
    digest = sha256()
    with path.open("rb") as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(part)
    return {"path": str(path.resolve()), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": digest.hexdigest()}


def _slug(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)) + "-" + sha256(str(value).encode()).hexdigest()[:8]


def _completed_network_scope(completion_path, input_paths, contrasts):
    """Authenticate a native lifecycle without relabelling an imported failure.

    The recovery scope accepts only byte-verified saved numerical tables plus
    a separately completed postprocessing job. It never calls discovery or
    infers that the failed source job completed its missing plotting stages.
    """
    proof = json.loads(completion_path.read_text(encoding="utf-8"))
    if proof.get("schema") == "katana.flow.network_view_inputs/v1":
        from katana_next.flow.presentation import validate_network_view_proof
        return validate_network_view_proof(completion_path, input_paths, contrasts)
    def native_text(key):
        value = proof.get(key, "")
        if isinstance(value, dict):
            if value.get("returncode") != 0 or not isinstance(value.get("stdout"), str):
                raise ValueError("Native accounting command did not complete successfully")
            return value["stdout"]
        if not isinstance(value, str):
            raise ValueError("Native accounting must retain its actual text or successful command receipt")
        return value
    scope = proof.get("completion_scope")
    if scope == "observed_network_baseline_only":
        required = ("network_observed", "associate_observed", "annotate_modules")
        if proof.get("status") != "completed" or any(proof.get("stage_snapshot", {}).get(key, {}).get("status") not in {"finished", "cache_hit"} for key in required):
            raise ValueError("All exact baseline stages must have completed")
        job = str(proof.get("baseline_job", ""))
        rows = [row.split("|") for row in native_text("actual_sacct").splitlines() if row.split("|")[0] == job]
        if not job or len(rows) != 1 or rows[0][1:3] != ["COMPLETED", "0:0"]:
            raise ValueError("Native baseline COMPLETED / 0:0 evidence is required")
        return scope, job, None, [completion_path]
    if scope != "observed_numerical_import_and_postprocessing_only" or proof.get("schema") != "katana.network_numeric_recovery_native_acceptance/v1" or proof.get("passed") is not True:
        raise ValueError("A completed observed-network baseline or explicit verified numerical-import receipt is required")
    job = str(proof.get("native_job_id", ""))
    rows = [row.split("|") for row in native_text("actual_sacct").splitlines() if row.split("|")[0] == job]
    if not job or len(rows) != 1 or rows[0][1:3] != ["COMPLETED", "0:0"]:
        raise ValueError("Native postprocessing COMPLETED / 0:0 evidence is required")
    if proof.get("source_failed_status_preserved") is not True or proof.get("fit_network_rerun") is not False or proof.get("module_bootstrap_accepted") is not False:
        raise ValueError("Recovery must retain its failed source and separate unaccepted bootstrap scope")
    root = completion_path.parent
    verified = {}
    for item in proof.get("files", []):
        path = (root / item["path"]).resolve()
        if not path.is_relative_to(root) or path in verified:
            raise ValueError("Unbounded or duplicate recovered network member")
        current = _fingerprint(path)
        if (current["bytes"], current["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Recovered network member SHA identity mismatch")
        verified[path] = current
    phase_path = root / "postprocessing/completion.receipt.json"
    import_path = root / "IMPORT_PROVENANCE.json"
    plan_path = root / "PLAN.json"
    if not {phase_path, import_path, plan_path, *input_paths} <= set(verified):
        raise ValueError("Network inputs and full completion/import/plan proofs must be authenticated members")
    phase = json.loads(phase_path.read_text(encoding="utf-8"))
    provenance = json.loads(import_path.read_text(encoding="utf-8"))
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if (phase.get("schema") != "katana.network_verified_numeric_postprocessing/v1"
            or phase.get("status") != "completed" or phase.get("completion_scope") != scope
            or str(phase.get("native_job_id")) != job or phase.get("fit_network_rerun") is not False
            or phase.get("source_failed_status_preserved") is not True
            or phase.get("frozen_kernels_unchanged") is not True):
        raise ValueError("The authenticated network postprocessing phase is incomplete or from another job")
    if (phase.get("source_import_manifest_sha256") != verified[import_path]["sha256"]
            or phase.get("plan_sha256") != verified[plan_path]["sha256"]
            or proof.get("plan_sha256") != verified[plan_path]["sha256"]):
        raise ValueError("Network phase is not bound to its exact import and plan")
    source_job = str(plan.get("source_job_id", ""))
    source_rows = [row.split("|") for row in native_text("source_actual_sacct").splitlines() if row.split("|")[0] == source_job]
    if not source_job or len(source_rows) != 1 or source_rows[0][1:3] != ["FAILED", "1:0"]:
        raise ValueError("The original failed network lifecycle must remain explicit")
    if (provenance.get("schema") != "katana.verified_failed_plot_numerical_import/v1"
            or provenance.get("status") != "verified_import"
            or str(provenance.get("source_job_id")) != source_job
            or provenance.get("fit_network_rerun") is not False):
        raise ValueError("An explicit verified failed-plot numerical import is required")
    if plan.get("design", {}).get("contrasts") != contrasts or phase.get("declared_contrasts") != [row["contrast_id"] for row in contrasts]:
        raise ValueError("Precomputed network tests must match every exact declared contrast and sample arm")
    for item in phase.get("files", []):
        path = (root / item["path"]).resolve()
        if path not in verified or (verified[path]["bytes"], verified[path]["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Network postprocessing file is absent from the native collection")
    names = {"full_tests": "trait_full_family.tsv", "selected_tests": "trait_selected.tsv",
             "trait_audit": "trait_audit.json", "tested": "module_enrichment.tsv.gz",
             "background": "module_background.tsv", "input_genes": "module_input_genes.tsv"}
    precomputed = {key: root / "postprocessing" / name for key, name in names.items()}
    if not set(precomputed.values()) <= set(verified):
        raise ValueError("The complete precomputed association and enrichment family is required")
    precomputed["provenance"] = provenance
    return scope, job, precomputed, [completion_path, *verified]


def _write_browser(out, members, matrix, scores, samples, enrichment, edges_path):
    """A full-member, paged module browser with no event or module truncation."""
    import pandas as pd
    rows = members.astype(object).where(pd.notna(members), None).to_dict("records")
    payload = {"members": rows, "samples": list(matrix.columns),
               "psi": matrix.to_numpy(dtype=float).tolist(),
               "modules": list(dict.fromkeys(members.module.astype(str))),
               "score_samples": list(scores.index), "score_modules": list(scores.columns),
               "scores": scores.to_numpy(dtype=float).tolist(),
               "conditions": samples.set_index("sample_id").reindex(scores.index).condition.astype(str).tolist()
                   if "condition" in samples else ["unspecified"] * len(scores)}
    payload["edges"] = []
    if edges_path:
        native_edges = pd.read_csv(edges_path, sep="\t")
        index = {key: i for i, key in enumerate(members.event_key)}
        if not {"event_a", "event_b", "adjacency"} <= set(native_edges):
            raise ValueError("Native exported network edges require event_a, event_b and adjacency")
        if not (set(native_edges.event_a) | set(native_edges.event_b)) <= set(index):
            raise ValueError("An exported network edge refers to an unknown event")
        payload["edges"] = [[index[row.event_a], index[row.event_b], row.adjacency] for row in native_edges.itertuples()]
    data = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    page = r'''<!doctype html><html lang="en"><meta charset="utf-8"><title>All network modules</title>
<style>body{font:14px Arial;color:#243444;margin:28px}select,input,button{padding:7px;margin:4px}table{border-collapse:collapse;width:100%;font-size:12px}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:left;max-width:460px;overflow-wrap:anywhere}thead{position:sticky;top:0;background:#fff}canvas{max-width:100%;border:1px solid #ddd}small{color:#526170}a{color:#2863a5}.scroll{max-height:530px;overflow:auto}h1,h2{font-weight:500}</style>
<h1>Katana-Splice-Network: all observed modules</h1>
<p>Every assigned and unassigned event is retained. Select a module, search an event or gene, and page through all members. PSI colours show original values from 0 (blue) to 1 (red); no row scaling. Module-PC signs are arbitrary. Connectivity is the frozen <strong>global adjacency sum</strong>, not intramodule kWithin. Bootstrap confidence is a separate phase.</p>
<p><a href="../tables/all_module_members.tsv">Complete members and continuous metrics</a> · <a href="../tables/all_event_psi.tsv.gz">Complete PSI matrix</a> · <a href="../tables/module_activity.tsv">Frozen activity</a> · <a href="../tables/module_summary.tsv">All module summaries</a> · <a href="../tables/module_event_type_counts.tsv">Event types</a> · <a href="../tables/module_enrichment.tsv.gz">All tested GO / KEGG terms</a> · <a href="../tables/enrichment_background.tsv">Measured background</a> · <a href="../tables/module_input_genes.tsv">Module input genes</a> · <a href="../tables/matched_module_tests_full_family.tsv">All matched test families</a></p>
<label>Module <select id="module"></select></label><label>Gene / event <input id="query" type="search"></label>
<button id="prev">Previous</button><button id="next">Next</button><span id="count"></span>
<h2>Complete member PSI, 60 rows per page</h2><canvas id="psi"></canvas>
<h2>Complete members and metric ranks</h2><div class="scroll"><table><thead id="head"></thead><tbody id="body"></tbody></table></div>
<h2>Selected module activity in all biological samples</h2><canvas id="activity"></canvas>
<h2>Selected module network structure</h2><p>All selected-module nodes and all available native exported within-module edges are drawn in a deterministic circular layout. The native output-edge list can have a display cap recorded by discovery; this is not the full dense adjacency. Click a node for its identity. Node positions do not encode new statistical evidence.</p><canvas id="network"></canvas><p id="node"></p>
<p id="extra"></p><script id="data" type="application/json">PAYLOAD</script>
<script>
const d=JSON.parse(document.querySelector('#data').textContent),sel=document.querySelector('#module'),query=document.querySelector('#query');let page=0;const size=60;
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
for(const m of ['all',...d.modules]){const o=document.createElement('option');o.value=m;o.textContent=m==='grey'?'grey — unassigned':m;sel.append(o)}
function draw(){let ids=d.members.map((r,i)=>i).filter(i=>(sel.value==='all'||d.members[i].module===sel.value)&&JSON.stringify(d.members[i]).toLowerCase().includes(query.value.toLowerCase()));page=Math.max(0,Math.min(page,Math.ceil(ids.length/size)-1));const shown=ids.slice(page*size,(page+1)*size),keys=Object.keys(d.members[0]||{});document.querySelector('#count').textContent=`${ids.length} events; page ${page+1} / ${Math.max(1,Math.ceil(ids.length/size))}`;
document.querySelector('#head').innerHTML='<tr>'+keys.map(k=>'<th>'+esc(k)+'</th>').join('')+'</tr>';document.querySelector('#body').innerHTML=shown.map(i=>'<tr>'+keys.map(k=>'<td>'+esc(d.members[i][k])+'</td>').join('')+'</tr>').join('');
const c=document.querySelector('#psi'),ctx=c.getContext('2d');c.width=1200;c.height=120+shown.length*14;ctx.fillStyle='#fff';ctx.fillRect(0,0,c.width,c.height);const left=260,w=(c.width-left-20)/Math.max(1,d.samples.length);ctx.font='10px Arial';d.samples.forEach((s,j)=>{ctx.save();ctx.translate(left+j*w+6,108);ctx.rotate(-.65);ctx.fillStyle='#243444';ctx.fillText(s,0,0);ctx.restore()});shown.forEach((i,k)=>{ctx.fillStyle='#243444';ctx.fillText(String(d.members[i].gene_symbol||d.members[i].gene_id||d.members[i].event_key).slice(0,36),2,130+k*14);d.psi[i].forEach((v,j)=>{ctx.fillStyle=v<=.5?`rgb(${Math.round(40+(245-40)*2*v)},${Math.round(99+(245-99)*2*v)},${Math.round(165+(245-165)*2*v)})`:`rgb(${Math.round(245+(198-245)*2*(v-.5))},${Math.round(245+(63-245)*2*(v-.5))},${Math.round(245+(74-245)*2*(v-.5))})`;ctx.fillRect(left+j*w,119+k*14,w,13)})});
const ac=document.querySelector('#activity'),ax=ac.getContext('2d');ac.width=1200;ac.height=300;ax.fillStyle='#fff';ax.fillRect(0,0,1200,300);const mi=d.score_modules.indexOf(sel.value);if(mi<0){ax.fillStyle='#526170';ax.fillText('Select an assigned module to view all sample activities.',20,30)}else{let v=d.scores.map(r=>r[mi]),mx=Math.max(...v.map(Math.abs),.001);ax.beginPath();ax.moveTo(40,120);ax.lineTo(1180,120);ax.strokeStyle='#bbb';ax.stroke();v.forEach((z,i)=>{let x=45+i*1120/v.length;ax.fillStyle=z>=0?'#c63f4a':'#2863a5';ax.fillRect(x,Math.min(120,120-z/mx*100),900/v.length,Math.abs(z/mx*100));ax.save();ax.translate(x,240);ax.rotate(-.7);ax.fillStyle='#243444';ax.font='10px Arial';ax.fillText(d.score_samples[i]+' | '+d.conditions[i],0,0);ax.restore()})}
document.querySelector('#extra').innerHTML=sel.value==='all'?'':`<a href="../tables/module_input_genes.tsv">Full input genes for ${esc(sel.value)}</a>. All measured network genes, including grey, are the enrichment background.`}
sel.onchange=query.oninput=()=>{page=0;draw()};document.querySelector('#prev').onclick=()=>{page--;draw()};document.querySelector('#next').onclick=()=>{page++;draw()};draw();
</script></html>'''
    graph = r'''
function drawNetwork(){const c=document.querySelector('#network'),ax=c.getContext('2d');c.width=1000;c.height=760;ax.fillStyle='#fff';ax.fillRect(0,0,c.width,c.height);c.points=[];
if(sel.value==='all'){ax.fillStyle='#526170';ax.font='14px Arial';ax.fillText('Select one module to view every node and its exported edges.',20,35);return}
let ids=d.members.map((r,i)=>i).filter(i=>d.members[i].module===sel.value),points=new Map();ids.forEach((i,k)=>{let a=k/ids.length*2*Math.PI;points.set(i,[500+340*Math.cos(a),380+340*Math.sin(a)])});let n=0;ax.strokeStyle='#2863a522';ax.lineWidth=.5;for(const [i,j,w]of d.edges){if(points.has(i)&&points.has(j)){let a=points.get(i),b=points.get(j);ax.beginPath();ax.moveTo(...a);ax.lineTo(...b);ax.stroke();n++}}ax.fillStyle='#c63f4a';for(const i of ids){let [x,y]=points.get(i);ax.beginPath();ax.arc(x,y,Math.max(.6,Math.min(4,150/Math.sqrt(ids.length))),0,2*Math.PI);ax.fill();c.points.push([x,y,i])}ax.fillStyle='#243444';ax.font='13px Arial';ax.fillText(`${sel.value}: ${ids.length} nodes; ${n} available within-module edges`,20,25);if(!d.edges.length)ax.fillText('No native edge list supplied; no substitute edges are inferred.',20,48)}
const originalDraw=draw;draw=function(){originalDraw();drawNetwork()};document.querySelector('#network').onclick=e=>{const c=e.currentTarget,r=c.getBoundingClientRect(),x=(e.clientX-r.left)*c.width/r.width,y=(e.clientY-r.top)*c.height/r.height;let best=null,dist=100;for(const p of c.points){let q=(p[0]-x)**2+(p[1]-y)**2;if(q<dist){best=p;dist=q}}document.querySelector('#node').textContent=best?JSON.stringify(d.members[best[2]]):'Click near a node to inspect its complete identity and metrics.'};drawNetwork();
'''
    page = page.replace("</script></html>", graph + "</script></html>")
    (out / "browser" / "ALL_MODULES.html").write_text(page.replace("PAYLOAD", data), encoding="utf-8")


def build_network_assets(*, assignments_path, scores_path, psi_path, annotation_path,
                         samples_path, contrasts, completion_manifest, output_dir,
                         libraries=None, project_id=None, edges_path=None):
    """Build complete observed-network views after a native completion receipt.

    ``libraries`` is an explicit name-to-GMT-path mapping; existing species and
    database versions are recorded. No network or upstream computation runs.
    """
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    completion_path = Path(completion_manifest).resolve()
    input_paths = [Path(p).resolve() for p in (assignments_path, scores_path, psi_path, annotation_path, samples_path)]
    if edges_path:
        input_paths.append(Path(edges_path).resolve())
    completion_scope, job, precomputed, proof_paths = _completed_network_scope(completion_path, input_paths, contrasts)
    paths = list(dict.fromkeys([*input_paths, *proof_paths]))
    libraries = {str(k): Path(v).resolve() for k, v in (libraries or {}).items()}
    paths += list(libraries.values()) + ([Path(edges_path).resolve()] if edges_path else [])
    out = Path(output_dir).resolve()
    if out.exists() or any(out == p or out in p.parents or p.parent == out for p in paths):
        raise ValueError("Use a new independent output directory outside source files")
    audit = [_fingerprint(p) for p in paths]
    assignments = pd.read_csv(assignments_path, sep="\t", dtype={"event_key": str, "module": str})
    annotation = pd.read_csv(annotation_path, sep="\t", usecols=lambda c: c in {"event_key", "gene_id", "gene_symbol", "gene_name", "event_type"}, dtype=str).fillna("")
    scores = pd.read_csv(scores_path, sep="\t", index_col=0)
    psi = pd.read_csv(psi_path, sep="\t", index_col=0)
    samples = pd.read_csv(samples_path, sep="\t", dtype=str)
    if not assignments.event_key.is_unique or not annotation.event_key.is_unique or not psi.index.is_unique:
        raise ValueError("Unique exact event keys are required")
    if set(assignments.event_key) != set(psi.index):
        raise ValueError("Assignments and frozen PSI must contain the same complete event population")
    if not set(assignments.event_key) <= set(annotation.event_key):
        raise ValueError("An assigned event is absent from its original event annotation")
    if set(scores.index) != set(psi.columns) or set(scores.index) != set(samples.sample_id):
        raise ValueError("Frozen activity, PSI and sample metadata must retain the same complete samples")
    if not np.isfinite(psi.to_numpy()).all() or not np.isfinite(scores.to_numpy()).all():
        raise ValueError("Finite frozen PSI and activity are required; no implicit imputation")
    if not ((psi.to_numpy() >= 0) & (psi.to_numpy() <= 1)).all():
        raise ValueError("Expected PSI units from zero to one")
    members = assignments.merge(annotation, on="event_key", how="left", validate="one_to_one")
    if "gene_symbol" not in members:
        members["gene_symbol"] = members.get("gene_name", pd.Series("", index=members.index))
    members["gene_symbol"] = members.gene_symbol.fillna("").astype(str)
    members["event_type"] = members.get("event_type", pd.Series("unrecorded", index=members.index)).fillna("unrecorded")
    members = members.sort_values(["module", "event_key"], kind="stable").reset_index(drop=True)
    psi = psi.reindex(index=members.event_key, columns=scores.index)
    modules = list(members.module.drop_duplicates())
    if set(scores.columns) != set(modules) - {"grey"}:
        raise ValueError("Frozen scores must cover every assigned module and exclude grey")
    summary = members.groupby("module", sort=False).agg(n_events=("event_key", "size"), n_genes=("gene_symbol", lambda x: x[x.ne("")].nunique())).reset_index()
    summary["assignment_status"] = np.where(summary.module.eq("grey"), "unassigned", "assigned")
    summary["event_fraction"] = summary.n_events / len(members)
    types = members.groupby(["module", "event_type"], sort=False).size().rename("n_events").reset_index()
    genes = members.loc[members.gene_symbol.ne(""), ["module", "gene_symbol"]].drop_duplicates().sort_values(["module", "gene_symbol"])
    if precomputed:
        tested = pd.read_csv(precomputed["tested"], sep="\t")
        background = pd.read_csv(precomputed["background"], sep="\t")
        full_tests = pd.read_csv(precomputed["full_tests"], sep="\t")
        selected_tests = pd.read_csv(precomputed["selected_tests"], sep="\t")
        trait_audit = json.loads(precomputed["trait_audit"].read_text(encoding="utf-8"))
        declared = [row["contrast_id"] for row in contrasts]
        if [row.get("contrast_id") for row in trait_audit] != declared or not set(full_tests.contrast_id) <= set(declared) or not set(selected_tests.contrast_id) <= set(declared):
            raise ValueError("Precomputed network association populations differ from declared comparisons")
        if set(background.gene_symbol) != set(members.loc[members.gene_symbol.ne(""), "gene_symbol"]):
            raise ValueError("Precomputed network background differs from all measured event genes")
        supplied_genes = pd.read_csv(precomputed["input_genes"], sep="\t").fillna("")
        if set(map(tuple, supplied_genes[["module", "gene_symbol"]].to_numpy())) != set(map(tuple, genes.to_numpy())):
            raise ValueError("Precomputed network module input genes differ from frozen memberships")
        if len(tested) and not set(tested.module.astype(str)) <= set(modules) - {"grey"}:
            raise ValueError("Precomputed enrichment refers to an unobserved module")
    else:
        from katana_next.workflows.module_functions import enrich_modules
        tested, background = enrich_modules(assignments, annotation, libraries) if libraries else (pd.DataFrame(columns=["module", "library", "term", "pvalue", "FDR"]), pd.DataFrame({"gene_symbol": sorted(members.loc[members.gene_symbol.ne(""), "gene_symbol"].unique())}))
        full_tests, selected_tests, trait_audit = matched_module_associations(scores, samples, contrasts)
    for name in ("tables", "figures", "browser"):
        (out / name).mkdir(parents=True, exist_ok=True)
    artifacts = []
    common = {"workflow_version": VERSION, "completion_verified": True,
              "completion_scope": completion_scope, "completion_manifest": str(completion_path),
              "native_job_id": job, "native_status": "COMPLETED", "native_exit_code": "0:0",
              "source_discovery_scope": "observed modules; bootstrap confidence is a separate phase",
              "network_input_policy": {"events": len(members), "samples": len(scores), "max_events": None},
              "network_label_policy": {"hard_filter": False, "apply_threshold_labels": False},
              "portable_root": str(out), "new_computation": "report views only; native-completed associations and ORA imported unchanged" if precomputed else "report views; explicit frozen trait tests and module ORA only"}
    flow_complete = completion_scope == "flow_verified_network_complete"
    if precomputed and not flow_complete:
        common["numerical_import_source_job"] = precomputed["provenance"]["source_job_id"]
        common["source_lifecycle"] = "FAILED/1:0 preserved; numerical outputs independently verified"
    elif flow_complete:
        common.pop("native_job_id", None)
        common.pop("native_status", None)
        common.pop("native_exit_code", None)
        common.update(execution_backend="verified portable Nextflow task receipts",
                      source_discovery_scope="Completed observed network and independently completed module bootstrap",
                      flow_execution_plan_id=precomputed["provenance"]["flow_plan_id"],
                      new_computation="Read-only report views; completed contrast and ORA families imported unchanged")

    def asset(path, title, caption, *, module=None, display=True, kind="file", contrast=None, extra=None):
        relative = path.relative_to(out).as_posix()
        meta = dict(common, caption=caption, display_in_main=display)
        if module is not None:
            meta["module_id"] = str(module)
        meta.update(extra or {})
        row = {"artifact_id": "network-view-" + sha256(relative.encode()).hexdigest()[:18],
               "analysis_id": "programs.cosplice", "run_fingerprint": out.name, "name": title,
               "kind": kind, "path": str(path), "status": "complete", "parameters": {},
               "contrast_ids": [contrast] if contrast else [], "metadata": meta}
        if project_id:
            row["project_id"] = project_id
        artifacts.append(row)

    def table(frame, stem, title, caption, *, compress=False, source=None):
        path = out / "tables" / (stem + (".tsv.gz" if compress else ".tsv"))
        if source:
            import shutil
            shutil.copyfile(source, path)
        else:
            frame.to_csv(path, sep="\t", index=False, compression="gzip" if compress else None)
        asset(path, title, caption, display=False, kind="table")

    table(members, "all_module_members", "All module members and continuous metrics", "Every observed event, including unassigned grey. network_connectivity is the global adjacency sum, not kWithin. Original continuous metrics are preserved.")
    table(psi.reset_index(names="event_key"), "all_event_psi", "Complete module-member PSI", "All events and samples in original zero-to-one units; no row scaling.", compress=True)
    table(scores.reset_index(names="sample_id"), "module_activity", "Frozen observed module activity", "Unchanged source module-PC activity, with arbitrary PC orientation.")
    table(summary, "module_summary", "All module sizes and gene composition", "All modules and grey; unique genes counted once within a module.")
    table(types, "module_event_type_counts", "Complete module event-type composition", "Every event contributes once to its assigned module and event type.")
    table(genes, "module_input_genes", "Module enrichment input genes", "Unique measured event genes for every module; grey retained for audit and background.")
    table(tested, "module_enrichment", "All module GO / KEGG tests", "Frozen hypergeometric ORA; all measured network genes including grey form background. BH across all modules and terms within each explicit library. No gene-set or term truncation.", compress=True, source=precomputed["tested"] if precomputed else None)
    table(background, "enrichment_background", "Measured network-gene background", "Unique genes mapped from the complete measured network-event population, including unassigned grey; not all reference genes.", source=precomputed["background"] if precomputed else None)
    table(full_tests, "matched_module_tests_full_family", "Full declared-contrast module test families", "All original frozen binary-adapter rows and their untouched P/FDR; separate correction family for each declared comparison.", source=precomputed["full_tests"] if precomputed else None)
    table(selected_tests, "matched_module_numerator_minus_denominator", "Matched module effects", "Exact declared sample arms. Display numerator-minus-denominator Welch rows; their FDR is retained from the complete original three-row-per-module family.", source=precomputed["selected_tests"] if precomputed else None)
    if edges_path:
        import shutil
        dst = out / "tables" / "native_network_edges.tsv"
        shutil.copyfile(edges_path, dst)
        asset(dst, "Native exported network edges", "Unchanged frozen exported edge list; its recorded output-edge cap does not change the full network used for discovery and connectivity.", display=False, kind="table")
    (out / "tables" / "matched_trait_audit.json").write_text(json.dumps(trait_audit, indent=2) + "\n", encoding="utf-8")
    cmap = LinearSegmentedColormap.from_list("katana_blue_red", [BLUE, "#F4F5F7", RED])

    def save(fig, stem, title, caption, *, data, unit, full_table=None, **options):
        plot_table = out / "tables" / (stem + "-plot-data.tsv")
        data.to_csv(plot_table, sep="\t", index=False)
        contrast = options.get("contrast") or "all"
        links = [{"role": "plot_data", "label": "Plot data", "path": str(plot_table),
                  "contrast_id": contrast, "policy": None, "unit": unit}]
        if full_table:
            links.append({"role": "full_results", "label": "Complete corresponding results",
                          "path": str(out / "tables" / full_table),
                          "contrast_id": "all", "policy": None, "scope": "full_collection", "unit": unit})
        if stem == "module_enrichment_overview":
            links += [{"role": role, "label": label, "path": str(out / "tables" / name),
                       "contrast_id": "all", "policy": None, "unit": "gene"}
                      for role, label, name in (("foreground", "Module input genes", "module_input_genes.tsv"),
                                                ("background", "Measured network-gene background", "enrichment_background.tsv"))]
        path = out / "figures" / (stem + ".png")
        fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
        svg = path.with_suffix(".svg")
        fig.savefig(svg, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        asset(path, title, caption, extra={"exports": {"png": str(path), "svg": str(svg)},
                                          "data_links": links}, **options)

    fig, axes = plt.subplots(1, 2, figsize=(max(12, len(modules) * .38), 5), layout="constrained")
    for ax, column, label in ((axes[0], "n_events", "Events"), (axes[1], "n_genes", "Unique measured genes")):
        ax.bar(summary.module, summary[column], color=["#A5ADB7" if m == "grey" else BLUE for m in summary.module])
        ax.set(ylabel=label, xlabel="Observed module (grey = unassigned)")
        ax.tick_params(axis="x", rotation=70, labelsize=8)
    save(fig, "all_module_sizes", "All module sizes and gene composition", "Every observed module and unassigned set; no module is excluded from the overview.", data=summary, unit="module", full_table="module_summary.tsv")
    pivot = types.pivot(index="module", columns="event_type", values="n_events").fillna(0).reindex(modules)
    fig, ax = plt.subplots(figsize=(max(11, len(modules) * .35), 5), layout="constrained")
    pivot.plot.bar(stacked=True, ax=ax, color=[BLUE, RED, "#218D87", "#9A64B0", "#D49A25"], width=.8)
    ax.set(xlabel="All observed modules", ylabel="Events")
    ax.legend(title="Event type", frameon=False, bbox_to_anchor=(1, 1))
    save(fig, "all_module_event_types", "Event types in every module", "Counts from complete native event identities; grey is the unassigned set.", data=pivot.reset_index(), unit="module_event_type", full_table="module_event_type_counts.tsv")
    ordered_samples = samples.set_index("sample_id").reindex(scores.index)
    if "condition" in ordered_samples:
        order = np.argsort(ordered_samples.condition.to_numpy(), kind="stable")
    else:
        order = np.arange(len(scores))
    if len(scores.columns):
        plot = scores.iloc[order].T
        limit = max(float(np.abs(plot.to_numpy()).max()), .001)
        fig, ax = plt.subplots(figsize=(12, max(4, len(plot) * .3)), layout="constrained")
        image = ax.imshow(plot, cmap=cmap, vmin=-limit, vmax=limit, aspect="auto", interpolation="nearest")
        ax.set(xticks=range(len(plot.columns)), xticklabels=plot.columns, yticks=range(len(plot)), yticklabels=plot.index)
        ax.tick_params(axis="x", rotation=70, labelsize=7)
        fig.colorbar(image, ax=ax, label="Frozen module-PC score")
        save(fig, "all_module_activity", "Module activity in all biological samples", "All modules and samples, ordered by declared condition for display. PC signs are arbitrary; positive activity is not increased exon inclusion.", data=plot.rename_axis("module").reset_index(), unit="module_sample_score", full_table="module_activity.tsv")
    metrics = [key for key in ("network_connectivity", "module_membership", "module_membership_margin", "stability_score") if key in members]
    if metrics:
        fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), max(4, len(modules) * .22)), layout="constrained", squeeze=False)
        for ax, key in zip(axes[0], metrics):
            values = [pd.to_numeric(members.loc[members.module.eq(m), key], errors="coerce").dropna().to_numpy() for m in modules]
            ax.boxplot(values, vert=False, showfliers=False)
            ax.set_yticks(range(1, len(modules) + 1), modules)
            ax.set(xlabel=key)
        save(fig, "all_module_continuous_metrics", "Continuous membership and connectivity distributions", "All members retained. Boxplot whiskers summarize the complete distribution; all event values and ranks are in the member table. Global network_connectivity is not intramodule kWithin; optional threshold labels do not hard-filter members.", data=members[["event_key", "module", *metrics]], unit="event", full_table="all_module_members.tsv")
    if len(selected_tests):
        effects = selected_tests.pivot(index="module", columns="contrast_id", values="effect").reindex(scores.columns)
        limit = max(float(np.abs(effects.to_numpy()).max()), .001)
        fig, ax = plt.subplots(figsize=(max(11, len(contrasts) * .6), max(4, len(effects) * .3)), layout="constrained")
        image = ax.imshow(effects, cmap=cmap, vmin=-limit, vmax=limit, aspect="auto")
        ax.set(xticks=range(len(effects.columns)), xticklabels=effects.columns, yticks=range(len(effects)), yticklabels=effects.index)
        ax.tick_params(axis="x", rotation=65, labelsize=7)
        fig.colorbar(image, ax=ax, label="Numerator − denominator module-PC mean")
        save(fig, "matched_module_effects", "Module activity: declared experimental comparisons", "Only exact contrastsheet arms are compared. Colour shows a frozen-PC mean difference; full original Welch P/FDR and the complete three-row-per-module correction families are exported. No biological rescue label is inferred.", data=effects.rename_axis("module").reset_index(), unit="module_contrast", full_table="matched_module_tests_full_family.tsv")
    if len(tested):
        # Main report display only. Full tested terms are written above.
        shown = tested.sort_values(["FDR", "pvalue", "term"], kind="stable").groupby(["module", "library"], sort=False).head(2).copy()
        fig, ax = plt.subplots(figsize=(12, max(5, len(shown) * .2)), layout="constrained")
        shown = shown.sort_values(["module", "library", "FDR", "term"], kind="stable")
        x = -np.log10(shown.FDR.clip(lower=1e-300))
        ax.scatter(x, range(len(shown)), s=20 + 12 * shown.overlap_count, c=[BLUE if q < .05 else "#A5ADB7" for q in shown.FDR], alpha=.8)
        ax.set(yticks=range(len(shown)), yticklabels=[f"{r.module} | {r.library} | {r.term}"[:125] for r in shown.itertuples()], xlabel="−log10(native BH FDR)")
        ax.tick_params(axis="y", labelsize=6)
        plotted_terms = shown.assign(plot_neglog10_fdr=x.to_numpy(), plot_point_area=20 + 12 * shown.overlap_count.to_numpy())
        save(fig, "module_enrichment_overview", "GO / KEGG coverage across all modules", "Display: two lowest-FDR terms per module and library, including nonsignificant terms in grey. Point area follows overlap count. Complete tests, all input and contributing genes, actual measured background and unchanged P/FDR are exported; this figure does not filter modules.", data=plotted_terms, unit="module_term_test", full_table="module_enrichment.tsv.gz")
    for module in modules:
        subset = members.loc[members.module.eq(module)]
        if module in scores.columns:
            fig, ax = plt.subplots(figsize=(12, 4), layout="constrained")
            values = scores.loc[:, module].iloc[order]
            ax.bar(range(len(values)), values, color=[BLUE if v < 0 else RED for v in values])
            ax.set(xticks=range(len(values)), xticklabels=values.index, ylabel="Frozen module-PC activity")
            ax.tick_params(axis="x", rotation=65, labelsize=7)
            save(fig, "module_activity_" + _slug(module), f"{module}: complete sample activity", "All biological samples; module-PC orientation remains unchanged from discovery.", data=values.rename_axis("sample_id").rename("activity").reset_index(), unit="module_sample_score", full_table="module_activity.tsv", module=module, display=False)
        if "network_connectivity" in subset:
            rank = subset.sort_values(["network_connectivity", "event_key"], ascending=[False, True], kind="stable").head(20)
            fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
            labels = [f"{r.gene_symbol or r.event_key} | {r.event_type}"[:65] for r in rank.itertuples()]
            ax.barh(range(len(rank)), rank.network_connectivity, color=BLUE)
            ax.set(yticks=range(len(rank)), yticklabels=labels, xlabel="Global adjacency connectivity (not intramodule kWithin)")
            ax.invert_yaxis()
            save(fig, "module_connectivity_" + _slug(module), f"{module}: connectivity-ranked members", "Display: first 20 by unchanged global connectivity, ties resolved by event key for display. Every member, original metric and rank remains in the complete table and paged browser; this is not a threshold-defined Hub claim.", data=rank, unit="event", full_table="all_module_members.tsv", module=module, display=False)
    _write_browser(out, members, psi, scores, samples, tested, edges_path)
    asset(out / "browser" / "ALL_MODULES.html", "All network modules: members, PSI and activity", "All modules and grey, all members, searchable and paginated full PSI, continuous metrics and complete table links. No event is removed for display.", extra={"report_asset_kind": "html"})
    for row in audit:
        current = _fingerprint(row["path"])
        if current != row:
            raise RuntimeError("Original completed network source changed during rendering")
    receipt = {"schema": "katana.observed_network_report/v1", "version": VERSION, "passed": True,
               "observed_utc": datetime.now(timezone.utc).isoformat(), "source_files": audit,
               "events": len(members), "samples": len(scores), "assigned_modules": len(scores.columns),
               "unassigned_events": int(members.module.eq("grey").sum()), "module_ora_rows": len(tested),
               "enrichment_background_genes": len(background), "declared_contrasts": len(contrasts),
               "matched_test_audit": trait_audit,
               "native_baseline_job": job if not precomputed else None,
               "native_postprocessing_job": job if precomputed else None,
               "original_source_native_job": precomputed["provenance"].get("source_job_id") if precomputed else job,
               "bootstrap_accepted": flow_complete, "network_discovery_runs": 0, "upstream_jobs": 0,
               "sources_preserved": True, "figures": sum(a["path"].endswith(".png") for a in artifacts),
               "completion_scope": completion_scope,
               "new_trait_fits": 0 if precomputed else len(contrasts), "new_module_ora": not bool(precomputed),
               "enrichment_method": "existing frozen enrich_modules: hypergeometric; BH all modules × terms within library",
               "trait_method": "existing frozen safe_trait_association per exact declared sample arms"}
    (out / "network-report.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    (out / "artifact-records.json").write_text(json.dumps(artifacts, indent=2) + "\n", encoding="utf-8")
    unavailable = [] if flow_complete else [{"chapter": 6, "title": "Bootstrap module confidence", "status": "not_computed",
                    "reason": "The observed baseline receipt does not establish completion of the separate bootstrap phase."}]
    if precomputed and not flow_complete:
        unavailable.append({"chapter": 6, "title": "Original network dendrogram", "status": "not_available",
                            "reason": "The source plot failed after complete numerical writes; original linkage was not saved. No replacement fit or tree is inferred."})
    unavailable += [{"chapter": 6, "title": "Matched module association: " + row["contrast_id"],
                     "status": "not_computed", "reason": row["reason"]} for row in trait_audit if row["status"] != "completed"]
    return {"artifacts": artifacts, "unavailable": unavailable, "receipt": receipt}
