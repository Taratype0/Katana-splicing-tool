"""Bounded auxiliary HTML for completed large MotifTrace statistical families.

This presentation fallback leaves every native data file and statistic intact.
The seven-chapter report uses its own approved renderer and is unaffected.
"""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path

SCIENTIFIC_FILES = {
    "ANALYSIS_PLAN_BEFORE_SCORING.json", "EXECUTION_BUDGET_BEFORE_SCORING.json", "FROZEN_SELECTED_EVENTS.tsv",
    "MATCHING_ALL_SELECTED_WITH_LOSSES.tsv", "ACTUAL_RNA_FORMS.tsv.gz", "ACTUAL_LOCAL_RNA.fa",
    "ACTUAL_LOCUS_HITS.tsv.gz", "PAIRED_FORM_INSTANCE_CHANGES.tsv.gz", "PER_REPLICATE_SEQUENCE_OPPORTUNITY.tsv.gz",
    "PAIRED_FORM_NORMALIZED_SCORES.tsv.gz", "NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz", "NEW_CLUSTER_FWER.tsv",
    "MAX_CLUSTER_NULL.npy", "RESULTS.json",
}
OVERSIZE_ERROR = "Portable event report exceeds 64MiB; complete tables retained, no truncation"


def render_large_family_index(summary, destination):
    destination = Path(destination)
    root = destination.parent
    required = {"ANALYSIS_PLAN_BEFORE_SCORING.json", "FROZEN_SELECTED_EVENTS.tsv", "MATCHING_ALL_SELECTED_WITH_LOSSES.tsv",
                "ACTUAL_RNA_FORMS.tsv.gz", "ACTUAL_LOCAL_RNA.fa", "ACTUAL_LOCUS_HITS.tsv.gz",
                "PAIRED_FORM_INSTANCE_CHANGES.tsv.gz", "PER_REPLICATE_SEQUENCE_OPPORTUNITY.tsv.gz",
                "PAIRED_FORM_NORMALIZED_SCORES.tsv.gz", "NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz",
                "NEW_CLUSTER_FWER.tsv", "MAX_CLUSTER_NULL.npy", "RESULTS.json"}
    if any(not (root / name).is_file() for name in required):
        raise ValueError("Large-family presentation requires all completed native scientific files")
    names = sorted(p.name for p in root.iterdir() if p.is_file() and p != destination and p.suffix != ".html")
    records = {}
    for name in names:
        path = root / name
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        records[name] = {"bytes": path.stat().st_size, "sha256": digest}
    payload = {"summary": {key: summary.get(key) for key in ("comparison_id", "selected_events", "selected_event_types",
                    "matched_independent_gene_pairs", "eligible_new_window_tests", "new_window_q_lt_0_05",
                    "cluster_P_FWER_lt_0_05", "inference_status")},
               "available_files": names + (["COMPLETE.json"] if "COMPLETE.json" not in names else []),
               "file_identities": records, "presentation": "complete_native_files_index",
               "native_statistics_modified": False, "statistical_family_split": False,
               "reason": "Complete test family exceeds the frozen inline HTML display budget; all native rows remain in files."}
    serialized = json.dumps(payload, ensure_ascii=True, allow_nan=False).replace("<", "\\u003c")
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Katana MotifTrace | Complete results</title><style>
body{font:16px/1.6 Arial,sans-serif;color:#1e2630;background:white;max-width:1120px;margin:40px auto;padding:0 24px}
h1,h2{color:#2863A5}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:10px;border-bottom:1px solid #d9dfe7;text-align:left}
pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f7fa;padding:20px}a{color:#2863A5}.note{border-left:3px solid #C63F4A;padding-left:16px}</style>
<h1>Katana MotifTrace: complete results</h1><h2>''' + html.escape(str(summary.get("comparison_id", ""))) + '''</h2>
<p class="note">This complete statistical family exceeds the inline HTML display budget. Every model, sequence, window test, original P/FDR, cluster FWER result and matching record is retained in the files below. No result was removed or recomputed. The one-click report provides the ranked motif views.</p>
<pre id="summary"></pre><h2>Complete native data</h2><table><thead><tr><th>File</th><th>Bytes</th><th>SHA-256</th></tr></thead><tbody id="files"></tbody></table>
<p>Window FDR and cluster FWER describe different statistical families. Motif sequence opportunity is not measured binding or a causal claim.</p>
<script>const D=''' + serialized + ''';
document.getElementById('summary').textContent=JSON.stringify(D.summary,null,2);
for(const n of D.available_files){const tr=document.createElement('tr'),a=document.createElement('a');a.href=encodeURIComponent(n);a.download=n;a.textContent=n;const td=document.createElement('td');td.append(a);tr.append(td);for(const value of [D.file_identities[n]?.bytes??'Written after report completion',D.file_identities[n]?.sha256??'See task completion receipt']){const cell=document.createElement('td');cell.textContent=value;tr.append(cell)}document.getElementById('files').append(tr)}
</script></html>'''
    destination.write_text(page, encoding="utf-8")


def render_budgeted_consequences(rna, changes, tests, clusters, summary, destination, *, paired_scores=None):
    """Keep the frozen native display when it fits; use full-file index on its exact size error."""
    from katana_next.science.motiftrace_consequences import render_consequences
    try:
        return render_consequences(rna, changes, tests, clusters, summary, destination, paired_scores=paired_scores)
    except ValueError as error:
        if str(error) != OVERSIZE_ERROR:
            raise
        return render_large_family_index(summary, destination)


def validate_recorded_consequences_inputs(source, complete, selected_ids, models, fasta, *,
        failure_receipt, permutations=999, seed=20261015, provenance=None,
        max_tests=None, max_matrix_cells=None):
    """Bind a presentation-only recovery to the current scientific request.

    Execution adapter source hashes may change to repair presentation. Every
    scientific kernel hash, input identity and statistical parameter must stay
    fixed. FASTA paths may relocate only when their bytes and index are equal.
    """
    from katana_next.flow.model import canonical_hash, file_sha256, read_json
    from katana_next.science.motiftrace_budget import execution_limits
    source, fasta = Path(source), Path(fasta)
    plan = read_json(source / "ANALYSIS_PLAN_BEFORE_SCORING.json")
    failure = read_json(failure_receipt)
    actual_table = hashlib.sha256(complete.to_json(orient="split", double_precision=15).encode()).hexdigest()
    if (actual_table != plan.get("input_event_table_sha256") or set(map(str, selected_ids)) != set(plan["selected_event_ids"]) or
        complete.comparison_id.nunique() != 1 or str(complete.comparison_id.iloc[0]) != plan["comparison_id"] or
        canonical_hash(models) != canonical_hash(plan["models"]) or permutations != plan["permutations"] or seed != plan["seed"] or
        plan.get("window_nt") != [10, 20, 50] or plan.get("step") != 1 or plan.get("minimum_independent_gene_pairs") != 20):
        raise ValueError("Presentation recovery request changes the completed scientific inputs or parameters")
    requested_budget = (plan.get("execution_budget") if max_tests is None and max_matrix_cells is None
                        else execution_limits(max_tests, max_matrix_cells))
    if requested_budget != plan.get("execution_budget"):
        raise ValueError("Presentation recovery must retain the recorded execution budgets")
    old, current = plan.get("provenance", {}), provenance or {}
    special = {"implementation_SHA256", "genome_FASTA", "flow_execution_adapter_sha256"}
    if {k: v for k, v in old.items() if k not in special} != {k: v for k, v in current.items() if k not in special}:
        raise ValueError("Presentation recovery changes scientific provenance")
    execution_files = {"motiftrace_execution.py", "motiftrace_budget.py", "motif_native_index.py"}
    old_sources, new_sources = old.get("implementation_SHA256", {}), current.get("implementation_SHA256", {})
    if not old_sources or ({k: v for k, v in old_sources.items() if k not in execution_files} !=
                           {k: v for k, v in new_sources.items() if k not in execution_files}):
        raise ValueError("Presentation recovery changes frozen scientific implementation identity")
    for original, actual in ((old.get("genome_FASTA"), fasta), (str(old.get("genome_FASTA")) + ".fai", Path(str(fasta) + ".fai"))):
        records = [r for r in failure.get("inputs", []) if r.get("path") == original]
        if len(records) != 1 or file_sha256(actual) != records[0]["sha256"]:
            raise ValueError("Presentation recovery FASTA/index differ from the actual failed-source inputs")
    interface_changes = {name: {"before": old_sources.get(name), "after": new_sources.get(name)}
                         for name in execution_files if old_sources.get(name) != new_sources.get(name)}
    if old.get("flow_execution_adapter_sha256") != current.get("flow_execution_adapter_sha256"):
        interface_changes["flow_execution_adapter_sha256"] = {
            "before": old.get("flow_execution_adapter_sha256"), "after": current.get("flow_execution_adapter_sha256")}
    return {"scientific_inputs_identical": True, "selected_events": len(selected_ids), "models": len(models),
            "permutations": permutations, "seed": seed, "scientific_calculations_rerun": 0,
            "execution_interface_changes": interface_changes}


def finalize_recorded_consequences(source, destination, source_evidence):
    """New presentation output from an exact, authenticated failed-render case.

    Caller must bind this operation to a new task lifecycle. The old FAILED
    receipt and original directory remain unchanged; no scientific function is
    called. This is not a general mechanism for accepting arbitrary partials.
    """
    import shutil
    import numpy as np
    import pandas as pd
    from katana_next.flow.model import canonical_hash, file_sha256, read_json, write_json
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if destination == source or destination.is_relative_to(source) or source.is_relative_to(destination):
        raise ValueError("Presentation recovery needs a new independent output directory")
    failure_path = Path(source_evidence["failure_receipt"]).resolve()
    if file_sha256(failure_path) != source_evidence["failure_receipt_sha256"]:
        raise ValueError("Failed-source receipt changed")
    failure = read_json(failure_path)
    if (failure.get("status") != "failed" or failure.get("error") != OVERSIZE_ERROR or
        failure.get("workflow_id") != "prioritize.motiftrace_consequences" or
        failure.get("task_id") != source_evidence["task_id"] or
        failure.get("task_contract_sha256") != source_evidence["task_contract_sha256"] or
        failure.get("receipt_sha256") != canonical_hash({k: v for k, v in failure.items() if k != "receipt_sha256"})):
        raise ValueError("Recovery requires the authentic exact native display-size failure")
    members = source_evidence["members"]
    if set(members) != SCIENTIFIC_FILES or (source / "COMPLETE.json").exists():
        raise ValueError("Recovery requires the complete expected scientific files of a failed presentation")
    before = {}
    for name, record in members.items():
        path = source / name
        if not path.is_file() or path.stat().st_size != record["bytes"] or file_sha256(path) != record["sha256"]:
            raise ValueError("Completed scientific file missing or changed: " + name)
        before[name] = (path.stat().st_mtime_ns, record["bytes"], record["sha256"])
    plan = read_json(source / "ANALYSIS_PLAN_BEFORE_SCORING.json")
    if file_sha256(source / "ANALYSIS_PLAN_BEFORE_SCORING.json") != source_evidence["scientific_plan_sha256"]:
        raise ValueError("The failed source's pre-scoring scientific plan changed")
    summary = read_json(source / "RESULTS.json")
    budget = read_json(source / "EXECUTION_BUDGET_BEFORE_SCORING.json")
    cid = plan["comparison_id"]
    if (failure.get("selection", {}).get("contrast_ids") != [cid] or summary.get("comparison_id") != cid or
        summary.get("policy") != plan or budget.get("limits") != plan.get("execution_budget") or
        plan.get("statistical_family_split") is not False):
        raise ValueError("Failed source does not bind one complete unchanged statistical family")
    read = lambda name: pd.read_csv(source / name, sep="\t", float_precision="round_trip", dtype={"event_id": str, "model_id": str})
    selected, matching = read("FROZEN_SELECTED_EVENTS.tsv"), read("MATCHING_ALL_SELECTED_WITH_LOSSES.tsv")
    if (selected.event_id.duplicated().any() or set(selected.event_id) != set(plan["selected_event_ids"]) or
        len(selected) != summary["selected_events"] or len(matching) != len(selected) or
        int(matching.status.eq("matched").sum()) != summary["matched_independent_gene_pairs"]):
        raise ValueError("Failed source selected-event or matching family is incomplete")
    tests = read("NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz")
    if (len(tests) != budget["estimate"]["scalar_and_window_tests"] or
        set(tests.model_id) != set(plan["models"]) or
        int(tests.raw_P_new_paired_swaps.notna().sum()) != summary["eligible_new_window_tests"] or
        int(tests.q_new_window_family.lt(.05).sum()) != summary["new_window_q_lt_0_05"]):
        raise ValueError("Failed source statistical table is not its complete declared family")
    for key in ("raw_P_new_paired_swaps", "q_new_window_family"):
        finite = tests[key].dropna()
        if not np.isfinite(finite).all() or not finite.between(0, 1).all():
            raise ValueError("Invalid native P/FDR in failed source")
    clusters = read("NEW_CLUSTER_FWER.tsv")
    if (len(clusters) != summary["cluster_rows"] or
        int(clusters.cluster_P_FWER_new_family.lt(.05).sum()) != summary["cluster_P_FWER_lt_0_05"]):
        raise ValueError("Failed source cluster table is incomplete")
    null = np.load(source / "MAX_CLUSTER_NULL.npy", allow_pickle=False)
    if null.ndim != 1 or len(null) != plan["permutations"] or not np.isfinite(null).all():
        raise ValueError("Failed source maximum-cluster null is incomplete")
    destination.mkdir(parents=True, exist_ok=False)
    for name in members:
        shutil.copyfile(source / name, destination / name)
    render_large_family_index(summary, destination / "MOTIFTRACE_EVENT_CONSEQUENCES.html")
    receipt = {"schema": "katana.motiftrace_presentation_recovery/v1", "source_task_status": "failed",
               "source_evidence": source_evidence, "comparison_id": cid, "scientific_calculations_rerun": 0,
               "native_statistics_modified": False, "scientific_family_split": False,
               "presentation_only": True}
    write_json(destination / "PRESENTATION_RECOVERY.json", receipt)
    for name, identity in before.items():
        path = source / name
        if (path.stat().st_mtime_ns, path.stat().st_size, file_sha256(path)) != identity:
            raise ValueError("Original failed source changed during presentation recovery")
    if file_sha256(failure_path) != source_evidence["failure_receipt_sha256"]:
        raise ValueError("Original failure receipt changed during presentation recovery")
    complete = {"complete": True, "statistical_family_split": False, "presentation_recovery": True,
                "quantification_rerun": False, "members": {}}
    for path in sorted(destination.iterdir()):
        if path.is_file():
            complete["members"][path.name] = {"bytes": path.stat().st_size, "sha256": file_sha256(path)}
    write_json(destination / "COMPLETE.json", complete)
    return summary
