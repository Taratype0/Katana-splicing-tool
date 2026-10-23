"""Generated-data execution-budget guards; frozen statistical kernels stay byte-identical."""
import ast
import copy
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_next.science import motiftrace_consequences as frozen
from katana_next.science import motiftrace_execution as adapter
from katana_next.science.motiftrace_budget import execution_limits, estimate_paired_family
from katana_next.science.motiftrace_clusters import build_paired_datasets
from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION


def _generated(tmp_path):
    sequence = "ACGT" * 100
    fasta = tmp_path / "generated.fa"
    fasta.write_text(">chr1\n" + sequence + "\n", newline="\n")
    Path(str(fasta) + ".fai").write_text("chr1\t400\t6\t400\t401\n")
    rows = []
    for i in range(24):
        for foreground in (True, False):
            rows.append(dict(event_id=f"{'f' if foreground else 'b'}{i:02}", gene_id=f"gene_{foreground}_{i}",
                comparison_id="B_vs_A", event_type="SE", strand="+", chr="chr1",
                dpsi=.3 if foreground else 0., fdr=.001 if foreground else .8,
                stable_significant=foreground, significant=foreground, measurement_eligible=True,
                read_support_policy=READ_SUPPORT_POLICY_VERSION, criterion_min_reads=10,
                numerator_mean_reads=100, denominator_mean_reads=100, IncLevel1=".6,.7,.65", IncLevel2=".3,.4,.35",
                upstreamES=0, upstreamEE=40, exonStart_0base=100, exonEnd=136, downstreamES=300, downstreamEE=340))
    models = {"consensus": {"kind": "consensus", "pattern": "ACG", "width": 3},
              "pwm": {"kind": "PWM", "matrix": [[.7,.1,.1,.1]] * 3, "width": 3}}
    frame = pd.DataFrame(rows)
    return frame, set(frame.loc[frame.stable_significant, "event_id"]), models, fasta


def test_adapter_ast_diff_is_execution_budget_only():
    baseline = ast.parse(inspect.getsource(frozen.analyze_consequences)).body[0]
    actual = ast.parse(inspect.getsource(adapter.analyze_consequences)).body[0]
    assert [a.arg for a in actual.args.kwonlyargs][-2:] == ["max_tests", "max_matrix_cells"]
    actual.args.kwonlyargs = actual.args.kwonlyargs[:-2]
    actual.args.kw_defaults = actual.args.kw_defaults[:-2]
    removable = {"limits", "estimate", "budget"}
    nodes = []
    for node in actual.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in removable for t in node.targets):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            if isinstance(node.value.func, ast.Name) and node.value.func.id == "enforce_paired_budget":
                continue
            if "EXECUTION_BUDGET_BEFORE_SCORING.json" in ast.unparse(node):
                continue
        nodes.append(node)
    actual.body = nodes
    for node in ast.walk(actual):
        if isinstance(node, ast.Dict):
            remove = {"execution_budget", "execution_adapter_version", "statistical_family_split"}
            pairs = [(k,v) for k,v in zip(node.keys, node.values) if not isinstance(k, ast.Constant) or k.value not in remove]
            node.keys, node.values = [p[0] for p in pairs], [p[1] for p in pairs]
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "build_paired_datasets":
            assert sum(k.arg is None and isinstance(k.value, ast.Name) and k.value.id == "limits" for k in node.keywords) == 1
            node.keywords = [k for k in node.keywords if k.arg is not None]
    assert ast.dump(actual, include_attributes=False) == ast.dump(baseline, include_attributes=False)


@pytest.mark.parametrize("value", [0, -1, 1.1, True, float("nan"), float("inf"), "100"])
def test_invalid_execution_budgets_rejected(value):
    with pytest.raises(ValueError):
        execution_limits(value, 8_000_000)


def test_preflight_counts_exact_full_family_and_fails_before_scoring(tmp_path, monkeypatch):
    frame, ids, models, fasta = _generated(tmp_path)
    def forbidden(*args, **kwargs):
        pytest.fail("Motif scoring began before insufficient-budget rejection")
    monkeypatch.setattr(adapter, "trace_instances", forbidden)
    with pytest.raises(ValueError, match="pre-scoring execution budget exceeded"):
        adapter.analyze_consequences(frame, ids, models, fasta, tmp_path / "rejected", max_tests=1, permutations=99)
    root = tmp_path / "rejected"
    budget = json.loads((root / "EXECUTION_BUDGET_BEFORE_SCORING.json").read_text())
    assert budget["estimate"]["models"] == 2
    assert budget["estimate"]["matched_gene_pairs"] == 24
    assert budget["estimate"]["scalar_and_window_tests"] > 1
    assert not (root / "ACTUAL_LOCUS_HITS.tsv.gz").exists()
    assert not (root / "COMPLETE.json").exists()
    matching = pd.read_csv(root / "MATCHING_ALL_SELECTED_WITH_LOSSES.tsv", sep="\t")
    rna = frozen.extract_event_rna(frame, fasta)
    datasets, _ = build_paired_datasets(rna, matching, models)
    assert budget["estimate"]["scalar_and_window_tests"] == sum(len(d["positions"]) for d in datasets)
    assert budget["estimate"]["paired_matrix_cells"] == sum(d["differences"].size for d in datasets)
    estimate = estimate_paired_family(rna, matching, models, widths=(10, 10, 500))
    datasets, _ = build_paired_datasets(rna, matching, models, widths=(10, 10, 500))
    assert estimate["scalar_and_window_tests"] == sum(len(d["positions"]) for d in datasets)


def test_two_sufficient_budgets_preserve_frozen_numeric_results(tmp_path):
    frame, ids, models, fasta = _generated(tmp_path)
    roots = [tmp_path / name for name in ("frozen", "bounded", "larger")]
    frozen.analyze_consequences(frame, ids, models, fasta, roots[0], permutations=99)
    adapter.analyze_consequences(frame, ids, models, fasta, roots[1], permutations=99, max_tests=100_000)
    adapter.analyze_consequences(frame, ids, models, fasta, roots[2], permutations=99, max_tests=1_000_000, max_matrix_cells=9_000_000)
    for source in roots[0].iterdir():
        if source.name.endswith((".tsv", ".tsv.gz")):
            expected = pd.read_csv(source, sep="\t", float_precision="round_trip")
            for destination in roots[1:]:
                pd.testing.assert_frame_equal(expected, pd.read_csv(destination/source.name, sep="\t", float_precision="round_trip"))
    for destination in roots[1:]:
        np.testing.assert_array_equal(np.load(roots[0]/"MAX_CLUSTER_NULL.npy"), np.load(destination/"MAX_CLUSTER_NULL.npy"))
    tests = pd.read_csv(roots[1]/"NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz", sep="\t")
    assert tests.raw_P_new_paired_swaps.notna().any()
    assert set(tests.model_id) == set(models)
    assert json.loads((roots[1]/"ANALYSIS_PLAN_BEFORE_SCORING.json").read_text())["execution_budget"]["max_tests"] == 100_000


def test_workflow_and_rbp_entry_expose_budget_api_parameters():
    from katana_next.workflows.registry import WorkflowRegistry
    for identity in ("prioritize.motiftrace_consequences", "prioritize.rbp_mechanism"):
        workflow = WorkflowRegistry().get(identity)
        assert workflow.normalize_parameters({})["max_tests"] == 1_000_000
        custom = workflow.normalize_parameters({"max_tests": 200_000, "max_matrix_cells": 2_000_000})
        assert custom["max_tests"] == 200_000 and custom["max_matrix_cells"] == 2_000_000
