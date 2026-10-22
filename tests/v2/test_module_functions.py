import json
import pandas as pd
import pytest
from katana_next.workflows.module_functions import enrich_modules
from katana_next.science.module_dashboard_report import dashboard_section


def test_module_gene_input_is_unique_and_background_includes_grey(tmp_path):
    assignments = pd.DataFrame({"event_key": ["e1", "e2", "e3", "e4"], "consensus_module": ["M1", "M1", "M1", "grey"]})
    events = pd.DataFrame({"event_key": assignments.event_key, "gene_symbol": ["A", "A", "B", "C"]})
    gmt = tmp_path/"BP.gmt"
    gmt.write_text("hit\tna\tA\tC\nzero_hit\tna\tC\n", encoding="utf-8")
    result, background = enrich_modules(assignments, events, {"BP": gmt})
    assert set(background.gene_symbol) == {"A", "B", "C"}
    assert result.input_genes.eq(2).all() and result.background_genes.eq(3).all()
    assert result.input_gene_list.eq("A;B").all()
    assert result.loc[result.term.eq("zero_hit"), "pvalue"].iloc[0] == 1
    assert result.FDR.notna().all()


def test_module_function_mapping_never_merges_event_identifiers(tmp_path):
    labels = pd.DataFrame({"event_key": ["duplicate", "duplicate"], "consensus_module": ["M1", "M2"]})
    with pytest.raises(ValueError, match="Unique"):
        enrich_modules(labels, pd.DataFrame({"event_key": ["duplicate"], "gene_symbol": ["A"]}), {})


def test_report_json_cannot_inject_script_and_contains_cached_measurements():
    payload = {"module": "</script><script>bad</script>", "event_measurements": {"events": ["e1"], "samples": [], "psi": []}}
    output = dashboard_section([payload])
    assert "</script><script>bad" not in output
    assert "event_measurements" in output and "\\u003c/script" in output
