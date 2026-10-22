"""Report labels follow the declared study and preserve completed result tables."""
import json
from pathlib import Path

import pandas as pd
import pytest

from katana_next.science.narrative_report import build_report


def test_report_uses_selected_contrast_count_and_explicit_design_metadata(tmp_path, monkeypatch):
    from katana_next.science import report_annotation, report_complete

    monkeypatch.setattr(report_complete, "polished_heatmaps", lambda base, output_dir=None: None)
    monkeypatch.setattr(report_complete, "complete_sections", lambda *args: "")
    monkeypatch.setattr(report_annotation, "annotate_exons", lambda frame, annotation: frame.copy())
    base = tmp_path / "completed-results"
    presentation = base / "20_Presentation_report"
    expression = base / "22_Expression_interaction_pathways" / "expression-design" / "results"
    exons = base / "16_Interaction" / "exon-design"
    for directory in (presentation, expression, exons):
        directory.mkdir(parents=True)
    pd.DataFrame({"contrast": ["selected-a", "selected-b"], "significant_events": [3, 7], "DEG": [2, 4]}).to_csv(presentation / "comparison_totals.tsv", sep="\t", index=False)
    pd.DataFrame({"contrast": ["selected-a"], "overlap": [1]}).to_csv(presentation / "DEG_DAS_gene_intersection.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": ["GENE_A", "GENE_B"], "gene_symbol": ["A", "B"], "padj": [0.01, 0.5]}).to_csv(expression / "interaction.tsv", sep="\t", index=False)
    pd.DataFrame({"contrast": ["interaction"], "study_family_FDR": [0.2], "pathway": ["KEGG_EXAMPLE"], "NES": [1.0], "leadingEdge": ["GENE_A"]}).to_csv(expression / "ranked_pathways.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": ["GENE_A", "GENE_B"], "bin_id": ["bin-a", "bin-b"], "pvalue": [0.005, 0.4], "interaction_FDR": [0.02, 0.8]}).to_csv(exons / "dexseq_interaction.tsv", sep="\t", index=False)
    manifest = tmp_path / "project.json"
    manifest.write_text(json.dumps({"samples": [{"sample_id": "sample-a", "condition": "control", "biological_sample_id": "donor-a"}], "metadata": {"report_design_labels": {"16_Interaction/exon-design": "Explicit donor model"}, "report_design_notes": ["Source projects were declared explicitly."]}}), encoding="utf-8")
    before = {p: p.read_bytes() for p in base.rglob("*.tsv")}
    output = base / "24_Narrative_report"
    build_report("Example study", base, manifest, output)
    page = (output / "REPORT_EN.html").read_text(encoding="utf-8")
    status = json.loads((output / "report_status.json").read_text())
    assert "The 2 reported contrasts identify 3–7 significant AS events" in page
    assert "Explicit donor model" in page
    assert "Source projects were declared explicitly." in page
    assert 'lang="en"' in page
    assert not (output / "REPORT_ZH.html").exists()
    assert status["exon_designs"][0]["Design"] == "Explicit donor model"
    assert status["exon_designs"][0]["Interaction FDR <= 0.05"] == 1
    assert all(path.read_bytes() == content for path, content in before.items())


@pytest.mark.parametrize("metadata", [{"report_design_labels": []}, {"report_design_notes": [False]}])
def test_invalid_design_display_metadata_is_rejected(tmp_path, metadata):
    manifest = tmp_path / "project.json"
    manifest.write_text(json.dumps({"metadata": metadata}), encoding="utf-8")
    with pytest.raises(ValueError, match="Report design labels"):
        build_report("Example", tmp_path / "results", manifest, tmp_path / "output")


def test_input_mode_report_does_not_assert_a_fixed_cohort_or_metric_winner(tmp_path):
    from katana_cosplice.input_mode_benchmark import _write_report

    qc = pd.DataFrame({"matrix_source": ["BAM competing junctions"], "n_stable_events": [3]})
    similarity = pd.DataFrame({"matrix_source": ["BAM competing junctions"], "same_condition": [True], "spearman_correlation": [0.5]})
    bam_only = pd.DataFrame({"event_key": pd.Series(dtype=str), "passes_effect_cutoff": pd.Series(dtype=bool)})
    benchmark = pd.DataFrame({"method": ["katana_sparse", "wgcna_exact"], "silhouette_cosine": [0.1, 0.8]})
    path = _write_report(tmp_path, qc, similarity, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), bam_only, benchmark, {"sample_count": 6, "study_label": "Synthetic cohort"})
    text = path.read_text(encoding="utf-8")
    assert "Synthetic cohort" in text
    assert "12 samples" not in text and "36 samples" not in text
    assert "slightly stronger" not in text
    assert "Their relative performance is dataset-specific" in text
