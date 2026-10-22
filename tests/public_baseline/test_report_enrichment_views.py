"""Completed ORA rendering never silently computes or changes a test family."""
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from katana_next.reports.enrichment_views import CONTRACT, build_enrichment_views, enrichment_attachment_tables


def _rows():
    return pd.DataFrame({"comparison_id": ["B_vs_A"] * 3, "evidence": ["DAS"] * 3,
        "report_policy": ["stable_significant"] * 3, "library_kind": ["GO"] * 3,
        "ontology_category": ["GO BP"] * 3, "term": ["GO / t1", "GO / t2", "GO / t3"],
        "pvalue": [5e-310, .02, 1.], "fdr": [1e-309, .03, 1.],
        "fdr_family": ["this comparison/evidence/policy; all eligible GO terms"] * 3,
        "overlap_genes": [3, 2, 0], "selected_genes": [5] * 3, "background_genes": [100] * 3})


def test_rendering_retains_all_tests_and_never_calls_statistical_kernel(tmp_path, monkeypatch):
    from katana_next.science import enrichment
    def forbidden(*args, **kwargs):
        raise AssertionError("Report rendering must not calculate ORA or BH.")
    monkeypatch.setattr(enrichment, "overrepresentation", forbidden)
    monkeypatch.setattr(enrichment, "library_overrepresentation", forbidden)
    source = tmp_path / "ora.parquet"
    original = _rows()
    original.to_parquet(source)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    result = build_enrichment_views({}, [{"path": str(source), "metadata": {"report_enrichment_contract": CONTRACT}}], tmp_path / "report")
    assert result["new_scientific_fits"] == 0
    assert len(result["panels"]) == 1
    assert result["audit"][0]["eligible_tests"] == 3
    assert result["audit"][0]["displayed_terms"] == 2
    assert result["panels"][0]["policy"] == "stable_significant"
    assert (tmp_path / "report" / result["panels"][0]["path"]).read_bytes().startswith(b"\x89PNG")
    copied = pd.read_parquet(tmp_path / "report" / result["tables"][0]["path"])
    pd.testing.assert_frame_equal(copied, original)
    assert copied.iloc[0].pvalue == 5e-310
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest


def test_uncontracted_historical_results_are_not_relabelled(tmp_path):
    source = tmp_path / "historical.parquet"
    _rows().to_parquet(source)
    result = build_enrichment_views({}, [{"path": str(source), "metadata": {}}], tmp_path / "report")
    assert not result["panels"] and not result["tables"]


def test_mixed_family_in_one_scope_is_rejected(tmp_path):
    source = tmp_path / "ora.parquet"
    rows = _rows()
    rows.loc[1, "fdr_family"] = "different family"
    rows.to_parquet(source)
    with pytest.raises(ValueError, match="exactly one correction family"):
        build_enrichment_views({}, [{"path": str(source), "metadata": {"report_enrichment_contract": CONTRACT}}], tmp_path / "report")


def test_complete_foreground_background_and_provenance_files_are_explicit_attachments(tmp_path):
    foreground = tmp_path / "foreground.tsv"
    background = tmp_path / "background.tsv"
    receipt = tmp_path / "provenance.json"
    foreground.write_text("gene_id\ng1\n", encoding="utf-8")
    background.write_text("gene_id\ng1\ng2\n", encoding="utf-8")
    receipt.write_text('{"resource_sha256":"fixed"}', encoding="utf-8")
    (tmp_path / "unrelated.tsv").write_text("unrelated\nprivate\n", encoding="utf-8")
    declarations = [{"path": str(p), "name": name} for p, name in
                    [(foreground, "Complete foreground"), (background, "Detected background"), (receipt, "Resource provenance")]]
    tables = enrichment_attachment_tables({"report_enrichment_attachments": declarations})
    assert len(tables) == 3
    assert tables[0]["rows"] == 1 and tables[1]["rows"] == 2
    assert tables[2]["kind"] == "file" and "rows" not in tables[2]
    assert all(Path(row["path"]).is_absolute() for row in tables)
    assert all(row["chapter"] == 5 for row in tables)
    assert not any("unrelated" in row["path"] for row in tables)
    receipt.unlink()
    with pytest.raises(FileNotFoundError):
        enrichment_attachment_tables({"report_enrichment_attachments": declarations})
