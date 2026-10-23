"""Presentation fixture only; numerical native-MAEI execution is tested separately."""
from pathlib import Path

import pandas as pd
import pytest

from katana_next.flow.mapp_views import build_mapp_views
from katana_next.flow.model import file_sha256, write_json
from katana_next.flow.runtime import canonical_hash


def inputs(tmp_path):
    root = tmp_path / "collection"
    write_json(root / "project.json", {"project_id": "synthetic-view", "name": "Synthetic native-view fixture"})
    table = root / "products/native.tsv"
    table.parent.mkdir(parents=True)
    pd.DataFrame({"activities_A1": [-.2, .4], "activities_B1": [.3, -.1],
                  "combined.Zscore": [2., 3.], "corr_pval_A1": [.1, .02],
                  "corr_pval_B1": [.04, .2]}, index=["window|AAA", "window|CCC"]).to_csv(table, sep="\t")
    source = str(tmp_path / "native-task/AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv")
    task = {"id": "native-1", "adapter": "mapp.maei", "scope": {"contrast_ids": ["B_vs_A"], "site": "3ss"}}
    receipt = {"status": "completed", "execution_evidence": {"native_fitting_performed": True},
               "outputs": [{"path": source, "sha256": file_sha256(table)}]}
    receipt["receipt_sha256"] = canonical_hash(receipt)
    return root, {"tasks": [task]}, {"native-1": (None, receipt)}, [{"source": {"path": source}, "collection_path": "products/native.tsv"}], table


def test_native_maei_view_keeps_bonferroni_and_full_table_without_calibration_claim(tmp_path):
    root, plan, receipts, transport, table = inputs(tmp_path)
    original = file_sha256(table)
    result = build_mapp_views(root, plan, receipts, transport)
    view = result["extra_views"][0]
    assert len(view["panels"]) == 1 and view["panels"][0]["chapter"] == 7
    assert "not a treatment-effect test" in view["panels"][0]["caption"]
    assert view["tables"][0]["path"] == str(table)
    plotted = pd.read_csv(view["panels"][0]["data_links"][0]["path"], sep="\t")
    assert plotted.set_index("native_model_id").loc["window|AAA", "corr_pval_A1"] == .1
    assert file_sha256(table) == original
    assert all(Path(path).is_file() for path in result["files"])


def test_no_native_completion_cannot_be_presented_as_fitted_mapp(tmp_path):
    root, plan, receipts, transport, table = inputs(tmp_path)
    receipts["native-1"][1]["status"] = "running"
    with pytest.raises(ValueError, match="actual completed"):
        build_mapp_views(root, plan, receipts, transport)


def test_native_fit_table_tampering_is_rejected(tmp_path):
    root, plan, receipts, transport, table = inputs(tmp_path)
    table.write_text("changed", encoding="utf-8")
    with pytest.raises(ValueError, match="changed"):
        build_mapp_views(root, plan, receipts, transport)
