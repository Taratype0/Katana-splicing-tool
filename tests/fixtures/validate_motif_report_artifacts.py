"""Read-only acceptance of actual generated pipeline MotifTrace/report outputs.

This does not recompute any sequence or statistic. Pass the finished report
directory produced by the generic Nextflow acceptance, not an existing study.
"""
from __future__ import annotations

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd
from PIL import Image


class _Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.figures = [], []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "figure":
            self.figures.append(attributes)
        if tag in {"a", "img", "script", "link"}:
            link = attributes.get("href", attributes.get("src"))
            if link:
                self.links.append(link)


def validate_motif_report(directory, *, model_count=12, comparison_count=2):
    root = Path(directory).resolve()
    manifest = json.loads((root / "report-manifest.json").read_text(encoding="utf-8"))
    provenance = json.loads((root / "provenance.json").read_text(encoding="utf-8"))
    page = (root / "REPORT.html").read_text(encoding="utf-8")
    assert len(manifest["contrasts"]) == comparison_count
    assert manifest["event_set"] == "stable_significant"
    assert all(f'id="chapter-{chapter}"' in page for chapter in range(1, 8))
    assert not manifest["figure_data_links"]["unbound_main_figures"]
    summary = []
    for comparison in manifest["contrasts"]:
        audits = [r for r in provenance["derived_views"]
                  if isinstance(r, dict) and r.get("schema", "").startswith("completed-motif-report") and r.get("comparison_id") == comparison]
        assert len(audits) == 1
        audit = audits[0]
        assert audit["proof_kind"] == "flow" and audit["models"] == model_count
        assert audit["new_scans"] == audit["new_statistical_tests"] == 0
        assert audit["native_P_q_FWER_unchanged"] and audit["originals_unchanged"]
        cards = [p for p in manifest["panels"] if p["contrast_id"] == comparison and p["title"].startswith("Top motif")]
        assert len(cards) == min(model_count, 10)
        full_link = next(link for link in cards[0]["data_links"] if link["label"] == "All ranked models and native test statistics")
        native_link = next(link for link in cards[0]["data_links"] if link["label"] == "Complete native test family")
        ranking = pd.read_csv(root / full_link["path"], sep="\t", float_precision="round_trip")
        native = pd.read_csv(root / native_link["path"], sep="\t", float_precision="round_trip")
        assert ranking.model_id.nunique() == model_count
        assert ranking.comparison_id.eq(comparison).all()
        estimable = native.loc[native.raw_P_new_paired_swaps.notna() & native.q_new_window_family.notna()]
        # This fixture deliberately includes 24 independently matched SE genes;
        # a completely NA family would be an execution/fixture regression.
        assert not estimable.empty
        assert estimable.n_matched_gene_pairs.ge(20).all()
        assert estimable.raw_P_new_paired_swaps.between(0, 1).all()
        assert estimable.q_new_window_family.between(0, 1).all()
        columns = ["event_type", "direction", "model_id", "region", "method", "window_nt", "window_start0"]
        for row in ranking.to_dict("records"):
            if pd.isna(row["q_new_window_family"]):
                assert row["significance"] == "not_estimable"
                continue
            matched = native.copy()
            for column in columns:
                matched = matched.loc[matched[column].eq(row[column])]
            assert len(matched) == 1
            for column in ("raw_P_new_paired_swaps", "q_new_window_family", "mean_paired_score_difference"):
                assert row[column] == matched.iloc[0][column]
            expected = "recorded_test_FDR_le_0.05" if row["q_new_window_family"] <= .05 else "not_significant"
            assert row["significance"] == expected
        ordered = ranking.sort_values(["q_new_window_family", "raw_P_new_paired_swaps", "model_id"], kind="stable", na_position="last")
        assert ordered.model_id.tolist() == ranking.model_id.tolist()
        for card in cards:
            assert card["policy"] is None and card["metadata"]["das_scope"]["recorded_event_set"] == "stable_significant"
            assert "P =" in card["caption"] and "FDR =" in card["caption"]
            assert all(link["contrast_id"] == comparison for link in card["data_links"])
            assert all((root / link["path"]).is_file() for link in card["data_links"])
            with Image.open(root / card["image"]) as image:
                assert image.width > 600 and image.height > 300
                assert np.asarray(image.convert("RGB")).std() > 5
            ET.parse(root / card["svg"])
        summary.append({"comparison_id": comparison, "declared_models": model_count, "main_motif_cards": len(cards),
                        "all_test_rows": len(native), "estimable_test_rows": len(estimable),
                        "selected_recorded_tests_FDR_le_005": int(ranking.q_new_window_family.le(.05).sum()),
                        "native_P_q_effect_same_row": True, "test_family_not_model_FDR": True,
                        "fixed_RS_scope": True, "all_models_and_native_test_rows_retained": True})
    for page_path in (root / "REPORT.html", root / "outputs/index.html"):
        parser = _Links()
        parser.feed(page_path.read_text(encoding="utf-8"))
        for link in parser.links:
            if "://" not in link and not link.startswith("#"):
                assert (page_path.parent / link.split("#", 1)[0]).is_file(), (page_path, link)
    with zipfile.ZipFile(root / "REPORT.zip") as archive:
        assert archive.testzip() is None
        for panel in manifest["panels"]:
            if panel["title"].startswith("Top motif"):
                assert archive.read(panel["image"]) == (root / panel["image"]).read_bytes()
    return {"passed": True, "scope": "actual synthetic pipeline output and report; no study or public data",
            "comparisons": summary, "seven_chapters_preserved": True, "links_exist": True,
            "top10_png_svg_decoded": True, "zip_top10_images_byte_identical": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("report")
    parser.add_argument("--receipt", required=True)
    args = parser.parse_args()
    receipt = Path(args.receipt)
    if receipt.exists():
        raise FileExistsError("Do not overwrite an earlier acceptance receipt")
    result = validate_motif_report(args.report)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
