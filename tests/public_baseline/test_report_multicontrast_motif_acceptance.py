"""Actual multi-comparison report artifacts from controlled generated inputs.

These controlled P/q records test display fidelity only. The separate generic
pipeline acceptance executes the native statistical kernel on generated RNA;
the controlled records here are never described as computed biological tests.
"""
from __future__ import annotations

import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd
from PIL import Image
import pytest

from katana_next.reports import build_report
from katana_next.reports.motif_views import SEQUENCE_FILES, STAT_FILES, build_completed_motif_views


PATTERNS = ("AAA", "CCC", "GGG", "UUU", "ACG", "CGU", "GUA", "UAC", "AAC", "CCA", "GGU", "UUA")


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")
    return Path(path)


def _controlled_family(root, cid, rotation):
    root.mkdir(parents=True)
    models = {f"m{i:02}": {"sequence": pattern, "RBPs": f"{cid}_FixtureRBP{i:02}"}
              for i, pattern in enumerate(PATTERNS)}
    _json(root / SEQUENCE_FILES[0], {"comparison_id": cid, "models": models, "permutations": 5,
          "synthetic_display_fixture": True, "provenance": {"statistical_family_split": False}})
    events, changes, forms, hits, tests = [], [], [], [], []
    for i, (token, model) in enumerate(models.items()):
        event = f"{cid}_E{i:02}"
        sequence = "CC" + model["sequence"] + "GG"
        events.append({"event_id": event, "event_type": "SE"})
        changes.append(dict(event_id=event, event_type="SE", model_id=token,
                            instance_state="inclusion_only", consequence_direction="increased"))
        forms.extend([dict(event_id=event, role="foreground", layer="mature_local_form", region=region,
                           sequence_RNA=word, length=len(word))
                      for region, word in (("inclusion", sequence), ("exclusion", "CCGG"))])
        hits.append(dict(event_id=event, model_id=token, form="inclusion", RNA_start0=2,
                         RNA_end0=5, hit_RNA=model["sequence"]))
        rank = (i + rotation) % 12
        q = [0.001, 0.004, 0.009, 0.015, 0.025, 0.04, 0.05, 0.051, 0.2, 0.5, 0.9, np.nan][rank]
        for start, p, adjusted, effect in ((0, q / 3, q, (i + 1) / 10),
                                           (7, 1e-9 if rank < 11 else np.nan,
                                            min(1., q + .08) if rank < 11 else np.nan, -(i + 1) / 10)):
            tests.append(dict(event_type="SE", direction="inclusion_up", model_id=token,
                              RBP=model["RBPs"], region="inclusion", method="controlled_fixture_not_computed",
                              window_nt=5, window_start0=start, n_matched_gene_pairs=24,
                              paired_rank_Z=i / 10, mean_paired_score_difference=effect,
                              raw_P_new_paired_swaps=p, q_new_window_family=adjusted,
                              status="estimable" if rank < 11 else "not_estimable"))
    for name, records in ((SEQUENCE_FILES[1], events), (SEQUENCE_FILES[2], changes),
                           (SEQUENCE_FILES[3], forms), (SEQUENCE_FILES[4], hits),
                           (SEQUENCE_FILES[5], [{"score": 1.25}]), (STAT_FILES[0], tests),
                           (STAT_FILES[1], [{"model_id": "m00", "cluster_P_FWER_new_family": .25}])):
        pd.DataFrame(records).to_csv(root / name, sep="\t", index=False, float_format="%.17g")
    np.save(root / STAT_FILES[2], np.array([0., 1., 1.5, 2., 2.5]))
    members = {name: {"sha256": _sha(root / name), "bytes": (root / name).stat().st_size}
               for name in (*SEQUENCE_FILES, *STAT_FILES)}
    marker = _json(root / "COMPLETE.json", {"complete": True, "statistical_family_split": False,
                                            "synthetic_display_fixture": True, "members": members})
    proof = _json(root.parent / (cid + ".fixture.receipt.json"), {
        "passed": True, "source_unchanged": True, "upstream_refits": 0,
        "scope": "controlled display fixture; no biological statistics claimed",
        "summary": {"comparison_id": cid}, "complete_manifest": {"sha256": _sha(marker)}})
    return root, proof


class _Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.figures, self.links, self.sections = [], [], []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "figure":
            self.figures.append(attrs)
        if tag in {"a", "img", "script", "link"}:
            href = attrs.get("href", attrs.get("src"))
            if href:
                self.links.append(href)
        if tag == "section":
            self.sections.append(attrs.get("id"))


def test_multicontrast_top10_report_has_complete_portable_statistics_and_figures(tmp_path):
    ids = ["B_vs_A", "C_vs_A", "C_vs_B"]
    project = {"project_id": "controlled-motif-report-fixture", "name": "Generated display acceptance fixture",
               "metadata": {"synthetic": True}, "annotation": {"genome_build": "fixture-only"},
               "samples": [], "contrasts": [{"contrast_id": cid,
                    "numerator": {"label": cid.split("_vs_")[0], "sample_ids": []},
                    "denominator": {"label": cid.split("_vs_")[1], "sample_ids": []}} for cid in ids]}
    manifest = _json(tmp_path / "project.json", project)
    views, source_identities = [], {}
    for index, cid in enumerate(ids):
        source, proof = _controlled_family(tmp_path / "sources" / cid, cid, index * 3)
        source_identities.update({str(p): _sha(p) for p in source.iterdir()})
        views.append(build_completed_motif_views(source, tmp_path / "views" / cid,
            comparison_id=cid, project_id=project["project_id"],
            source_project_manifest_sha256=_sha(manifest), completion_proof=proof, proof_kind="verified_local"))
    out = tmp_path / "report"
    report = build_report({"project_manifest": str(manifest), "artifacts": [], "extra_views": views,
                           "dummy": True, "report_require_plot_data": True}, out)
    page = _Page()
    page.feed((out / "REPORT.html").read_text(encoding="utf-8"))
    assert all(f"chapter-{chapter}" in page.sections for chapter in range(1, 8))
    assert report["event_set"] == "stable_significant"
    assert report["contrasts"] == ids
    assert not report["figure_data_links"]["unbound_main_figures"]
    card_counts = {}
    statistic_checks = []
    for view, cid in zip(views, ids, strict=True):
        assert all(table["contrast_id"] == cid and table["name"].endswith(" | " + cid) for table in view["tables"])
        cards = [p for p in report["panels"] if p["contrast_id"] == cid and p["title"].startswith("Top motif")]
        assert len(cards) == 10
        card_counts[cid] = len(cards)
        assert view["audit"]["top_motif_RNA_examples"] == 10
        ranking = pd.read_csv(tmp_path / "views" / cid / "all-motif-ranked-recorded-tests.tsv", sep="\t", float_precision="round_trip")
        native = pd.read_csv(tmp_path / "sources" / cid / STAT_FILES[0], sep="\t", float_precision="round_trip")
        assert len(ranking) == 12 and ranking.model_id.nunique() == 12
        assert ranking["comparison_id"].eq(cid).all()
        assert ranking["significance"].tolist().count("recorded_test_FDR_le_0.05") == 7
        assert ranking.iloc[-1].significance == "not_estimable"
        for row in ranking.iloc[:-1].itertuples():
            source_row = native.loc[native.model_id.eq(row.model_id) & native.window_start0.eq(row.window_start0)].iloc[0]
            assert row.window_start0 == 0
            for field in ("raw_P_new_paired_swaps", "q_new_window_family", "mean_paired_score_difference"):
                assert getattr(row, field) == source_row[field]
        for panel in cards:
            assert panel["metadata"]["das_scope"]["recorded_event_set"] == "stable_significant"
            assert panel["policy"] is None  # completed RS motif scope is explicit, never relabelled by UI
            assert all(link["contrast_id"] == cid for link in panel["data_links"])
            assert "P =" in panel["caption"] and "FDR =" in panel["caption"]
            assert "not a whole-model hypothesis" in panel["caption"]
            with Image.open(out / panel["image"]) as image:
                pixels = np.asarray(image.convert("RGB"))
                assert image.width > 600 and image.height > 300 and pixels.std() > 5
            ET.parse(out / panel["svg"])
        assert any("not significant" in p["caption"] for p in cards)
        main = [p for p in report["panels"] if p["contrast_id"] == cid and not p.get("output_only")]
        assert not any(any(term in p["title"] for term in ("window profile", "maximum-cluster null", "sequence opportunity", "All-model completed")) for p in main)
        statistic_checks.append({"comparison_id": cid, "models_retained": 12,
                                "top10": 10, "paired_native_values_preserved": True})
    assert len({p["id"] for p in report["panels"]}) == len(report["panels"])
    assert len([x for x in page.figures if "top10-motif-" in x.get("id", "")]) == 30
    for link in page.links:
        if "://" not in link and not link.startswith("#"):
            assert (out / link.split("#", 1)[0]).is_file(), link
    with zipfile.ZipFile(out / "REPORT.zip") as archive:
        assert archive.testzip() is None
        names = set(archive.namelist())
        for path in out.rglob("*"):
            if path.is_file() and path.name != "REPORT.zip":
                relative = path.relative_to(out).as_posix()
                assert relative in names and archive.read(relative) == path.read_bytes()
    assert source_identities == {name: _sha(name) for name in source_identities}
    _json(tmp_path / "report-acceptance.receipt.json", {
        "passed": True, "scope": "controlled generated display inputs; not a biological statistical validation",
        "comparisons": statistic_checks, "motif_cards": card_counts, "main_top10_figures": 30,
        "native_P_q_effect_same_row": True, "source_files_unchanged": True,
        "all_report_links_exist": True, "all_png_and_svg_decoded": True,
        "seven_chapters_preserved": True, "complete_zip_members_byte_identical": True,
        "ordinary_selector_does_not_relabel_fixed_RS_motif": True})


def test_degenerate_consensus_remains_visible_without_fabricating_pwm_probabilities():
    from matplotlib.figure import Figure
    from katana_next.reports.motif_views import _logo, _matrix, _motif_ranking
    model = {"kind": "consensus", "pattern": "[CT]CT[CT]", "RBPs": "FixtureRBP"}
    assert _matrix(model) is None
    ranking = _motif_ranking(pd.DataFrame([{"model_id": "degenerate"}]), {}, {"degenerate": model})
    assert ranking.iloc[0].motif_sequence == "[CU]CU[CU]"
    assert ranking.iloc[0].sequence_representation == "recorded IUPAC/fixed-pattern consensus"
    assert ranking.iloc[0].significance == "not_estimable"
    fig = Figure()
    ax = fig.subplots()
    _logo(ax, model)
    assert any("[CU]CU[CU]" in text.get_text() for text in ax.texts)
    assert any("no base probabilities supplied" in text.get_text() for text in ax.texts)


@pytest.mark.parametrize("column", ["raw_P_new_paired_swaps", "q_new_window_family"])
@pytest.mark.parametrize("invalid", [-.01, 1.01, np.inf])
def test_invalid_recorded_significance_cannot_be_presented_as_a_valid_top_motif(column, invalid):
    from katana_next.reports.motif_views import _best_recorded_tests
    row = {"raw_P_new_paired_swaps": .01, "q_new_window_family": .05}
    row[column] = invalid
    with pytest.raises(ValueError, match="finite in"):
        _best_recorded_tests(pd.DataFrame([row]), {})
