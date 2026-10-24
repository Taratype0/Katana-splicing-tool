"""Fresh completed-result views retain the accepted graph and padding rules."""
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from katana_next.flow.presentation import _producer_path
from katana_next.reports.network_views import (
    _authenticated_display_sources, _pad_enrichment_markers, _scatter_marker_bounds,
)


def _identity_fixture(root):
    native = root / "products/fit/compute/engine/04_network/observed"
    native.mkdir(parents=True)
    assignments = native / "event_module_assignments.tsv"
    assignments.write_text("event_key\tmodule\ne1\tM01\n")
    (native / "parameters.json").write_text('{"settings":{"signed":true}}')
    (native / "network_diagnostics.tsv").write_text("method\nbicor\n")
    original = "/cfs/owned/one-fit/compute/engine/04_network/observed"
    files = [native / "parameters.json", native / "network_diagnostics.tsv"]
    receipt = {"task_id": "fit", "outputs": [{"path": original + "/" + p.name,
        "size_bytes": p.stat().st_size, "sha256": sha256(p.read_bytes()).hexdigest()} for p in files]}
    receipt_path = root / "evidence/task.receipt.json"
    receipt_path.parent.mkdir()
    receipt_path.write_text(json.dumps(receipt))
    proof = {"files": [{"path": assignments.relative_to(root).as_posix(), "task_id": "fit",
                       "source_path": original + "/event_module_assignments.tsv"}],
             "task_receipts": [{"path": receipt_path.relative_to(root).as_posix(),
                                "sha256": sha256(receipt_path.read_bytes()).hexdigest()}]}
    return root / "presentation/view-inputs.json", proof, assignments, files


def test_missing_old_proof_members_can_use_only_same_fit_output_hashes(tmp_path):
    completion, proof, assignments, files = _identity_fixture(tmp_path)
    before = {p: p.read_bytes() for p in files}
    result = _authenticated_display_sources(completion, proof, {"provenance": {"source_task_id": "fit"}}, [assignments])
    assert result == files
    assert before == {p: p.read_bytes() for p in files}
    files[0].write_text('{"settings":{"signed":false}}')
    with pytest.raises(ValueError, match="exact completed fit output ledger"):
        _authenticated_display_sources(completion, proof, {"provenance": {"source_task_id": "fit"}}, [assignments])


def test_similar_suffix_from_another_fit_cannot_supply_display_parameters(tmp_path):
    completion, proof, assignments, _ = _identity_fixture(tmp_path)
    proof["files"][0]["source_path"] = "/another-fit/compute/engine/04_network/observed/event_module_assignments.tsv"
    with pytest.raises(ValueError, match="exact completed fit output ledger"):
        _authenticated_display_sources(completion, proof, {"provenance": {"source_task_id": "fit"}}, [assignments])


def test_producer_paths_keep_their_own_os_identity_without_consumer_resolution():
    assert _producer_path("/cfs/fit/04_network/observed/file.tsv") == PurePosixPath("/cfs/fit/04_network/observed/file.tsv")
    assert _producer_path(r"C:\Study\fit\file.tsv") == PureWindowsPath("c:/study/fit/file.tsv")
    assert _producer_path(r"\\cluster\share\fit\file.tsv") == PureWindowsPath("//cluster/share/fit/file.tsv")
    assert _producer_path("/cfs/fit/file.tsv") != _producer_path("E:/cfs/fit/file.tsv")
    assert _producer_path("/other/fit/file.tsv") != _producer_path("/cfs/fit/file.tsv")


def _svg_markers_inside_clip(svg, expected):
    root = ET.parse(svg).getroot()
    collection = next(e for e in root.iter() if e.get("id") == "PathCollection_1")
    parents = {child: parent for parent in root.iter() for child in parent}
    markers = [e for e in collection.iter() if (e.tag.endswith("path") and not e.get("id")) or e.tag.endswith("use")]
    assert len(markers) == expected
    for element in markers:
        source, x, y = element, 0., 0.
        if element.tag.endswith("use"):
            marker_id = element.get("{http://www.w3.org/1999/xlink}href").lstrip("#")
            source = next(e for e in root.iter() if e.get("id") == marker_id)
            x, y = float(element.get("x", 0)), float(element.get("y", 0))
        values = [float(value) for value in re.findall(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?", source.get("d"))]
        style = dict(item.strip().split(":", 1) for item in element.get("style", "").split(";") if ":" in item)
        stroke = float(style.get("stroke-width", 1)) / 2 if "stroke" in style else 0
        bounds = [min(values[::2]) + x - stroke, min(values[1::2]) + y - stroke,
                  max(values[::2]) + x + stroke, max(values[1::2]) + y + stroke]
        clip = element
        while not clip.get("clip-path"):
            clip = parents[clip]
        clip_id = re.fullmatch(r"url\(#([^)]*)\)", clip.get("clip-path")).group(1)
        rectangle = next(e for e in root.iter() if e.get("id") == clip_id)[0]
        left, top, width, height = [float(rectangle.get(key)) for key in ("x", "y", "width", "height")]
        assert min(bounds[0] - left, bounds[1] - top, left + width - bounds[2], top + height - bounds[3]) >= -1e-5


def test_accepted_marker_padding_keeps_full_bubbles_in_png_and_svg_without_changing_values(tmp_path):
    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained", dpi=150)
    x, areas = np.array([1.3, 1.301, 1.302]), np.array([2000., 8000., 4000.])
    points = ax.scatter(x, range(3), s=areas, color="#2863A5", alpha=.8)
    labels = ["M01 | GO BP | exact original term " + str(i) for i in range(3)]
    ax.set(yticks=range(3), yticklabels=labels, xlabel="−log10(native BH FDR)")
    ax.tick_params(axis="y", labelsize=6)
    assert _scatter_marker_bounds(fig, ax, points)[2] < 0  # Reproduce the clipping defect.
    receipt = _pad_enrichment_markers(fig, ax, points)
    assert receipt["all_point_extents_inside_axes"] and receipt["minimum_margin_px"] > 4
    np.testing.assert_array_equal(points.get_offsets()[:, 0], x)
    np.testing.assert_array_equal(points.get_offsets()[:, 1], np.arange(3))
    np.testing.assert_array_equal(points.get_sizes(), areas)
    assert [tick.get_text() for tick in ax.get_yticklabels()] == labels
    png, svg = tmp_path / "points.png", tmp_path / "points.svg"
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    assert png.stat().st_size > 1000 and _scatter_marker_bounds(fig, ax, points)[2] > 0
    _svg_markers_inside_clip(svg, 3)
    plt.close(fig)


def test_no_significant_terms_retains_empty_display_without_added_points():
    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained", dpi=150)
    points = ax.scatter([], [], s=[], color="#2863A5", alpha=.8)
    ax.text(.5, .5, "No module GO / KEGG terms meet native FDR <= 0.05", transform=ax.transAxes)
    ax.axis("off")
    receipt = _pad_enrichment_markers(fig, ax, points)
    assert receipt["point_count"] == 0 and receipt["minimum_margin_px"] is None
    assert receipt["all_point_extents_inside_axes"] and receipt["scientific_values_recomputed"] is False
    plt.close(fig)
