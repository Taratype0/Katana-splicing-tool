"""Readable scientific exports without changing plotted measurements."""
from katana_next.workflows.bulk import _plot


def readable_plot(*args, **kwargs):
    chart = _plot(*args, **kwargs)
    chart["background"] = "#ffffff"
    chart["config"]["font"] = "Georgia"
    chart["config"].setdefault("title", {})["font"] = "Georgia"
    for section in ("axis", "legend", "header"):
        chart["config"].setdefault(section, {}).update(labelFont="Georgia", titleFont="Georgia")
    return chart
