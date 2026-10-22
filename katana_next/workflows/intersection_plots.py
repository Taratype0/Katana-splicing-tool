from __future__ import annotations

from collections import Counter
from pathlib import Path
import textwrap

import pandas as pd
from matplotlib.figure import Figure

from katana_next.workflows.figure_exports import _RENDER_LOCK


def exact_intersections(sets: dict[str, set[str]]) -> pd.DataFrame:
    names = list(sets)
    membership = Counter(tuple(name for name in names if event in sets[name])
                         for event in set().union(*sets.values()))
    rows = [{"row_id": str(index), "included_comparisons": " | ".join(group),
             "excluded_comparisons": " | ".join(name for name in names if name not in group),
             "event_count": count, **{name: name in group for name in names}}
            for index, (group, count) in enumerate(sorted(membership.items(), key=lambda item: (-item[1], item[0])))]
    return pd.DataFrame(rows, columns=["row_id", "included_comparisons", "excluded_comparisons", "event_count", *names])


def render_intersections(sets: dict[str, set[str]], labels: dict[str, str], directory: Path, *, mode: str = "exact_event") -> dict[str, str]:
    """True Venn for 2-3 sets; exact UpSet membership for larger selections."""
    names = list(sets)
    directory.mkdir(parents=True, exist_ok=True)
    with _RENDER_LOCK:
        figure = Figure(figsize=(11, 7), layout="constrained")
        if mode == "event_and_direction":
            figure.suptitle("Same event + same dPSI sign required for overlap", fontsize=13)
        if len(names) <= 3:
            from matplotlib_venn import venn2, venn3
            axis = figure.subplots()
            renderer = venn2 if len(names) == 2 else venn3
            renderer([sets[name] for name in names],
                     set_labels=[textwrap.fill(labels.get(name, name), 42) + f"\nn = {len(sets[name]):,}" for name in names],
                     set_colors=("#0072B2", "#D55E00", "#009E73")[:len(names)], alpha=0.48, ax=axis)
            axis.set_title("Significant AS events shared between contrasts", fontsize=16, pad=25)
        else:
            table = exact_intersections(sets)
            # The complete membership table is exported even when the figure is capped.
            table = table.head(40)
            figure.set_size_inches(max(11, 0.36 * len(table) + 4), max(7, 0.38 * len(names) + 4))
            top, bottom = figure.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, max(1.3, len(names) * 0.25)]})
            top.bar(range(len(table)), table["event_count"], color="#0072B2")
            top.bar_label(top.containers[0], fontsize=8, padding=3)
            top.set_ylabel("Events in this exact intersection")
            top.set_title("Exact AS intersections across selected contrasts", fontsize=16)
            for x, row in enumerate(table.to_dict("records")):
                members = [i for i, name in enumerate(names) if row[name]]
                bottom.scatter([x] * len(names), range(len(names)), color="#dce1e6", s=22)
                bottom.scatter([x] * len(members), members, color="#0072B2", s=35)
                if members:
                    bottom.plot([x, x], [min(members), max(members)], color="#0072B2", lw=1.4)
            bottom.set_yticks(range(len(names)), [textwrap.fill(labels.get(name, name), 44) + f" (n={len(sets[name]):,})" for name in names])
            bottom.invert_yaxis()
            bottom.set_xticks([])
            bottom.set_xlabel("Each column is one exact membership pattern; top 40 shown, all exported to TSV.")
            for axis in (top, bottom):
                axis.spines[["top", "right"]].set_visible(False)
        exports = {}
        for extension in ("png", "svg"):
            path = directory / f"selected_contrast_intersections.{extension}"
            figure.savefig(path, dpi=170, bbox_inches="tight", facecolor="white")
            exports[extension] = str(path)
        figure.clear()
    return exports
def event_membership_sets(events, mode="exact_event"):
    """Direction splits a set member, not the genomic identity of an event."""
    import numpy as np
    import pandas as pd
    if mode not in {"exact_event", "event_and_direction"}:
        raise ValueError("Choose exact_event or event_and_direction membership")
    sets = {}
    rows = []
    for comparison, frame in events.groupby("comparison_id", sort=False):
        chosen = frame.loc[frame.stable_significant.fillna(False)].copy()
        chosen = chosen[np.isfinite(pd.to_numeric(chosen.dpsi, errors="coerce"))]
        if mode == "event_and_direction":
            chosen = chosen[chosen.dpsi.ne(0)]
        keys = set()
        for event in chosen.itertuples():
            direction = "positive" if event.dpsi > 0 else "negative" if event.dpsi < 0 else "zero"
            key = str(event.event_id) + (f"::direction={direction}" if mode == "event_and_direction" else "")
            keys.add(key)
            rows.append({"row_id": f"{comparison}::{key}", "comparison_id": comparison,
                         "event_id": event.event_id, "member_id": key, "direction": direction,
                         "dpsi": event.dpsi, "membership_mode": mode})
        sets[str(comparison)] = keys
    return sets, pd.DataFrame(rows, columns=["row_id", "comparison_id", "event_id", "member_id", "direction", "dpsi", "membership_mode"])
