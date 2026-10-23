"""Comparison-specific displays of completed molecular results.

This module only reads completed tables. It does not fit models, change native
P/FDR/effects, open a project database, or reuse historical RS statistics under
a new label. Both report foregrounds are derived explicitly from original
counts; historical flags remain in the complete exported table.
"""
from __future__ import annotations

from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .sources import read_table, table_schema, report_contrast_selection

BLUE, RED, GRAY = "#2863A5", "#C63F4A", "#CDD1D6"
EVENT_COLORS = [BLUE, RED, "#218D87", "#9A64B0", "#D49A25"]
EVENTS = ("SE", "A3SS", "A5SS", "MXE", "RI")
POLICIES = ("stable_significant", "significant")
LABELS = {"stable_significant": "RS-DAS", "significant": "DAS"}


def _dict(value):
    return value.model_dump(mode="json") if hasattr(value, "model_dump") else dict(value)


def _safe(value):
    original = str(value)
    clean = re.sub(r"[^A-Za-z0-9._-]+", "-", original).strip("-")[:140] or "result"
    # Distinct valid scientific comparison IDs can normalize to the same HTML
    # slug (e.g. 'B:A' and 'B A'). Preserve simple established IDs, hashing only
    # when normalization/truncation changed their spelling.
    return clean if clean == original else clean + "-" + hashlib.sha256(original.encode()).hexdigest()[:10]


def _num(frame, key):
    return pd.to_numeric(frame[key], errors="coerce") if key in frame else pd.Series(np.nan, index=frame.index)


def _genes(frame):
    return set(frame.get("gene_id", pd.Series(dtype=str)).dropna().astype(str)) - {"", "nan", "None"}


def _rank(frame, fdr, effect=None, identity="row_id"):
    """Stable ordering; missing effects are not fabricated for usage results."""
    result = frame.copy()
    result["_report_fdr"] = _num(result, fdr)
    result["_report_abs_effect"] = _num(result, effect).abs() if effect else np.nan
    if identity not in result:
        result[identity] = result.index.astype(str)
    result["_report_identity"] = result[identity].astype(str)
    return result.sort_values(["_report_fdr", "_report_abs_effect", "_report_identity"],
                              ascending=[True, False, True], na_position="last", kind="stable").drop(
                                  columns=["_report_fdr", "_report_abs_effect", "_report_identity"])


def _counts(value, expected):
    try:
        values = np.asarray([float(item.strip()) for item in str(value).split(",")], dtype=float)
    except (TypeError, ValueError):
        return None
    if len(values) != expected or not np.isfinite(values).all() or (values < 0).any():
        return None
    return values


def _rs_eligibility(frame, contrast):
    """Actual count arrays in both declared arms, never a mean/old stable flag."""
    n1, n2 = len(contrast["numerator"]["sample_ids"]), len(contrast["denominator"]["sample_ids"])
    required = ("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2")
    if not all(column in frame for column in required) or min(n1, n2) == 0:
        return pd.Series(False, index=frame.index), "missing_original_count_arrays_or_declared_samples"
    valid = []
    for values in frame[list(required)].itertuples(index=False, name=None):
        arrays = [_counts(values[i], n1 if i < 2 else n2) for i in range(4)]
        valid.append(all(x is not None for x in arrays)
                     and bool(((arrays[0] + arrays[1]) >= 10).all())
                     and bool(((arrays[2] + arrays[3]) >= 10).all()))
    usable = _num(frame, "fdr").between(0, 1) & _num(frame, "dpsi").between(-1, 1)
    for psi_key, expected in (("IncLevel1", n1), ("IncLevel2", n2)):
        if psi_key in frame:
            usable &= frame[psi_key].map(lambda value: (values := _counts(value, expected)) is not None
                                         and bool((values <= 1).all()))
    return pd.Series(valid, index=frame.index) & usable, "original_counts_all_declared_samples_ge10"


def _native_significance(frame, fdr, effect=None, *, fdr_cutoff=.05, effect_cutoff=0):
    if "significant" in frame:
        observed = frame["significant"]
        if not observed.dropna().isin([True, False, 0, 1]).all():
            raise ValueError("Completed significant flags must be explicit booleans/0/1.")
        return observed.fillna(False).astype(bool)
    selected = _num(frame, fdr).le(fdr_cutoff)
    if effect:
        selected &= _num(frame, effect).abs().ge(effect_cutoff)
    return selected


class _Views:
    def __init__(self, project, artifacts, destination):
        self.project = project
        self.artifacts = artifacts
        self.root = Path(destination).resolve()
        self.panels, self.tables, self.top10, self.audit = [], [], [], []
        self.gene_sets, self.policy_data = [], {key: {} for key in POLICIES}
        self.snapshots = {}

    def input(self, path):
        path = Path(path)
        if self.root == path.resolve() or self.root in path.resolve().parents:
            raise ValueError("Report destination must not contain a source input.")
        stat = path.stat()
        self.snapshots[path] = (stat.st_size, stat.st_mtime_ns)
        return path

    def source(self, analysis, name):
        found = [item for item in self.artifacts if item.get("analysis_id") == analysis and item.get("name") == name
                 and Path(item.get("path", "")).is_file()]
        if not found:
            self.audit.append({"analysis": analysis, "source": name, "status": "not_computed"})
            return None
        record = found[0]
        self.input(record["path"])
        return record

    def frame(self, artifact, comparison):
        if artifact is None:
            return pd.DataFrame()
        schema = table_schema(artifact["path"])
        if "comparison_id" in schema["columns"]:
            if Path(artifact["path"]).suffix == ".parquet":
                return read_table(artifact["path"], filters=[("comparison_id", "=", comparison)])
            frame = read_table(artifact["path"])
            return frame.loc[frame.comparison_id.eq(comparison)].copy()
        native_ids = self.project.get('report_native_contrast_ids', [row['contrast_id'] for row in self.project.get('contrasts', [])])
        if artifact.get("contrast_ids") == [comparison] or len(native_ids) == 1:
            return read_table(artifact["path"]).assign(comparison_id=comparison)
        self.audit.append({"source": artifact.get("name"), "comparison_id": comparison,
                           "status": "unavailable", "reason": "No authenticated comparison column/scope."})
        return pd.DataFrame()

    def table(self, frame, name, comparison, unit, *, policy=None, suffix="parquet", chapter=3):
        directory = self.root / "outputs/comparisons" / _safe(comparison)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (_safe(name) + "." + suffix)
        if suffix == "parquet":
            frame.to_parquet(path, index=False)
        else:
            frame.to_csv(path, sep="\t", index=False)
        record = {"id": _safe(comparison + "-" + name), "name": name, "title": name,
                  "path": path.relative_to(self.root).as_posix(), "source": path.relative_to(self.root).as_posix(),
                  "contrast_id": comparison, "comparison_id": comparison, "rows": len(frame),
                  "columns": list(frame), "unit": unit, "policy": policy, "chapter": chapter}
        self.tables.append(record)
        return record["path"]

    def figure(self, fig, identity, chapter, title, comparison, caption, policy=None):
        directory = self.root / "figures/comparisons"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (_safe(identity) + ".png")
        fig.savefig(path, dpi=145, bbox_inches="tight", facecolor="white")
        fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.clear()
        self.panels.append({"id": _safe(identity), "chapter": chapter, "title": title,
                            "path": path.relative_to(self.root).as_posix(),
                            "caption": caption, "contrast_id": comparison or "all", "scope": "comparison" if comparison else "global",
                            "policy": policy, "event_set": policy})

    def bind_data(self, plotted, unit, *, full_results=None, comparison_ids=None):
        """Attach each panel's exact display values, independently of caption."""
        panel = self.panels[-1]
        comparison, policy = panel.get('contrast_id') or 'all', panel.get('policy')
        name = panel['id'] + '-plot-data'
        path = self.table(plotted, name, comparison, unit, policy=policy,
                          suffix='tsv', chapter=panel['chapter'])
        links = [{'role': 'plot_data', 'label': 'Exact plotted data', 'path': path,
                  'contrast_id': comparison, 'policy': policy, 'unit': unit}]
        for label, source, source_unit in full_results or []:
            links.append({'role': 'full_results', 'label': label, 'path': source,
                          'contrast_id': comparison, 'policy': policy, 'unit': source_unit})
        if comparison_ids is not None:
            panel['comparison_ids'] = list(comparison_ids)
            for link in links:
                link['comparison_ids'] = list(comparison_ids)
        panel['data_links'] = links
        panel['unit'] = unit
        return path

    def top(self, frame, analysis, comparison, fdr, effect, unit, *, policy=None, full_table=None):
        selected = frame.loc[frame["report_significant"]].copy() if "report_significant" in frame else frame
        ranked = _rank(selected, fdr, effect).head(10)
        prefix = LABELS.get(policy, analysis) if analysis == "DAS" else analysis
        self.top10.append({"id": _safe("top10-" + comparison + "-" + analysis + "-" + str(policy or "common")),
                           "analysis": analysis, "title": f"Top {len(ranked)} {prefix} | {comparison}",
                           "contrast_id": comparison, "comparison_id": comparison, "kind": analysis,
                           "policy": policy, "event_set": policy,
                           "unit": unit, "columns": list(ranked),
                           "rows": ranked.astype(object).where(ranked.notna(), None).to_dict(orient="records"),
                           "full_table": full_table, "full_path": full_table, "total_rows": len(selected),
                           "caption": "Native FDR ascending, absolute available effect descending, then stable object ID. Complete result is linked; original IDs and caller statistics are retained.",
                           "ranking": "native FDR ascending; absolute available effect descending; stable object ID ascending"})
        return ranked

    def gene_set(self, frame, foreground, background, analysis, comparison, policy=None):
        if not _genes(frame):
            return
        names = frame[[c for c in ("gene_id", "gene_symbol") if c in frame]].drop_duplicates("gene_id")
        fore = names.loc[names.gene_id.isin(foreground)]
        back = names.loc[names.gene_id.isin(background)]
        if not set(fore.gene_id).issubset(set(back.gene_id)):
            raise ValueError("Enrichment foreground is not a subset of its measured/testable background.")
        stem = analysis + "-" + str(policy or "common")
        self.gene_sets.append({"analysis": analysis, "contrast_id": comparison, "policy": policy,
                              "unit": "unique_gene", "foreground_genes": len(fore), "background_genes": len(back),
                              "foreground": self.table(fore, stem + "-foreground-genes", comparison, "unique_gene", policy=policy, suffix="tsv", chapter=5),
                              "background": self.table(back, stem + "-detected-background", comparison, "unique_gene", policy=policy, suffix="tsv", chapter=5),
                              "background_scope": "detected/testable genes in this completed comparison and source layer; never all annotated genes"})


def _axes(title, xlabel="", ylabel=""):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    fig = Figure(figsize=(9, 5.5))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#ECEEF1", linewidth=.6, zorder=0)
    return fig, ax


def _volcano(view, frame, comparison, analysis, effect, fdr, selected, *, policy=None, xlabel=None, fdr_cutoff=.05, full_results=None):
    title = f"{LABELS.get(policy, analysis) if analysis == 'DAS' else analysis} | {comparison}"
    fig, ax = _axes(title, xlabel or effect, "−log10(native FDR)")
    x, q = _num(frame, effect), _num(frame, fdr)
    valid = x.notna() & np.isfinite(x) & q.between(0, 1)
    # Zero P/FDR is placed at a finite plotting floor only; the saved table and
    # Top 10 retain literal zero and native subnormal positive values.
    y = -np.log10(q.clip(lower=np.nextafter(0., 1.)))
    colors = np.where(selected & x.ge(0), RED, np.where(selected & x.lt(0), BLUE, GRAY))
    ax.scatter(x[valid], y[valid], c=colors[valid], s=8, alpha=.5, linewidths=0, rasterized=True)
    if "fdr_cutoff" in frame and _num(frame, "fdr_cutoff").dropna().nunique() == 1:
        fdr_cutoff = float(_num(frame, "fdr_cutoff").dropna().iloc[0])
    ax.axhline(-np.log10(fdr_cutoff), color="#85909B", linewidth=.7, linestyle="--")
    for _, row in _rank(frame.loc[valid & selected], fdr, effect).head(8).iterrows():
        label = row.get("gene_symbol") or row.get("gene_id") or row.get("event_id")
        ax.annotate(str(label), (row[effect], -np.log10(max(float(row[fdr]), np.nextafter(0., 1.)))),
                    fontsize=7, xytext=(3, 4), textcoords="offset points")
    view.figure(fig, comparison + "-" + analysis + "-volcano-" + str(policy or "common"), 3, title,
                comparison, "Original caller effect and FDR; red positive, blue negative. All finite tested rows are plotted; only foreground coloring changes. Plotting floor for native zero FDR does not change saved values.", policy)
    view.panels[-1]['plot_parameters'] = {'effect_column': effect, 'fdr_column': fdr,
        'fdr_reference_line': fdr_cutoff, 'zero_fdr_display_floor': float(np.nextafter(0., 1.))}
    plotted = frame.loc[valid].copy()
    plotted['plot_x'] = x.loc[valid]
    plotted['plot_y_minus_log10_fdr'] = y.loc[valid]
    plotted['plot_foreground'] = selected.loc[valid].astype(bool)
    plotted['plot_color'] = colors[valid.to_numpy()]
    view.bind_data(plotted, 'gene' if analysis == 'DEG' else ('exact_splicing_event' if analysis == 'DAS' else 'exon_bin'),
                   full_results=[('Complete native results', full_results, 'native_result_row')] if full_results else [])


def _heat(view, data, rowlabels, collabels, comparison, name, caption, *, chapter=3, policy=None, bounded=False, full_results=None, row_ids=None):
    from matplotlib.colors import LinearSegmentedColormap
    if data.size == 0:
        view.audit.append({"comparison_id": comparison, "figure": name, "status": "not_computed", "reason": "No eligible matrix rows."})
        return
    fig, ax = _axes(name + " | " + comparison)
    fig.set_size_inches(max(8, len(collabels) * .4), max(4, min(14, len(rowlabels) * .28 + 1.8)))
    values = np.asarray(data, dtype=float)
    limit = 1. if bounded else max(1., float(np.nanpercentile(np.abs(values), 98)))
    image = ax.imshow(values, aspect="auto", cmap=LinearSegmentedColormap.from_list("katana_blue_red", [BLUE, "white", RED]),
                      vmin=-limit if not bounded else 0, vmax=limit)
    ax.set_xticks(range(len(collabels)), collabels, rotation=55, ha="right", fontsize=7)
    ax.set_yticks(range(len(rowlabels)), rowlabels, fontsize=8)
    ax.grid(False)
    fig.colorbar(image, ax=ax, shrink=.7)
    view.figure(fig, comparison + "-" + name + "-" + str(policy or "common"), chapter,
                name + " | " + comparison, comparison, caption, policy)
    view.panels[-1]['plot_parameters'] = {'color_limits': [0 if bounded else -limit, limit], 'missing_values': 'unchanged'}
    plotted = pd.DataFrame(values, columns=list(collabels))
    plotted.insert(0, 'display_row', list(rowlabels))
    plotted.insert(0, 'display_row_index', range(len(rowlabels)))
    if row_ids is not None:
        plotted.insert(0, 'object_id', list(row_ids))
    view.bind_data(plotted, 'PSI' if bounded else 'display_matrix_value',
                   full_results=[('Complete matrix source', full_results, 'native_matrix_row')] if full_results else [])


def _asset(project, kind, comparison=None):
    return next((a for a in project.get("assets", []) if a.get("kind") == kind
                 and (comparison is None or a.get("metadata", {}).get("contrast_id") == comparison)
                 and Path(a.get("path", "")).is_file()), None)


def _usage_effect(view, frame, comparison):
    """Attach exact native DEXSeq fields only when its declared direction exists."""
    asset = _asset(view.project, "dexseq_result", comparison)
    if asset is None or asset.get("metadata", {}).get("native_fold_change_direction") != "control/treatment":
        return frame, None
    path = view.input(asset["path"])
    header = list(pd.read_csv(path, sep="\t", nrows=0))
    effects = [c for c in header if c.startswith("log2fold_")]
    if len(effects) != 1 or not {"groupID", "featureID"}.issubset(header):
        return frame, None
    raw = pd.read_csv(path, sep="\t", usecols=[c for c in header if c in {"groupID", "featureID", "pvalue", "padj", "exonBaseMean", *effects}])
    raw["gene_id"] = raw.groupID.astype(str).str.replace(r"\.\d+$", "", regex=True)
    raw["transcript_id"] = raw.featureID.astype(str)
    raw = raw.rename(columns={"pvalue": "native_pvalue", "padj": "native_padj", "exonBaseMean": "native_exonBaseMean"})
    merged = frame.merge(raw.drop(columns=["groupID", "featureID"]), on=["gene_id", "transcript_id"], how="left", validate="many_to_one")
    merged["native_effect_direction"] = "control/treatment (original column; not numerator/denominator)"
    return merged, effects[0]


def _cross_patterns(view, members, identity, comparisons, title, policy=None, *, group=None):
    """Exact observed UpSet patterns, with a complete searchable member view."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    if not comparisons or len(comparisons) > 4 or len(comparisons) != len(set(comparisons)):
        raise ValueError('Grouped exact intersections require one to four distinct comparison axes.')
    members = members.loc[members.comparison_id.isin(comparisons)].copy()
    scope = str((group or {}).get('group_id', 'all'))
    label = str((group or {}).get('label', 'Selected comparisons'))
    scoped_title = title + '-' + scope if scope != 'all' else title
    masks = {}
    for number, comparison in enumerate(comparisons):
        for key in set(members.loc[members.comparison_id.eq(comparison), identity].astype(str)):
            masks[key] = masks.get(key, 0) | (1 << number)
    counts = Counter(masks.values())
    ranked = sorted(counts, key=lambda mask: (-counts[mask], mask))
    pattern_ids = {mask: "P" + str(i + 1).zfill(4) for i, mask in enumerate(ranked)}
    labels = {mask: "; ".join(c for i, c in enumerate(comparisons) if mask & (1 << i)) for mask in ranked}
    complete = pd.DataFrame([{"pattern_id": pattern_ids[mask], "included_comparisons": labels[mask],
                              "comparison_count": mask.bit_count(), "object_count": counts[mask],
                              'membership_mask': mask,
                              **{comparison: bool(mask & (1 << i)) for i, comparison in enumerate(comparisons)}}
                             for mask in ranked], columns=["pattern_id", "included_comparisons", "comparison_count", "object_count", 'membership_mask', *comparisons])
    count_path = view.table(complete, scoped_title + "-all-exact-intersection-counts-" + str(policy or "common"),
                            "all", "exact_membership_pattern", policy=policy, suffix="tsv", chapter=4)
    members = members.copy()
    members["pattern_id"] = members[identity].astype(str).map(lambda key: pattern_ids[masks[key]])
    members["included_comparisons"] = members[identity].astype(str).map(lambda key: labels[masks[key]])
    member_path = view.table(members, scoped_title + "-all-intersection-members-" + str(policy or "common"),
                             "all", "object_by_comparison", policy=policy, chapter=4)
    member_tsv = view.table(members, scoped_title + '-all-members-native-values-' + str(policy or 'common'),
                           'all', 'object_by_comparison', policy=policy, suffix='tsv', chapter=4)
    set_sizes = pd.DataFrame([{'comparison_id': comparison,
        'object_count': len(set(members.loc[members.comparison_id.eq(comparison), identity].astype(str)))} for comparison in comparisons])
    set_path = view.table(set_sizes, scoped_title + '-set-sizes-' + str(policy or 'common'),
                          'all', 'unique_exact_object', policy=policy, suffix='tsv', chapter=4)
    # The complete identity browser is output-only. The native per-comparison
    # effect/P/FDR rows remain available in the full linked Parquet table.
    data = {"columns": list(members), "rows": members.astype(object).where(members.notna(), None).values.tolist()}
    template = Path(__file__).with_name("assets").joinpath("table-template.html").read_text(encoding="utf-8")
    browser = view.root / "outputs/comparisons/all" / (_safe(scoped_title + "-complete-members-" + str(policy or "common")) + ".html")
    browser.parent.mkdir(parents=True, exist_ok=True)
    page = template.replace("__TITLE__", title + ' | ' + label + " complete exact memberships")
    page = page.replace("../../REPORT.html", "../../../REPORT.html").replace("../index.html", "../../index.html")
    page = page.replace("__RAW__", Path(member_tsv).name).replace("__DATA__", json.dumps(data, allow_nan=False).replace("<", "\\u003c"))
    browser.write_text(page, encoding="utf-8")
    view.tables.append({"name": scoped_title + " complete searchable memberships | " + str(policy or "common"),
                        "path": browser.relative_to(view.root).as_posix(), "chapter": 4, "kind": "html", "rows": len(members), "policy": policy})
    visible = ranked
    fig = Figure(figsize=(14, 7.5))
    FigureCanvasAgg(fig)
    grid = fig.add_gridspec(2, 2, width_ratios=(3, 8), height_ratios=(3, 2), hspace=.05, wspace=.08)
    bar, matrix, sizes = fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1]), fig.add_subplot(grid[1, 0])
    note = fig.add_subplot(grid[0, 0]); note.axis('off')
    note.text(.05, .9, label + '\n' + str(len(masks)) + ' unique objects\n' + str(len(ranked)) + ' observed exact patterns\nAll patterns shown', va='top', fontsize=10)
    bar.bar(range(len(visible)), [counts[mask] for mask in visible], color=BLUE)
    bar.set(ylabel="Exact intersection objects", title=title + ' | ' + label + " | " + (LABELS.get(policy) or "native feature significance"))
    bar.spines[["top", "right"]].set_visible(False)
    bar.tick_params(axis="x", bottom=False, labelbottom=False)
    for x, mask in enumerate(visible):
        bar.text(x, counts[mask], str(counts[mask]), ha='center', va='bottom', fontsize=8)
    for x, mask in enumerate(visible):
        rows = [i for i in range(len(comparisons)) if mask & (1 << i)]
        matrix.scatter([x] * len(comparisons), range(len(comparisons)), color="#E2E5E9", s=17)
        matrix.scatter([x] * len(rows), rows, color=BLUE, s=24)
        if rows:
            matrix.plot([x, x], [min(rows), max(rows)], color=BLUE, lw=1)
    matrix.set_xticks(range(len(visible)), [pattern_ids[mask] for mask in visible], rotation=60, fontsize=7)
    matrix.set_yticks(range(len(comparisons)), [''] * len(comparisons))
    sizes.barh(range(len(comparisons)), set_sizes.object_count, color=RED)
    sizes.set_yticks(range(len(comparisons)), comparisons, fontsize=7)
    sizes.invert_xaxis(); sizes.invert_yaxis(); sizes.set_xlabel('Set size (unique objects)')
    largest = max(set_sizes.object_count.max(), 1)
    for i, n in enumerate(set_sizes.object_count):
        if n >= largest * .15:
            sizes.text(n - largest * .015, i, str(n), ha='left', va='center', fontsize=8, color='white')
        else:
            sizes.text(n + largest * .015, i, str(n), ha='right', va='center', fontsize=8)
    sizes.spines[['top', 'right', 'left']].set_visible(False)
    matrix.invert_yaxis(); matrix.spines[["top", "right", "bottom", "left"]].set_visible(False)
    view.figure(fig, "all-" + scoped_title + "-upset-" + str(policy or "common"), 4, title + ' | ' + label + " exact intersections",
                None, 'Memberships are recomputed using only these selected axes. Every observed nonempty pattern is shown (at most 15 for four axes), with intersection counts and set-size bars. All pattern counts, every member and unchanged native effects/P/FDR are linked; no interaction model or overlap P/FDR is fitted.', policy)
    view.panels[-1]['contrast_group_id'] = scope
    view.bind_data(complete, 'exact_membership_pattern', comparison_ids=comparisons,
        full_results=[('Set-size counts', set_path, 'unique_exact_object'),
                      ('Complete native member values', member_tsv, 'object_by_comparison'),
                      ('Search every member', browser.relative_to(view.root).as_posix(), 'object_by_comparison')])
    for link in view.panels[-1]['data_links']:
        if link['label'] == 'Set-size counts':
            link['role'] = 'plot_data'
    return {"intersection_counts": count_path, "intersection_members": member_path,
            'intersection_members_tsv': member_tsv, 'set_sizes': set_path,
            'comparison_ids': list(comparisons), 'group_id': scope,
            "complete_browser": browser.relative_to(view.root).as_posix(), "patterns": len(ranked), "objects": len(masks)}


def _direction_summary(members, comparisons):
    effects = {c: members.loc[members.comparison_id.eq(c)].drop_duplicates("event_id").set_index("event_id").dpsi
               for c in comparisons}
    records = []
    for a in comparisons:
        for b in comparisons:
            joined = pd.concat([pd.to_numeric(effects[a], errors="coerce").rename("a"),
                                pd.to_numeric(effects[b], errors="coerce").rename("b")], axis=1, join="inner")
            valid = joined.loc[np.isfinite(joined.a) & np.isfinite(joined.b) & joined.a.ne(0) & joined.b.ne(0)]
            same = int((np.sign(valid.a) == np.sign(valid.b)).sum())
            records.append({"comparison_a": a, "comparison_b": b, "shared_selected_events": len(joined),
                            "direction_comparable_events": len(valid), "same_direction_events": same,
                            "opposite_direction_events": len(valid) - same,
                            "same_direction_fraction": same / len(valid) if len(valid) else np.nan})
    return pd.DataFrame(records)


def _grouped_cross_views(view, members, identity, available, title, policy=None, *, directions=False):
    """Recalculate descriptive sets from raw selected-object rows per group."""
    result = {}
    for group in view.project['report_contrast_groups']:
        axes = group['contrast_ids']
        if len(axes) < 2:
            entry = {'status': 'not_applicable', 'contrast_ids': list(axes),
                     'reason': 'Cross-comparison analysis requires at least two selected comparisons.'}
            view.audit.append({'analysis': title, 'group_id': group['group_id'], **entry})
            result[group['group_id']] = entry
            continue
        missing = [cid for cid in axes if cid not in available]
        if missing:
            view.audit.append({'analysis': title, 'group_id': group['group_id'], 'status': 'not_computed',
                               'missing_comparisons': missing, 'reason': 'No zero sets are inferred for missing comparison results.'})
            result[group['group_id']] = {'status': 'not_computed', 'missing_comparisons': missing}
            continue
        selected = members.loc[members.comparison_id.isin(axes)].copy()
        sets = {cid: set(selected.loc[selected.comparison_id.eq(cid), identity].astype(str)) for cid in axes}
        overlap = np.asarray([[len(sets[a] & sets[b]) for b in axes] for a in axes])
        summary = pd.DataFrame([{'comparison_a': a, 'comparison_b': b,
            'shared_exact_objects': int(overlap[i, j]), 'set_size_a': len(sets[a]), 'set_size_b': len(sets[b])}
            for i, a in enumerate(axes) for j, b in enumerate(axes)])
        stem = title + '-' + group['group_id'] + '-' + str(policy or 'common')
        count_path = view.table(summary, stem + '-exact-overlap-counts', 'all',
                                 'exact_object_set_pair', policy=policy, suffix='tsv', chapter=4)
        unit = 'exact_splicing_event' if identity == 'event_id' else 'canonical_gene_and_exon_bin'
        membership_path = view.table(selected, stem + '-all-native-memberships', 'all', unit, policy=policy, chapter=4)
        if len(axes) > 1:
            fig, ax = _axes(title + ' | ' + group['label'])
            fig.set_size_inches(10, 7)
            image = ax.imshow(overlap, cmap='Blues', aspect='auto')
            ax.set_xticks(range(len(axes)), axes, rotation=45, ha='right', fontsize=8)
            ax.set_yticks(range(len(axes)), axes, fontsize=8)
            for i in range(len(axes)):
                for j in range(len(axes)):
                    ax.text(j, i, str(overlap[i, j]), ha='center', va='center', fontsize=9)
            ax.grid(False); fig.colorbar(image, ax=ax, label='Shared exact significant objects')
            view.figure(fig, stem + '-overlap', 4, title + ' overlaps | ' + group['label'], None,
                'Descriptive exact-object counts within this selected comparison group. Original IDs and significance are retained; no interaction test or biological rescue interpretation is inferred.', policy)
            view.panels[-1]['contrast_group_id'] = group['group_id']
            view.bind_data(summary, 'exact_object_set_pair', comparison_ids=axes,
                            full_results=[('Complete exact native members', membership_path, unit)])
        entry = {'status': 'available', 'memberships': membership_path, 'counts': count_path,
                 'inference': 'not_computed; descriptive exact memberships only'}
        entry.update(_cross_patterns(view, selected, identity, axes, title, policy, group=group))
        if directions:
            summary = _direction_summary(selected, axes)
            entry['direction_counts'] = view.table(summary, stem + '-exact-direction-counts', 'all',
                'exact_event_comparison_pair', policy=policy, suffix='tsv', chapter=4)
            if len(axes) > 1:
                from matplotlib.colors import LinearSegmentedColormap
                fig, ax = _axes('Exact event directions | ' + group['label'])
                fig.set_size_inches(10, 7)
                values = summary.pivot(index='comparison_a', columns='comparison_b', values='same_direction_fraction').reindex(index=axes, columns=axes)
                cmap = LinearSegmentedColormap.from_list('direction_blue_red', [BLUE, 'white', RED]); cmap.set_bad('#D9DDE2')
                image = ax.imshow(np.ma.masked_invalid(values.to_numpy()), cmap=cmap, vmin=0, vmax=1, aspect='auto')
                ax.set_xticks(range(len(axes)), axes, rotation=45, ha='right', fontsize=8)
                ax.set_yticks(range(len(axes)), axes, fontsize=8); ax.grid(False)
                fig.colorbar(image, ax=ax, label='Same-sign fraction among comparable shared events')
                view.figure(fig, stem + '-direction', 4, 'Exact event direction concordance | ' + group['label'], None,
                    'Original numerator-minus-denominator ΔPSI signs for the same selected events. Gray/NA means no shared finite nonzero effects. No new P/FDR, interaction fit or rescue classification.', policy)
                view.panels[-1]['contrast_group_id'] = group['group_id']
                view.bind_data(summary, 'exact_event_comparison_pair', comparison_ids=axes,
                                full_results=[('Complete exact native members', membership_path, unit)])
        result[group['group_id']] = entry
    return {'groups': result, 'inference': 'descriptive exact selected-group membership only; no new fits'}


def build_comparison_views(project, artifacts, destination, event_set="stable_significant"):
    """Build all requested comparison views in both DAS policies, read-only.

    Returned paths are relative to ``destination``. Large Parquet inputs are
    filtered by comparison before loading; no complete 7.8M-row usage table is
    materialized. Full output tables retain original scalar values and IDs.
    """
    from katana_next.science.das_selection import normalize_event_set
    event_set = normalize_event_set(event_set)
    project = _dict(project)
    project = _dict(project.get("project", project))
    selection = report_contrast_selection(project, project)
    project = copy.deepcopy(project)
    indexed = {row['contrast_id']: row for row in project.get('contrasts', [])}
    project['contrasts'] = [indexed[cid] for cid in selection['contrast_ids']]
    project['report_contrast_groups'] = selection['groups']
    project['report_native_contrast_ids'] = project.get('report_native_contrast_ids', selection['native_contrast_ids'])
    view = _Views(project, [_dict(a) for a in artifacts], destination)
    sources = {key: view.source(analysis, name) for key, analysis, name in (
        ("DAS", "differential.as", "Differential AS events"),
        ("DEG", "differential.deg", "Differentially expressed genes"),
        ("USAGE", "differential.dtu", "DTU and exon-usage evidence"))}
    exact_memberships = {policy: [] for policy in POLICIES}
    dexseq_memberships = []
    dexseq_available = set()
    counts_asset = _asset(project, "normalized_counts")
    counts = None
    if counts_asset:
        counts = read_table(view.input(counts_asset["path"]))
        if "gene_id" in counts:
            counts["gene_id"] = counts.gene_id.astype(str).str.replace(r"\.\d+$", "", regex=True)
            counts = counts.set_index("gene_id")
            if counts.index.has_duplicates:
                counts = counts.groupby(level=0).sum(numeric_only=True)
        else:
            counts = None
    for contrast in project.get("contrasts", []):
        comparison = contrast["contrast_id"]
        as_frame, deg, usage = (view.frame(sources[key], comparison) for key in ("DAS", "DEG", "USAGE"))
        sets, backgrounds = {}, {}
        if not deg.empty:
            deg["report_significant"] = _native_significance(deg, "deg_fdr", "log2fc", effect_cutoff=float((sources["DEG"].get("parameters") or {}).get("min_abs_log2fc", 1.)))
            full = view.table(deg, "DEG-full", comparison, "gene")
            ranked = view.top(deg, "DEG", comparison, "deg_fdr", "log2fc", "gene", full_table=full)
            _volcano(view, deg, comparison, "DEG", "log2fc", "deg_fdr", deg.report_significant,
                     xlabel="log2 fold change (numerator − denominator)", fdr_cutoff=float((sources["DEG"].get("parameters") or {}).get("fdr", .05)), full_results=full)
            sets["DEG"] = _genes(deg.loc[deg.report_significant])
            if counts is not None:
                sample_ids = contrast["numerator"]["sample_ids"] + contrast["denominator"]["sample_ids"]
                available = [s for s in sample_ids if s in counts]
                backgrounds["DEG"] = set(counts.index[counts[available].sum(axis=1).gt(0)]) & _genes(deg) if available else _genes(deg.loc[_num(deg, "deg_fdr").notna()])
                genes = _rank(deg.loc[deg.report_significant], "deg_fdr", "log2fc").head(30)
                available_genes = [g for g in genes.gene_id if g in counts.index]
                matrix = np.log2(counts.loc[available_genes, available].apply(pd.to_numeric, errors="coerce").to_numpy() + 1)
                if matrix.size:
                    matrix = (matrix - np.nanmean(matrix, axis=1, keepdims=True)) / np.maximum(np.nanstd(matrix, axis=1, keepdims=True), 1e-9)
                mapping = deg.drop_duplicates("gene_id").set_index("gene_id").gene_symbol.to_dict() if "gene_symbol" in deg else {}
                _heat(view, matrix, [str(mapping.get(g, g)) for g in available_genes], available, comparison, "DEG expression heatmap",
                      "Top 30 significant genes by native FDR/absolute effect. log2(normalized count + 1), row Z score; sample order follows the declared contrast. Full normalized counts remain in outputs.", full_results=str(counts_asset['path']), row_ids=available_genes)
            else:
                backgrounds["DEG"] = _genes(deg.loc[_num(deg, "deg_fdr").notna()])
            view.gene_set(deg, sets["DEG"], backgrounds["DEG"], "DEG", comparison)
        for kind, analysis, unit in (("dtu_result", "DTU", "transcript"), ("dexseq_result", "DEXSeq", "exon_bin")):
            features = usage.loc[usage.source_kind.eq(kind)].copy() if not usage.empty and "source_kind" in usage else pd.DataFrame()
            if features.empty:
                view.audit.append({"analysis": analysis, "comparison_id": comparison, "status": "not_computed"})
                continue
            effect = None
            if analysis == "DEXSeq":
                features, effect = _usage_effect(view, features, comparison)
            features["report_significant"] = _native_significance(features, "dtu_fdr", fdr_cutoff=float((sources["USAGE"].get("parameters") or {}).get("fdr", .05)))
            if analysis == "DEXSeq":
                dexseq_available.add(comparison)
                part = features.loc[features.report_significant,
                                    [c for c in ("comparison_id", "gene_id", "gene_symbol", "transcript_id", "dtu_fdr", effect) if c and c in features]].copy()
                part["exact_feature_id"] = part.gene_id.astype(str) + "::" + part.transcript_id.astype(str)
                dexseq_memberships.append(part)
            full = view.table(features, analysis + "-full", comparison, unit)
            ranked = view.top(features, analysis, comparison, "dtu_fdr", effect, unit, full_table=full)
            sets[analysis] = _genes(features.loc[features.report_significant])
            background_mask = _num(features, "dtu_fdr").notna()
            if analysis == "DTU":
                background_mask |= _num(features, "gene_screen_fdr").notna()
            backgrounds[analysis] = _genes(features.loc[background_mask])
            view.gene_set(features, sets[analysis], backgrounds[analysis], analysis, comparison)
            if effect:
                _volcano(view, features, comparison, "DEXSeq", effect, "dtu_fdr", features.report_significant,
                         xlabel="Native " + effect + " (control/treatment)", full_results=full)
            fig, ax = _axes(analysis + " top feature evidence | " + comparison, "−log10(native adjusted P)", unit)
            shown = ranked.iloc[::-1]
            ax.barh(range(len(shown)), -np.log10(_num(shown, "dtu_fdr").clip(lower=np.nextafter(0., 1.))), color=BLUE)
            labels = [str(r.get("gene_symbol", r.get("gene_id", ""))) + " | " + str(r.get("transcript_id", "")) for _, r in shown.iterrows()]
            ax.set_yticks(range(len(shown)), labels, fontsize=7)
            view.figure(fig, comparison + "-" + analysis + "-top-feature-evidence", 3,
                        analysis + " feature evidence | " + comparison, comparison,
                        ("Transcript-level stageR confirmation FDR; gene screening FDR remains a separate column. No inferential transcript-usage effect was provided by this completed source." if analysis == "DTU" else
                         "Exon-bin native adjusted P, not per-gene q. Exact native fold-change field and control/treatment direction are retained when declared; source transcripts and exon bins are never conflated."))
            plotted = shown.copy()
            plotted['plot_minus_log10_fdr'] = -np.log10(_num(shown, 'dtu_fdr').clip(lower=np.nextafter(0., 1.)))
            plotted['display_label'] = labels
            view.bind_data(plotted, unit, full_results=[('Complete native feature results', full, unit)])
        if not as_frame.empty:
            original = as_frame.copy(deep=False)
            ordinary = _native_significance(as_frame, "fdr", "dpsi", effect_cutoff=.1)
            eligible, rs_rule = _rs_eligibility(as_frame, contrast)
            as_frame["report_das"] = ordinary
            as_frame["report_rs_eligible"] = eligible
            as_frame["report_rs_das"] = ordinary & eligible
            as_frame["report_rs_policy"] = rs_rule
            full = view.table(as_frame, "DAS-complete-original-and-report-flags", comparison, "exact_splicing_event")
            # Complete native columns survive even when historical stable flags
            # disagree with the independently reconstructed every-sample gate.
            pd.testing.assert_frame_equal(as_frame[list(original)], original)
            raw_counts = as_frame.event_type.value_counts().reindex(EVENTS, fill_value=0)
            rs_counts = as_frame.loc[eligible].event_type.value_counts().reindex(EVENTS, fill_value=0)
            fig, ax = _axes("AS landscape | " + comparison, "Event type", "Events")
            x = np.arange(len(EVENTS))
            ax.bar(x - .2, raw_counts, .4, label="Raw completed events", color=GRAY)
            ax.bar(x + .2, rs_counts, .4, label="Read support only (every sample ≥10)", color=BLUE)
            ax.set_xticks(x, EVENTS); ax.legend(frameon=False, fontsize=8)
            view.figure(fig, comparison + "-AS-raw-versus-RS", 2, "Raw and read-supported AS | " + comparison,
                        comparison, "Read support only does not apply a significance cutoff. Raw and supported counts are from this comparison's complete exact-event table.")
            view.bind_data(pd.DataFrame({'event_type': EVENTS, 'raw_events': raw_counts.to_numpy(),
                'read_supported_events': rs_counts.to_numpy()}), 'exact_splicing_event_count',
                full_results=[('Complete events and read-support evidence', full, 'exact_splicing_event')])
            for label, values in (("Raw", raw_counts), ("Read-supported", rs_counts)):
                fig, ax = _axes(label + " event types | " + comparison)
                ax.grid(False)
                if values.sum():
                    ax.pie(values, labels=EVENTS, colors=EVENT_COLORS, autopct="%1.1f%%", startangle=90)
                else:
                    ax.text(.5, .5, "No eligible events", transform=ax.transAxes, ha="center")
                view.figure(fig, comparison + "-AS-composition-" + label, 2, label + " event composition | " + comparison,
                            comparison, "Composition of all raw or read-supported events, independently of significance.")
                view.bind_data(pd.DataFrame({'event_type': EVENTS, 'events': values.to_numpy(),
                    'percent': values.to_numpy() / values.sum() * 100 if values.sum() else np.zeros(len(EVENTS)),
                    'display_color': EVENT_COLORS}), 'exact_splicing_event_count',
                    full_results=[('Complete event population', full, 'exact_splicing_event')])
            for policy in POLICIES:
                if policy == "stable_significant" and rs_rule.startswith("missing_"):
                    view.policy_data[policy][comparison] = {"available": False, "events": None, "genes": None,
                        "raw_events": len(as_frame), "read_supported_events": None, "reason": rs_rule}
                    view.audit.append({"comparison_id": comparison, "analysis": "RS-DAS", "status": "not_computed", "reason": rs_rule})
                    continue
                selected = as_frame.report_rs_das if policy == "stable_significant" else as_frame.report_das
                display = as_frame.copy(deep=False).assign(report_significant=selected)
                member_columns = [c for c in ("comparison_id", "event_id", "gene_id", "gene_symbol", "event_type", "dpsi", "fdr", "PValue") if c in as_frame]
                exact_memberships[policy].append(as_frame.loc[selected, member_columns].copy())
                view.top(display, "DAS", comparison, "fdr", "dpsi", "exact_splicing_event", policy=policy, full_table=full)
                _volcano(view, as_frame, comparison, "DAS", "dpsi", "fdr", selected, policy=policy, xlabel="ΔPSI (numerator − denominator)", full_results=full)
                fg, bg = _genes(as_frame.loc[selected]), _genes(as_frame.loc[eligible if policy == "stable_significant" else _num(as_frame, "fdr").notna()])
                view.gene_set(as_frame, fg, bg, "DAS", comparison, policy)
                current_sets = {**sets, "DAS": fg}
                population = set().union(*current_sets.values()) if current_sets else set()
                memberships = []
                for gene in sorted(population):
                    flags = {key: gene in current_sets.get(key, set()) for key in ("DEG", "DAS", "DTU", "DEXSeq")}
                    memberships.append({"gene_id": gene, **flags, "evidence_layers": sum(flags.values()),
                                        **{key + "_available": key in current_sets for key in ("DEG", "DAS", "DTU", "DEXSeq")},
                                        "pattern": " + ".join(key for key, present in flags.items() if present)})
                matrix = pd.DataFrame(memberships, columns=["gene_id", "DEG", "DAS", "DTU", "DEXSeq", "evidence_layers", "pattern",
                                                          "DEG_available", "DAS_available", "DTU_available", "DEXSeq_available"])
                membership_file = view.table(matrix, "four-layer-memberships-" + policy, comparison, "unique_gene", policy=policy)
                fig, ax = _axes("Four-layer gene coverage | " + comparison, "Analysis layer", "Significant unique genes")
                layer_counts = [len(current_sets.get(key, set())) for key in ("DEG", "DAS", "DTU", "DEXSeq")]
                ax.bar([(LABELS[policy] if key == "DAS" else key) + ("\nnot computed" if key not in current_sets else "")
                        for key in ("DEG", "DAS", "DTU", "DEXSeq")], layer_counts, color=[RED, BLUE, "#218D87", "#9A64B0"])
                view.figure(fig, comparison + "-four-layer-coverage-" + policy, 3, "Four-layer coverage | " + comparison,
                            comparison, "Unique-gene sets from completed native significance; a gene may contribute multiple events/features but is counted once per layer. Missing layers are recorded separately, not interpreted as biological zero.", policy)
                view.bind_data(pd.DataFrame({'analysis': ['DEG', 'DAS', 'DTU', 'DEXSeq'],
                    'significant_unique_genes': layer_counts,
                    'available': [key in current_sets for key in ('DEG', 'DAS', 'DTU', 'DEXSeq')]}),
                    'unique_gene_count', full_results=[('Complete gene memberships', membership_file, 'unique_gene')])
                patterns = matrix.pattern.value_counts() if not matrix.empty else pd.Series(dtype=int)
                fig, ax = _axes("Exact four-layer gene overlaps | " + comparison, "Unique genes", "Exact evidence membership")
                ax.barh(range(len(patterns)), patterns.iloc[::-1], color=BLUE)
                ax.set_yticks(range(len(patterns)), patterns.index[::-1], fontsize=7)
                view.figure(fig, comparison + "-four-layer-overlaps-" + policy, 3, "Exact layer overlap | " + comparison,
                            comparison, "All nonempty exact membership patterns; no automatic biological interpretation or rescue inference. Complete gene membership table is linked in outputs.", policy)
                view.bind_data(pd.DataFrame({'pattern': patterns.iloc[::-1].index,
                    'unique_genes': patterns.iloc[::-1].to_numpy()}), 'unique_gene_count',
                    full_results=[('Every gene and its evidence layers', membership_file, 'unique_gene')])
                view.policy_data[policy][comparison] = {"available": True, "events": int(selected.sum()), "genes": len(fg),
                    "raw_events": len(as_frame), "read_supported_events": int(eligible.sum()),
                    "counts_by_type": {str(k): int(v) for k, v in as_frame.loc[selected].event_type.value_counts().items()},
                    "gene_sets": {key: len(value) for key, value in current_sets.items()},
                    "memberships": membership_file, "missing_layers": [key for key in ("DEG", "DTU", "DEXSeq") if key not in sets]}
            view.audit.append({"comparison_id": comparison, "analysis": "DAS", "status": "display_derived",
                               "caller_values_changed": False, "rs_rule": rs_rule,
                               "historical_rs_flags_reused": False, "rs_eligible": int(eligible.sum()),
                               "native_stable_significant": int(as_frame.stable_significant.sum()) if "stable_significant" in as_frame else None})
    # New display-only exact memberships are recomputed from the declared
    # report foreground. Historical overlap P/FDR is never relabelled as this
    # collection's inference; only exact counts and direction signs are shown.
    cross = {}
    for policy, parts in exact_memberships.items():
        if not parts:
            continue
        members = pd.concat(parts, ignore_index=True)
        available = [c["contrast_id"] for c in project.get("contrasts", [])
                     if view.policy_data[policy].get(c["contrast_id"], {}).get("available")]
        cross[policy] = _grouped_cross_views(view, members, "event_id", available,
                                           "DAS-events", policy, directions=True)
    if dexseq_memberships:
        members = pd.concat(dexseq_memberships, ignore_index=True)
        cross["DEXSeq"] = _grouped_cross_views(view, members, "exact_feature_id",
                                               dexseq_available, "DEXSeq-features")
    for path, initial in view.snapshots.items():
        now = path.stat()
        if initial != (now.st_size, now.st_mtime_ns):
            raise RuntimeError("Completed source changed during report construction: " + str(path))
    return {"panels": view.panels, "tables": view.tables, "top10": view.top10,
            "policy_data": view.policy_data, "gene_sets": view.gene_sets, "cross_comparison": cross, "audit": view.audit,
            "default_event_set": event_set, "new_scientific_fits": 0,
            "source_invariance": {"files_checked": len(view.snapshots), "unchanged": True}}
