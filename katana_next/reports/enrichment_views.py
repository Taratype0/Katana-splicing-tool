"""Draw completed, explicitly scoped ORA tables without calculating enrichment.

Executive consulting editorial visual overhaul:
- Strips internal database prefixes (`mouse::GO_MGI_mouse / `) on display ticks only;
  100% of raw identifiers, libraries, foreground/background sets, and BH families remain
  intact in `plotted.tsv`, `complete-tests.parquet`, and `data_links`.
- Dynamic per-row vertical allocation + max 2-line wrapping eliminates all text overlap.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import re
import shutil
import textwrap

import numpy as np
import pandas as pd

from .comparisons import _safe
from .sources import read_table, table_schema

CONTRACT = "katana.report_scoped_ora/v1"

MCK_NAVY = "#051C2C"
MCK_BLUE = "#2251FF"
MCK_CERULEAN = "#0077C8"
MCK_PALE = "#93C5FD"
MCK_CORAL = "#B83A4B"
MCK_SLATE = "#64748B"

_GO_ID_RE = re.compile(r"(GO:\d{7})")


@lru_cache(maxsize=1)
def _load_go_category_map() -> dict[str, str]:
    res = Path(__file__).resolve().parents[1] / "resources" / "go_categories.json.gz"
    if res.is_file():
        raw = json.loads(gzip.decompress(res.read_bytes()).decode("utf-8"))
        return raw.get("categories", raw)
    return {}


def _infer_ontology_category(row_cat: str, term_text: str, lib_text: str = "") -> str:
    norm_map = {
        "BP": "GO Biological Process", "GO BP": "GO Biological Process", "GO:BP": "GO Biological Process", "biological_process": "GO Biological Process",
        "CC": "GO Cellular Component", "GO CC": "GO Cellular Component", "GO:CC": "GO Cellular Component", "cellular_component": "GO Cellular Component",
        "MF": "GO Molecular Function", "GO MF": "GO Molecular Function", "GO:MF": "GO Molecular Function", "molecular_function": "GO Molecular Function",
        "KEGG": "KEGG Pathways", "KEGG Pathways": "KEGG Pathways",
    }
    rc = str(row_cat or "").strip()
    if rc in norm_map:
        return norm_map[rc]
    comb = f"{lib_text} {term_text}".lower()
    if "kegg" in comb or re.search(r"\b(?:hsa|mmu)\d{5}\b", comb):
        return "KEGG Pathways"
    m = _GO_ID_RE.search(str(term_text or ""))
    if m:
        mapped = _load_go_category_map().get(m.group(1), "")
        if mapped in norm_map:
            return norm_map[mapped]
    if rc and rc not in {"Unclassified", "nan", "None"}:
        return rc
    return "Enriched Terms"


def clean_ora_term_label(raw_term: str, wrap_width: int = 46) -> str:
    """Format a clean human-readable GO/KEGG/Reactome label with at most 2 wrapped lines."""
    text = str(raw_term or "").strip()
    if " / " in text:
        _, text = text.split(" / ", 1)
    text = re.sub(r"^[A-Za-z0-9_.-]+::[A-Za-z0-9_.-]+\s*/?\s*", "", text).strip()
    go_match = re.match(r"^(GO:\d{7})\s+(.+)$", text)
    if go_match:
        go_id, desc = go_match.group(1), go_match.group(2).strip()
        desc = desc[0].upper() + desc[1:] if desc else desc
        clean = f"{desc} ({go_id})"
    else:
        # Remove trailing species suffix like ` - Mus musculus (house mouse)`
        clean = re.sub(r"\s+-\s+(?:Mus musculus|Homo sapiens).*$", "", text, flags=re.I).strip()
        clean = clean[0].upper() + clean[1:] if clean else clean
    lines = textwrap.wrap(clean, width=wrap_width, break_long_words=False, break_on_hyphens=False)
    if not lines:
        return clean[:wrap_width]
    if len(lines) > 2:
        second = lines[1]
        if len(second) > wrap_width - 2:
            second = second[: wrap_width - 2].rstrip() + "…"
        else:
            second = second + "…"
        lines = [lines[0], second]
    return "\n".join(lines)


def render_mckinsey_ora_facets(
    plotted_df: pd.DataFrame,
    png_path: str | Path,
    svg_path: str | Path | None = None,
    title: str = "GO & KEGG Functional Enrichment",
    subtitle: str = "Horizontal bars show −log10(BH FDR); circle size encodes contributing genes; color encodes gene ratio",
    complete_df: pd.DataFrame | None = None,
) -> None:
    """Render a fixed multi-ontology (BP / CC / MF / KEGG) executive lollipop chart.

    Guarantees KEGG is never silently dropped when tested in `complete_df`:
    categories with 0 BH-FDR ≤ 0.05 hits render their top nominal candidates (`overlap ≥ 2`)
    as clearly labeled hollow exploratory lollipops to the left of the FDR = 0.05 threshold.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    def _prepare_frame(src: pd.DataFrame | None, is_exploratory_pool: bool = False) -> pd.DataFrame:
        if src is None or src.empty:
            return pd.DataFrame()
        w = src.copy()
        cat_col = None
        for cand_c in ("category", "ontology_category", "source", "facet", "group"):
            if cand_c in w.columns:
                vals_u = set(w[cand_c].dropna().astype(str)) - {"Unclassified", ""}
                if vals_u:
                    cat_col = cand_c
                    break
        if cat_col is None:
            cat_col = next((c for c in ("category", "ontology_category", "source") if c in w.columns), None)
        term_col = next((c for c in ("term", "term_name", "Description", "ID", "pathway") if c in w.columns), w.columns[0])
        lib_col = next((c for c in ("library", "library_kind", "source") if c in w.columns), None)
        fdr_col = next((c for c in ("fdr", "FDR", "p.adjust", "padj", "qvalue") if c in w.columns), None)
        p_col = next((c for c in ("pvalue", "p_value", "P") if c in w.columns), fdr_col)
        ratio_col = next((c for c in ("gene_ratio", "plot_gene_ratio", "GeneRatio") if c in w.columns), None)
        count_col = next((c for c in ("overlap_genes", "intersection_size", "Count", "hits", "count") if c in w.columns), None)

        def _parse_ratio(v):
            if isinstance(v, str) and "/" in v:
                a, b = v.split("/", 1)
                try:
                    return float(a) / max(1.0, float(b))
                except ValueError:
                    return np.nan
            return pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]

        w["_term_raw"] = w[term_col].astype(str)
        w["_fdr"] = pd.to_numeric(w[fdr_col], errors="coerce").clip(lower=1e-300) if fdr_col else 1.0
        w["_pval"] = pd.to_numeric(w[p_col], errors="coerce").clip(lower=1e-300) if p_col else w["_fdr"]
        w["_neglog10"] = -np.log10(w["_fdr"])
        if ratio_col:
            w["_ratio"] = w[ratio_col].map(_parse_ratio)
        elif count_col and "selected_genes" in w.columns:
            w["_ratio"] = pd.to_numeric(w[count_col], errors="coerce") / pd.to_numeric(w["selected_genes"], errors="coerce").replace(0, np.nan)
        elif count_col and "query_size" in w.columns:
            w["_ratio"] = pd.to_numeric(w[count_col], errors="coerce") / pd.to_numeric(w["query_size"], errors="coerce").replace(0, np.nan)
        else:
            w["_ratio"] = 0.15
        w["_count"] = pd.to_numeric(w[count_col], errors="coerce").fillna(0.0) if count_col else 2.0

        terms_list = w["_term_raw"].tolist()
        cats_raw = w[cat_col].astype(str).tolist() if cat_col else [""] * len(w)
        libs_raw = w[lib_col].astype(str).tolist() if lib_col else [""] * len(w)
        w["_cat"] = [_infer_ontology_category(rc, tm, lb) for rc, tm, lb in zip(cats_raw, terms_list, libs_raw)]
        w["_exploratory"] = is_exploratory_pool | (w["_fdr"] > 0.05)
        return w

    df = _prepare_frame(plotted_df, is_exploratory_pool=False)
    comp = _prepare_frame(complete_df, is_exploratory_pool=True)

    all_standard = ["GO Biological Process", "GO Cellular Component", "GO Molecular Function", "KEGG Pathways"]
    tested_cats = set(df["_cat"].unique()) if not df.empty else set()
    if not comp.empty:
        tested_cats |= set(comp["_cat"].unique())

    if df.empty and comp.empty:
        fig, ax = plt.subplots(figsize=(10.6, 2.8), dpi=220)
        fig.patch.set_facecolor("white")
        ax.axis("off")
        ax.add_patch(plt.Rectangle((0.01, 0.12), 0.98, 0.72, facecolor="#F8FAFC", edgecolor="#E2E8F0", lw=1.0, transform=ax.transAxes))
        ax.text(0.03, 0.58, title, transform=ax.transAxes, fontsize=11, fontweight="bold", color=MCK_NAVY)
        ax.text(
            0.03, 0.30,
            "No ontology or pathway terms passed the significance and display-overlap criteria (FDR ≤ 0.05, overlap ≥ 2). Complete tested sets remain linked below.",
            transform=ax.transAxes, fontsize=9.0, color=MCK_SLATE,
        )
        fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
        if svg_path:
            fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return

    ordered_cats = [c for c in all_standard if c in tested_cats] + [
        c for c in dict.fromkeys(list(df.get("_cat", [])) + list(comp.get("_cat", [])))
        if c not in set(all_standard)
    ]

    per_cat_limit = 5 if len(ordered_cats) > 1 else 12
    blocks: list[tuple[str, pd.DataFrame, bool]] = []
    for cat in ordered_cats:
        sub_sig = df.loc[(df["_cat"] == cat) & (df["_fdr"] <= 0.05) & (df["_count"] >= 2)].sort_values(["_fdr", "_pval", "_term_raw"], kind="stable").head(per_cat_limit).copy() if not df.empty else pd.DataFrame()
        if sub_sig.empty and not df.empty and comp.empty:
            sub_sig = df.loc[df["_cat"] == cat].sort_values(["_fdr", "_pval", "_term_raw"], kind="stable").head(per_cat_limit).copy()
        if not sub_sig.empty:
            sub_sig["_clean_label"] = sub_sig["_term_raw"].map(lambda t: clean_ora_term_label(t, wrap_width=46))
            blocks.append((cat, sub_sig, False))
        elif not comp.empty:
            pool = comp.loc[(comp["_cat"] == cat) & (comp["_count"] >= 2)].sort_values(["_pval", "_fdr", "_term_raw"], kind="stable")
            if pool.empty:
                pool = comp.loc[(comp["_cat"] == cat) & (comp["_count"] >= 1)].sort_values(["_pval", "_fdr", "_term_raw"], kind="stable")
            sub_exp = pool.loc[pool["_pval"] <= 0.05].head(4).copy()
            if sub_exp.empty:
                sub_exp = pool.head(3).copy()
            n_tot = int((comp["_cat"] == cat).sum())
            if not sub_exp.empty:
                sub_exp["_exploratory"] = True
                sub_exp["_clean_label"] = sub_exp["_term_raw"].map(lambda t: clean_ora_term_label(t, wrap_width=46))
                blocks.append((f"{cat}  [0 at BH FDR ≤ 0.05 · Top nominal candidates of {n_tot} tested]", sub_exp, True))
            else:
                placeholder = pd.DataFrame([{
                    "_term_raw": f"No pathway with ≥2 overlapping genes ({n_tot} pathways tested)",
                    "_clean_label": f"0 terms with ≥2 gene overlap ({n_tot} tested)",
                    "_fdr": 1.0, "_pval": 1.0, "_neglog10": 0.0, "_ratio": 0.0, "_count": 0.0, "_exploratory": True,
                }])
                blocks.append((f"{cat}  [0 at BH FDR ≤ 0.05]", placeholder, True))

    if not blocks:
        return

    height_ratios = [max(1.05, len(sub) * 0.56 + 0.38) for _, sub, _ in blocks]
    fig_h = max(4.2, sum(height_ratios) + 1.70)

    fig, axes = plt.subplots(
        len(blocks), 1,
        figsize=(12.0, fig_h),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.56},
        sharex=False,
        dpi=220,
    )
    if len(blocks) == 1:
        axes = [axes]
    fig.patch.set_facecolor("white")

    cmap = LinearSegmentedColormap.from_list("mck_seq", ["#60A5FA", "#2251FF", "#051C2C"])
    r_vals = pd.concat([b[1]["_ratio"] for b in blocks]).dropna()
    vmin = float(r_vals.min()) if len(r_vals) else 0.0
    vmax = float(r_vals.max()) if len(r_vals) else 0.5
    if abs(vmax - vmin) < 1e-6:
        vmax = vmin + 0.1
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    max_cnt = max(1.0, float(max(b[1]["_count"].max() for b in blocks)))

    last_scatter = None
    for ax, (cat_title, sub, is_exp_panel) in zip(axes, blocks):
        sub_rev = sub.iloc[::-1].reset_index(drop=True)
        ys = np.arange(len(sub_rev))
        pvals = sub_rev["_pval"].to_numpy(dtype=float)
        xs = (-np.log10(np.clip(pvals, 1e-300, 1.0))) if is_exp_panel else sub_rev["_neglog10"].to_numpy(dtype=float)
        cnts = sub_rev["_count"].to_numpy(dtype=float)
        ratios = sub_rev["_ratio"].fillna(vmin).to_numpy(dtype=float)
        panel_max = float(np.max(xs)) if len(xs) else 2.0
        x_upper = max(2.35, panel_max * 1.28)

        ax.grid(axis="x", color="#E2E8F0", linewidth=0.6, zorder=0)
        ax.axvline(-np.log10(0.05), color="#94A3B8", linestyle="--" if not is_exp_panel else ":", linewidth=0.95, zorder=1)
        ax.hlines(ys, 0, xs, color="#CBD5E1" if not is_exp_panel else "#BAE6FD", linewidth=2.1, linestyle="-" if not is_exp_panel else "--", zorder=2)
        sizes = 48.0 + 220.0 * (cnts / max_cnt)
        if is_exp_panel:
            ax.scatter(
                xs, ys, s=np.where(cnts > 0, sizes * 0.85, 0.0),
                facecolors="white", edgecolors=MCK_CERULEAN, linewidths=1.6, zorder=4
            )
            for y_i, x_i, c_i, p_i in zip(ys, xs, cnts, pvals):
                if c_i > 0:
                    ax.text(
                        x_i + x_upper * 0.016, y_i,
                        f"n={int(round(c_i))} · nominal P={p_i:.3g} (ns)",
                        va="center", ha="left", fontsize=7.6, color=MCK_SLATE, style="italic",
                    )
        else:
            last_scatter = ax.scatter(
                xs, ys, s=sizes, c=ratios, cmap=cmap, norm=norm,
                edgecolors="white", linewidths=0.8, zorder=4
            )
            for y_i, x_i, c_i, r_i in zip(ys, xs, cnts, ratios):
                ax.text(
                    x_i + x_upper * 0.016, y_i,
                    f"n={int(round(c_i))} ({r_i * 100:.1f}%)",
                    va="center", ha="left", fontsize=7.8, color=MCK_SLATE,
                )

        ax.set_yticks(ys)
        ax.set_yticklabels(sub_rev["_clean_label"].tolist(), fontsize=8.5, color="#1E293B" if not is_exp_panel else "#475569", linespacing=1.12)
        ax.set_ylim(-0.62, len(sub_rev) - 0.38)
        ax.set_xlim(0, x_upper)
        ax.tick_params(axis="x", labelsize=7.8)
        ax.spines["left"].set_color("#94A3B8")
        ax.spines["bottom"].set_color("#94A3B8")
        ax.set_title(cat_title, loc="left", fontsize=9.5, fontweight="bold", color=MCK_NAVY if not is_exp_panel else "#334155", pad=5)

    axes[-1].set_xlabel("Significance: −log10(BH FDR) for solid markers; −log10(nominal P) for hollow exploratory panels   [line = 0.05]", fontsize=8.8, color=MCK_NAVY, labelpad=7)

    fig.subplots_adjust(left=0.40, right=0.86, top=1.0 - (1.10 / fig_h), bottom=0.64 / fig_h)

    if last_scatter is not None:
        cax = fig.add_axes([0.89, 0.34, 0.018, min(0.42, 1.8 / fig_h + 0.12)])
        cb = fig.colorbar(last_scatter, cax=cax)
        cb.set_label("Gene ratio", fontsize=8.5, color=MCK_NAVY)
        cb.ax.tick_params(labelsize=7.8)
        cb.outline.set_edgecolor("#CBD5E1")

    fig.text(0.02, 1.0 - (0.18 / fig_h), title, fontsize=12.0, fontweight="bold", color=MCK_NAVY, ha="left", va="top")
    fig.text(0.02, 1.0 - (0.48 / fig_h), subtitle, fontsize=8.4, color=MCK_SLATE, ha="left", va="top")

    fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
    if svg_path:
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _scope_inputs(metadata, scope, frame):
    """Resolve the producer's exact set declarations, never infer from names."""
    declared = metadata.get("report_enrichment_scope_audit")
    if not declared:
        return []
    root = Path(metadata["report_enrichment_root"]).resolve()
    audit = pd.read_csv(declared, sep="\t", dtype=str)
    libraries = frame.library.dropna().astype(str).unique() if "library" in frame else []
    if len(libraries) != 1:
        raise ValueError("A scoped ORA panel requires one explicit library identity")
    chosen = audit
    for key in ("comparison_id", "evidence", "report_policy"):
        chosen = chosen.loc[chosen[key].eq(str(scope[key]))]
    chosen = chosen.loc[chosen.library.eq(libraries[0])]
    if len(chosen) != 1 or chosen.iloc[0].fdr_family != frame.fdr_family.iloc[0]:
        raise ValueError("ORA foreground/background declaration does not match its comparison and FDR family")
    row = chosen.iloc[0]
    links = []
    for key, role, label in (("foreground_path", "foreground", "Input genes"),
                             ("background_path", "background", "Detected / testable genes"),
                             ("actual_foreground_path", "foreground", "Input genes mapped to this library"),
                             ("actual_background_path", "background", "Actual enrichment background")):
        path = (root / row[key]).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError("Missing or unbounded declared enrichment set: " + str(path))
        links.append({"role": role, "label": label, "path": str(path),
                      "contrast_id": str(scope["comparison_id"]),
                      "policy": None if scope["report_policy"] == "common" else scope["report_policy"],
                      "unit": "gene", "fdr_family": row.fdr_family,
                      "library": libraries[0]})
    for key, label in (("report_enrichment_mapping", "Gene mapping and library identifiers"),
                       ("report_enrichment_methods", "Methods and resource versions")):
        if metadata.get(key):
            path = Path(metadata[key]).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            links.append({"role": "methods", "label": label, "path": str(path),
                          "contrast_id": "all", "policy": None, "scope": "shared", "unit": "resource"})
    return links


def enrichment_attachment_tables(metadata):
    rows = []
    for raw in metadata.get("report_enrichment_attachments", []):
        record = dict(raw)
        path = Path(record.get("path", "")).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError("Declared enrichment attachment is unavailable: " + str(path))
        identity = hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:16]
        record.update({"id": "ora-attachment-" + identity,
                       "path": str(path), "source": str(path),
                       "name": record.get("name", path.name), "chapter": record.get("chapter", 5)})
        record["kind"] = record.get("kind", "file" if path.suffix.lower() in {".json", ".html"} else "table")
        if record["kind"] != "file":
            record.update(table_schema(path))
        rows.append(record)
    return rows


def build_enrichment_views(project, artifacts, destination):
    """Read only certified report-scoped completed tables; do not run ORA/BH."""
    output = Path(destination).resolve()
    panels, tables, audit = [], [], []
    for raw in artifacts:
        artifact = raw.model_dump(mode="json") if hasattr(raw, "model_dump") else dict(raw)
        meta = artifact.get("metadata", {})
        if meta.get("report_enrichment_contract") != CONTRACT:
            continue
        tables.extend(enrichment_attachment_tables(meta))
        source = Path(artifact.get("path", ""))
        if not source.is_file() or artifact.get("status", "complete") not in ("complete", "completed"):
            audit.append({"name": artifact.get("name"), "status": "not_computed", "reason": "Completed enrichment table unavailable."})
            continue
        required = {"comparison_id", "evidence", "report_policy", "library_kind", "ontology_category",
                    "term", "pvalue", "fdr", "fdr_family", "overlap_genes", "selected_genes", "background_genes"}
        schema = table_schema(source)
        if not required.issubset(schema["columns"]):
            raise ValueError("Scoped ORA table lacks explicit comparison, policy, resource, or correction family.")
        stamp = (source.stat().st_size, source.stat().st_mtime_ns)
        with source.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        directory = output / "outputs/enrichment"
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / (digest[:12] + "-" + source.name)
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        tables.append({"id": "scoped-ora-" + digest[:12], "name": "Complete measured-background ORA tests",
                       "path": target.relative_to(output).as_posix(), "source": target.relative_to(output).as_posix(),
                       "rows": schema["rows"], "columns": schema["columns"], "unit": "library_term_test", "chapter": 5})
        keys = ["comparison_id", "evidence", "report_policy", "library_kind", "ontology_category"]
        scopes = read_table(source, columns=keys).drop_duplicates()
        for scope in scopes.to_dict("records"):
            if str(scope["comparison_id"]) not in {c["contrast_id"] for c in project.get("contrasts", [])} and project.get("contrasts"):
                continue
            filters = [(key, "=", value) for key, value in scope.items()]
            frame = read_table(source, filters=filters) if source.suffix == ".parquet" else read_table(source)
            if source.suffix != ".parquet":
                for key, value in scope.items():
                    frame = frame.loc[frame[key].eq(value)]
            families = list(frame.fdr_family.dropna().astype(str).unique())
            if len(families) != 1:
                raise ValueError("A displayed ORA scope must declare exactly one correction family.")
            visible = frame.loc[pd.to_numeric(frame.overlap_genes, errors="coerce").ge(2)].copy()
            visible = visible.sort_values(["fdr", "pvalue", "term"], kind="stable").head(12)
            category = str(scope["ontology_category"])
            evidence = str(scope["evidence"])
            comparison = str(scope["comparison_id"])
            policy = scope["report_policy"] if scope["report_policy"] != "common" else None
            chapter = 7 if scope["library_kind"] == "RBP_identity" else 5
            title = f"{category} | {evidence} | {comparison}"

            identity = _safe("ORA-" + "-".join(str(scope[key]) for key in keys))
            data_directory = output / "outputs/enrichment/plot-data"
            data_directory.mkdir(parents=True, exist_ok=True)
            plotted = visible.iloc[::-1].copy()
            plotted["plot_order"] = np.arange(len(plotted))
            plotted["plot_gene_ratio"] = pd.to_numeric(plotted.overlap_genes).div(pd.to_numeric(plotted.selected_genes).replace(0, np.nan))
            plotted["plot_neglog10_fdr"] = -np.log10(pd.to_numeric(plotted.fdr).clip(lower=np.nextafter(0., 1.)))
            plot_path = data_directory / (identity + "-plotted.tsv")
            plotted.to_csv(plot_path, sep="\t", index=False)
            full_path = data_directory / (identity + "-complete-tests.parquet")
            frame.to_parquet(full_path, index=False)

            links = [{"role": "plot_data", "label": "Plotted terms and contributing genes", "path": str(plot_path),
                      "contrast_id": comparison, "policy": policy, "unit": "library_term_test", "fdr_family": families[0]},
                     {"role": "full_results", "label": "All tested terms in this panel", "path": str(full_path),
                      "contrast_id": comparison, "policy": policy, "unit": "library_term_test", "fdr_family": families[0]},
                     *_scope_inputs(meta, scope, frame)]
            image = output / "figures/enrichment" / (identity + ".png")
            render_mckinsey_ora_facets(plotted, image, image.with_suffix(".svg"), title=title)

            rbp_note = (" RBP identity membership only; this is not an RBP-target, motif-binding, or causal-regulation test."
                        if chapter == 7 else " GO categories are display subsets of the original GO-library correction family.")
            panels.append({"id": identity, "chapter": chapter, "title": title,
                           "path": image.relative_to(output).as_posix(), "contrast_id": comparison,
                           "scope": "comparison", "policy": policy,
                           "data_links": links,
                           "caption": "Completed downstream hypergeometric ORA; foreground and native detected/testable background files are in outputs. "
                                      + "BH family: " + families[0] + ". Top 12 displayed by FDR then P; overlap ≥2 is display only. "
                                      + "All eligible zero-overlap and nonsignificant tests remain in the complete table." + rbp_note})
            audit.append({**scope, "status": "rendered_completed_tests", "eligible_tests": len(frame),
                          "displayed_terms": len(visible), "fdr_family": families[0], "new_inference": False})
        if stamp != (source.stat().st_size, source.stat().st_mtime_ns):
            raise RuntimeError("Completed enrichment input changed during rendering.")
    return {"panels": panels, "tables": tables, "audit": audit, "new_inference": False, "new_scientific_fits": 0}
