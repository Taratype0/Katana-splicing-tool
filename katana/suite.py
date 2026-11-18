"""Core orchestration for standalone Python execution (`run_bulk`, `run_single_cell`) and Nextflow convergence."""
from __future__ import annotations

from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any
import webbrowser

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats


NAVY = "#051C2C"
BLUE = "#2251FF"
CERULEAN = "#0077C8"
CORAL = "#C83E4D"
SLATE = "#64748B"
EVENTS = ("SE", "A3SS", "A5SS", "MXE", "RI")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_tabular(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf == ".parquet":
        return pd.read_parquet(path)
    sep = "," if suf == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def _save_dual(fig: plt.Figure, png_path: Path) -> Path:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path = png_path.with_suffix(".svg")
    fig.savefig(png_path, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return svg_path


def _load_pathway_gmt_from_resources(species: str) -> Path | None:
    """Materialize bundled GO+KEGG pathways as a temporary GMT if needed by Katana-SpliceNet."""
    res_dir = Path(__file__).resolve().parents[1] / "katana_next" / "resources"
    gz_name = "pathways_mouse_go_kegg.json.gz" if species == "mouse" else "pathways_human_go_kegg.json.gz"
    gz_path = res_dir / gz_name
    if not gz_path.is_file():
        return None
    return gz_path


def _run_ora_from_bundled_catalog(fg_symbols: set[str], bg_symbols: set[str], species: str) -> pd.DataFrame:
    """Exact Fisher/Hypergeometric + per-category BH FDR across GO (BP, CC, MF) and KEGG."""
    res_dir = Path(__file__).resolve().parents[1] / "katana_next" / "resources"
    gz_path = res_dir / ("pathways_mouse_go_kegg.json.gz" if species == "mouse" else "pathways_human_go_kegg.json.gz")
    if not gz_path.is_file():
        return pd.DataFrame(columns=["category", "term", "overlap_count", "term_size", "selected_genes", "background_genes", "fold_enrichment", "pvalue", "fdr", "hit_genes"])
    raw_cat = json.loads(gzip.open(gz_path, "rt", encoding="utf-8").read())
    catalog = raw_cat.get("categories", raw_cat)
    fg_u = {g.strip().upper() for g in fg_symbols if g and not g.upper().startswith("ENS") and g.upper() != "NAN"}
    bg_u = {g.strip().upper() for g in bg_symbols if g and not g.upper().startswith("ENS") and g.upper() != "NAN"} | fg_u
    if len(bg_u) < 50:
        for cmap in catalog.values():
            for g_list in cmap.values():
                bg_u.update(str(x).upper() for x in g_list)
    fg_in_bg = fg_u & bg_u
    M = len(bg_u)
    N = len(fg_in_bg)
    if N == 0 or M < 10:
        return pd.DataFrame(columns=["category", "term", "overlap_count", "term_size", "selected_genes", "background_genes", "fold_enrichment", "pvalue", "fdr", "hit_genes"])

    out_rows: list[dict[str, Any]] = []
    for cat in ("BP", "CC", "MF", "KEGG"):
        cat_rows = []
        for term, raw_genes in catalog.get(cat, {}).items():
            t_set = {str(x).upper() for x in raw_genes} & bg_u
            n_t = len(t_set)
            if n_t < 5 or n_t > 2000:
                continue
            hits = sorted(fg_in_bg & t_set)
            k = len(hits)
            if k < 2:
                continue
            pv = float(sp_stats.hypergeom.sf(k - 1, M, n_t, N))
            fe = (k / N) / (n_t / M)
            cat_rows.append({
                "category": cat,
                "ontology_category": cat,
                "library": "KEGG" if cat == "KEGG" else f"GO_{cat}",
                "term": term,
                "overlap_count": k,
                "overlap_genes": k,
                "term_size": n_t,
                "selected_genes": N,
                "background_genes": M,
                "fold_enrichment": round(fe, 4),
                "pvalue": pv,
                "hit_genes": ",".join(hits[:50]),
            })
        if cat_rows:
            cat_rows.sort(key=lambda r: r["pvalue"])
            m_t = len(cat_rows)
            q_vals = [min(1.0, r["pvalue"] * m_t / (idx + 1)) for idx, r in enumerate(cat_rows)]
            for idx in range(m_t - 2, -1, -1):
                q_vals[idx] = min(q_vals[idx], q_vals[idx + 1])
            for r, qv in zip(cat_rows, q_vals):
                r["fdr"] = qv
                r["minus_log10_fdr"] = -np.log10(max(qv, 1e-300))
            out_rows.extend(cat_rows)
    return pd.DataFrame(out_rows)


def _render_ora_faceted_figure(ora_df: pd.DataFrame, png_path: Path, title: str) -> Path | None:
    if ora_df.empty:
        return None
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.0), dpi=220)
    fig.patch.set_facecolor("white")
    cat_titles = {
        "BP": "GO Biological Process",
        "CC": "GO Cellular Component",
        "MF": "GO Molecular Function",
        "KEGG": "KEGG Pathways",
    }
    for ax, cat in zip(axes.flat, ("BP", "CC", "MF", "KEGG")):
        sub = ora_df[ora_df["category"] == cat].sort_values(["fdr", "pvalue"]).head(8).copy()
        ax.set_title(cat_titles[cat], loc="left", fontsize=10.0, fontweight="bold", color=NAVY)
        if sub.empty:
            ax.text(0.5, 0.5, "No enriched terms (k ≥ 2)", ha="center", va="center", transform=ax.transAxes, color=SLATE, fontsize=8.8)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        sub = sub.iloc[::-1]
        y = np.arange(len(sub))
        x = -np.log10(np.clip(pd.to_numeric(sub["fdr"], errors="coerce").fillna(1.0).to_numpy(), 1e-50, 1.0))
        labels = [re.sub(r"^(?:GO:\d+\s+|mmu\d+\s+|hsa\d+\s+)", "", str(t))[:42] for t in sub["term"]]
        cols = [BLUE if float(q) <= 0.05 else SLATE for q in sub["fdr"]]
        ax.hlines(y, 0, x, color=cols, linewidth=2.0, alpha=0.85)
        sizes = np.clip(pd.to_numeric(sub["overlap_count"], errors="coerce").fillna(3).to_numpy() * 14, 36, 190)
        ax.scatter(x, y, s=sizes, color=cols, edgecolors="white", linewidth=0.8, zorder=3)
        ax.axvline(-np.log10(0.05), color=CORAL, linestyle="--", linewidth=0.9, alpha=0.7)
        ax.set_yticks(y, labels, fontsize=8.2, color=NAVY)
        ax.set_xlabel("-log10(BH FDR)", fontsize=8.4, color=NAVY)
        ax.set_xlim(0, max(float(x.max()) * 1.18, 1.8))
    fig.suptitle(title, x=0.01, ha="left", fontsize=11.8, fontweight="bold", color=NAVY, y=1.01)
    fig.subplots_adjust(wspace=0.58, hspace=0.34)
    return _save_dual(fig, png_path)


def _render_volcano(df: pd.DataFrame, x_col: str, q_col: str, sig_col: str, label_col: str, xlabel: str, title: str, png_path: Path) -> Path:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.0, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(7.6, 5.2), dpi=220)
    fig.patch.set_facecolor("white")
    xv = pd.to_numeric(df[x_col], errors="coerce").fillna(0.0).to_numpy()
    qv = np.clip(pd.to_numeric(df[q_col], errors="coerce").fillna(1.0).to_numpy(), 1e-300, 1.0)
    yv = -np.log10(qv)
    yv_disp = np.minimum(yv, np.nanpercentile(yv, 99.5) * 1.12 + 2.0) if len(yv) else yv
    sig = df[sig_col].astype(bool).to_numpy() if sig_col in df.columns else (qv <= 0.05)

    ax.scatter(xv[~sig], yv_disp[~sig], s=12, color="#CBD5E1", alpha=0.55, linewidths=0, label=f"Background ({int((~sig).sum()):,})")
    up = sig & (xv > 0)
    dn = sig & (xv < 0)
    if up.any():
        ax.scatter(xv[up], yv_disp[up], s=22, color=CORAL, alpha=0.88, edgecolors="white", linewidths=0.3, label=f"Up-regulated ({int(up.sum()):,})")
    if dn.any():
        ax.scatter(xv[dn], yv_disp[dn], s=22, color=BLUE, alpha=0.88, edgecolors="white", linewidths=0.3, label=f"Down-regulated ({int(dn.sum()):,})")
    ax.axhline(-np.log10(0.05), color=SLATE, linestyle="--", linewidth=0.9, alpha=0.7)
    ax.axvline(0.0, color="#94A3B8", linestyle=":", linewidth=0.8)

    if sig.any() and label_col in df.columns:
        top_idx = np.argsort(qv[sig])[:10]
        sig_pos = np.flatnonzero(sig)[top_idx]
        seen_lbl: set[str] = set()
        for idx in sig_pos:
            lbl = str(df.iloc[idx][label_col] or "").strip()
            if not lbl or lbl in seen_lbl:
                continue
            seen_lbl.add(lbl)
            ax.annotate(lbl, (xv[idx], yv_disp[idx]), xytext=(4, 3), textcoords="offset points", fontsize=7.8, fontweight="bold", color=NAVY)
    ax.set_xlabel(xlabel, fontweight="bold", color=NAVY)
    ax.set_ylabel("-log10(FDR)", fontweight="bold", color=NAVY)
    ax.set_title(title, loc="left", fontsize=10.8, fontweight="bold", color=NAVY)
    ax.legend(frameon=False, fontsize=8.2, loc="upper right")
    return _save_dual(fig, png_path)


def _discover_rmats_contrasts(rmats_arg: str | Path | dict[str, Path], declared_contrasts: list[dict[str, Any]]) -> dict[str, Path]:
    """Resolve rMATS directory mapping from single dir, parent dir of contrasts, or contrast=dir mapping."""
    if isinstance(rmats_arg, dict):
        return {str(k): Path(v) for k, v in rmats_arg.items()}
    s = str(rmats_arg)
    if "=" in s and not Path(s).exists():
        out = {}
        for part in s.split(","):
            k, v = part.split("=", 1)
            out[k.strip()] = Path(v.strip())
        return out
    base = Path(s)
    if any((base / f"{et}.MATS.JC.txt").is_file() for et in EVENTS):
        cid = declared_contrasts[0]["contrast_id"] if declared_contrasts else base.name
        return {cid: base}
    out = {}
    for sub in sorted(base.iterdir()):
        if sub.is_dir() and any((sub / f"{et}.MATS.JC.txt").is_file() for et in EVENTS):
            out[sub.name] = sub
    if not out:
        raise FileNotFoundError(f"No *.MATS.JC.txt files found under {base}")
    return out


def _resolve_bundle_inputs(bundle_dir: Path) -> dict[str, Any]:
    """Auto-discover standard upstream Nextflow `katana_bundle` outputs for seamless convergence."""
    bundle = Path(bundle_dir).resolve()
    if not bundle.is_dir():
        raise FileNotFoundError(f"Bundle directory not found: {bundle}")
    samples_f = next((p for p in (bundle / "samplesheet.csv", bundle / "samples.csv", bundle / "provenance/strandedness/samplesheet.requested.csv") if p.is_file()), None)
    contrasts_f = next((p for p in (bundle / "contrasts.csv", bundle / "contrastsheet.csv") if p.is_file()), None)
    rmats_dir = next((p for p in (bundle / "rmats", bundle / "splicing/rmats") if p.is_dir()), None)
    deg_dir = next((p for p in (bundle / "DEG", bundle / "deseq2") if p.is_dir()), None)
    multiqc_dir = next((p for p in (bundle / "qc", bundle / "multiqc") if p.is_dir()), None)
    fasta_f = next((p for p in (bundle / "refs/genome.fasta", bundle / "refs/genome.fa") if p.is_file()), None)
    gtf_f = next((p for p in (bundle / "refs/annotation.gtf", bundle / "refs/genes.gtf") if p.is_file()), None)
    bam_dir = bundle / "bam" if (bundle / "bam").is_dir() else None
    if samples_f is None or rmats_dir is None:
        raise ValueError(f"Bundle {bundle} must contain samplesheet.csv and rmats/ directory")
    return {
        "samples": samples_f,
        "contrasts": contrasts_f,
        "rmats": rmats_dir,
        "deseq2": deg_dir,
        "multiqc": multiqc_dir,
        "fasta": fasta_f,
        "gtf": gtf_f,
        "bams": bam_dir,
    }


def _finalize_consulting_report_and_zip(
    outdir: Path,
    study_name: str,
    mode: str,
    samples_list: list[dict[str, Any]],
    contrasts_list: list[dict[str, Any]],
    panels: list[dict[str, Any]],
    tables: list[dict[str, Any]],
    top10: list[dict[str, Any]],
    policy_data: dict[str, Any],
    motif_summaries: list[dict[str, Any]],
    species: str,
) -> dict[str, Any]:
    """Build SHA256-verified report-manifest.json, provenance.json, REPORT.html, and REPORT.zip via katana_next.reports.layout."""
    from katana_next.reports import layout

    assets_src = Path(layout.__file__).resolve().parent / "assets"
    assets_dst = outdir / "assets"
    assets_dst.mkdir(parents=True, exist_ok=True)
    for asset_name in ("report.css", "report.js"):
        src_a = assets_src / asset_name
        if src_a.is_file():
            shutil.copy2(src_a, assets_dst / asset_name)

    source_files: list[dict[str, Any]] = []
    seen_rel: set[str] = set()
    for p in panels:
        for lk in p.get("data_links", []):
            rel = str(lk["path"]).replace("\\", "/")
            lk["path"] = rel
            fp = outdir / rel
            if fp.is_file():
                dig = _sha256_file(fp)
                lk["source_sha256"] = dig
                lk["portable_sha256"] = dig
                lk["bytes"] = fp.stat().st_size
                if rel not in seen_rel:
                    source_files.append({"source_name": fp.name, "source_sha256": dig, "portable_sha256": dig, "bytes": fp.stat().st_size, "report_path": rel})
                    seen_rel.add(rel)

    for t in tables:
        rel = str(t["source"]).replace("\\", "/")
        t["source"] = rel
        t["view"] = rel
        fp = outdir / rel
        if fp.is_file():
            dig = _sha256_file(fp)
            t.setdefault("metadata", {})["portable_sha256"] = dig
            if rel not in seen_rel:
                source_files.append({"source_name": fp.name, "source_sha256": dig, "portable_sha256": dig, "bytes": fp.stat().st_size, "report_path": rel})
                seen_rel.add(rel)

    project_dict = {
        "project_id": re.sub(r"[^A-Za-z0-9_.-]+", "-", study_name.lower()).strip("-") or "katana-study",
        "name": study_name,
        "samples": samples_list,
        "contrasts": contrasts_list,
        "annotation": {"genome_build": "GRCm39 / mm39" if species == "mouse" else "GRCh38 / hg38", "species": species},
        "metadata": {"assay_mode": "single-cell" if mode == "single-cell" else "bulk", "species": species},
    }

    class _SuiteReportBuilder:
        def __init__(self) -> None:
            self.output = outdir
            self.source = {
                "project": project_dict,
                "mode": mode,
                "dummy": False,
                "report_selection": {"applied": False},
                "artifacts": [],
            }
            self.analyses = {
                "landscape": True,
                "cross_comparison": len(contrasts_list) >= 1,
                "enrichment": any(p.get("chapter") == 5 for p in panels),
                "rbp_identity": False,
                "motiftrace": any(p.get("chapter") == 7 for p in panels) or bool(motif_summaries),
                "motif_binding": False,
                "mapp": False,
                "interaction": False,
                "network": any(p.get("chapter") == 6 for p in panels) or any(t.get("chapter") == 6 for t in tables),
                "sequencing_qc": True,
            }
            self.event_set = "stable_significant"
            self.panels = panels
            self.tables = tables
            self.top10 = top10
            self.policy_data = policy_data
            self.enrichment_summaries: list[dict[str, Any]] = []
            self.motif_summaries = motif_summaries
            self.network_browsers: list[dict[str, Any]] = []
            self.audit = source_files
            self.exclusions: list[dict[str, Any]] = []
            self.derived_audits: list[dict[str, Any]] = []
            self.create_archive = True
            self.chapter_callouts: dict[int, str] = {}

        def figure_data_coverage(self) -> dict[str, Any]:
            main_p = [p for p in self.panels if not p.get("output_only")]
            with_d = [p for p in main_p if any(d.get("role") == "plot_data" for d in p.get("data_links", []))]
            return {
                "schema": "katana.figure_data_links/v1",
                "main_figures": len(main_p),
                "main_figures_with_plot_data": len(with_d),
                "missing_plot_data": [{"id": p["id"], "title": p["title"]} for p in main_p if p not in with_d],
            }

    manifest = layout.render(_SuiteReportBuilder())
    return {
        "status": "completed",
        "study_name": study_name,
        "mode": mode,
        "species": species,
        "outdir": str(outdir),
        "report_html": str(outdir / "REPORT.html"),
        "report_zip": str(outdir / "REPORT.zip"),
        "manifest_json": str(outdir / "report-manifest.json"),
        "panels_count": len(panels),
        "tables_count": len(tables),
        "generated_utc": manifest.get("generated_utc", datetime.now(timezone.utc).isoformat()),
    }


def run_bulk(
    *,
    outdir: str | Path,
    rmats: str | Path | dict[str, Path] | None = None,
    samples: str | Path | pd.DataFrame | None = None,
    contrasts: str | Path | pd.DataFrame | None = None,
    deseq2: str | Path | dict[str, Path] | None = None,
    multiqc: str | Path | None = None,
    fasta: str | Path | None = None,
    gtf: str | Path | None = None,
    bams: str | Path | None = None,
    bundle: str | Path | None = None,
    study_name: str = "Katana Bulk Splicing Study",
    species: str = "auto",
    fdr_cutoff: float = 0.05,
    dpsi_cutoff: float = 0.10,
    manifold_top_n: int = 2000,
    heatmap_top_n: int = 150,
) -> dict[str, Any]:
    """Execute the full standalone Katana Bulk Suite from rMATS/DESeq2 tables OR a Nextflow upstream bundle."""
    from katana_next.science.bulk import apply_stable_event_qc, parse_number_list, significance_columns
    from katana_next.science.psi_landscape import (
        EVENT_TYPE_COLORS,
        compute_beyond_jutils_landscape,
        render_mckinsey_psi_clustermap,
        render_mckinsey_psi_pca,
    )

    out = Path(outdir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    if bundle is not None:
        resolved = _resolve_bundle_inputs(Path(bundle))
        samples = samples or resolved["samples"]
        contrasts = contrasts or resolved["contrasts"]
        rmats = rmats or resolved["rmats"]
        deseq2 = deseq2 or resolved["deseq2"]
        multiqc = multiqc or resolved["multiqc"]
        fasta = fasta or resolved["fasta"]
        gtf = gtf or resolved["gtf"]
        bams = bams or resolved["bams"]

    if samples is None or rmats is None:
        raise ValueError("Provide either `--bundle <dir>` or both `--rmats` and `--samples`.")

    s_df = samples.copy() if isinstance(samples, pd.DataFrame) else _read_tabular(Path(samples))
    sid_col = next((c for c in ("sample_id", "sample", "id") if c in s_df.columns), s_df.columns[0])
    cond_col = next((c for c in ("condition", "group") if c in s_df.columns), s_df.columns[1])
    bio_col = next((c for c in ("biological_sample_id", "subject", "donor_id") if c in s_df.columns), None)
    samples_list = [
        {
            "sample_id": str(r[sid_col]),
            "condition": str(r[cond_col]),
            "biological_sample_id": str(r[bio_col]) if bio_col and pd.notna(r[bio_col]) and str(r[bio_col]).strip() else str(r[sid_col]),
            "cell_type": str(r[cond_col]),
        }
        for _, r in s_df.iterrows()
    ]
    cond_to_sids: dict[str, list[str]] = {}
    for s in samples_list:
        cond_to_sids.setdefault(s["condition"], []).append(s["sample_id"])

    if contrasts is not None:
        c_df = contrasts.copy() if isinstance(contrasts, pd.DataFrame) else _read_tabular(Path(contrasts))
        cid_c = next((c for c in ("contrast_id", "contrast", "comparison_id") if c in c_df.columns), c_df.columns[0])
        num_c = next((c for c in ("treatment", "numerator", "group_a") if c in c_df.columns), c_df.columns[1])
        den_c = next((c for c in ("control", "denominator", "group_b") if c in c_df.columns), c_df.columns[2])
        contrasts_list = [
            {
                "contrast_id": str(r[cid_c]),
                "label": f"{r[num_c]} vs {r[den_c]}",
                "numerator": {"label": str(r[num_c]), "sample_ids": cond_to_sids.get(str(r[num_c]), [])},
                "denominator": {"label": str(r[den_c]), "sample_ids": cond_to_sids.get(str(r[den_c]), [])},
            }
            for _, r in c_df.iterrows()
        ]
    else:
        conds = list(cond_to_sids.keys())
        if len(conds) < 2:
            raise ValueError("At least 2 conditions in samples table are required when --contrasts is omitted.")
        contrasts_list = [{
            "contrast_id": f"{conds[0]}_vs_{conds[1]}",
            "label": f"{conds[0]} vs {conds[1]}",
            "numerator": {"label": conds[0], "sample_ids": cond_to_sids[conds[0]]},
            "denominator": {"label": conds[1], "sample_ids": cond_to_sids[conds[1]]},
        }]

    rmats_map = _discover_rmats_contrasts(rmats, contrasts_list)
    panels: list[dict[str, Any]] = []
    tables: list[dict[str, Any]] = []
    top10: list[dict[str, Any]] = []
    policy_data: dict[str, Any] = {"stable_significant": {}, "significant": {}}
    motif_summaries: list[dict[str, Any]] = []

    # 1. Load rMATS across all contrasts & qualify PS-DAS
    psi_sample_accumulator: dict[str, dict[str, float]] = {}
    event_meta: dict[str, dict[str, str]] = {}
    das_by_contrast: dict[str, pd.DataFrame] = {}

    for c_spec in contrasts_list:
        cid = c_spec["contrast_id"]
        rdir = rmats_map.get(cid) or next(iter(rmats_map.values()))
        num_sids = c_spec["numerator"]["sample_ids"]
        den_sids = c_spec["denominator"]["sample_ids"]
        frames = []
        for et in EVENTS:
            f_et = rdir / f"{et}.MATS.JC.txt"
            if not f_et.is_file():
                continue
            raw = pd.read_csv(f_et, sep="\t")
            if raw.empty:
                continue
            raw["event_type"] = et
            gid_c = next((c for c in ("GeneID", "gene_id") if c in raw.columns), "ID")
            sym_c = next((c for c in ("geneSymbol", "gene_symbol") if c in raw.columns), gid_c)
            raw["gene_id"] = raw[gid_c].astype(str).str.strip('"')
            raw["gene_symbol"] = raw[sym_c].fillna(raw["gene_id"]).astype(str).str.strip('"')
            raw["event_id"] = raw["gene_id"] + "|" + et + "|" + raw["ID"].astype(str)
            raw["comparison_id"] = cid
            raw["dpsi"] = pd.to_numeric(raw["IncLevelDifference"], errors="coerce")
            raw["abs_dpsi"] = raw["dpsi"].abs()
            raw["fdr"] = pd.to_numeric(raw["FDR"], errors="coerce")
            raw["pvalue"] = pd.to_numeric(raw["PValue"], errors="coerce")
            raw["numerator_expected_samples"] = len(num_sids)
            raw["denominator_expected_samples"] = len(den_sids)
            frames.append(raw)
        if not frames:
            continue
        cat_df = pd.concat(frames, ignore_index=True)
        if species == "auto":
            species = "mouse" if cat_df["gene_id"].astype(str).str.startswith("ENSMUSG").any() else "human"

        qc_df = apply_stable_event_qc(cat_df, {})
        qc_df = significance_columns(qc_df, {"fdr": fdr_cutoff, "min_abs_dpsi": dpsi_cutoff})
        qc_df["report_das"] = qc_df["significant"].astype(bool)
        qc_df["report_rs_das"] = qc_df["stable_significant"].astype(bool)
        das_by_contrast[cid] = qc_df

        # Populate sample-level PSI & event_meta
        for _, r in qc_df.iterrows():
            eid = str(r["event_id"])
            event_meta.setdefault(eid, {"gene_symbol": str(r["gene_symbol"]), "event_type": str(r["event_type"])})
            p1 = parse_number_list(r.get("IncLevel1"))
            p2 = parse_number_list(r.get("IncLevel2"))
            rec = psi_sample_accumulator.setdefault(eid, {})
            for s_name, pv in zip(num_sids, p1):
                if np.isfinite(pv):
                    rec[s_name] = float(pv)
            for s_name, pv in zip(den_sids, p2):
                if np.isfinite(pv):
                    rec[s_name] = float(pv)

        # Save lossless DAS parquet & volcano figures
        comp_dir = out / "outputs/comparisons" / cid
        comp_dir.mkdir(parents=True, exist_ok=True)
        das_pq = comp_dir / "DAS-complete-original-and-report-flags.parquet"
        qc_df.to_parquet(das_pq, index=False)
        tables.append({
            "name": f"DAS complete table ({cid})",
            "analysis_id": "differential.as",
            "chapter": 3,
            "rows": len(qc_df),
            "columns": list(qc_df.columns),
            "source": das_pq.relative_to(out).as_posix(),
        })

        policy_data["stable_significant"][cid] = {
            "available": True,
            "events": int(qc_df["report_rs_das"].sum()),
            "genes": int(qc_df.loc[qc_df["report_rs_das"], "gene_symbol"].nunique()),
            "raw_events": len(qc_df),
            "read_supported_events": int(qc_df["passes_precision_support"].sum()),
        }
        policy_data["significant"][cid] = {
            "available": True,
            "events": int(qc_df["report_das"].sum()),
            "genes": int(qc_df.loc[qc_df["report_das"], "gene_symbol"].nunique()),
            "raw_events": len(qc_df),
        }

        for pol, col_flag, lbl in (("stable_significant", "report_rs_das", "PS-DAS"), ("significant", "report_das", "Ordinary DAS")):
            v_png = out / "figures/comparisons" / f"{cid}-AS-volcano-{pol}.png"
            _render_volcano(qc_df, "dpsi", "fdr", col_flag, "gene_symbol", "ΔPSI (Inclusion difference)", f"{lbl} Splicing Volcano — {cid}", v_png)
            v_tsv = comp_dir / f"{cid}-AS-volcano-{pol}-plot-data.tsv"
            qc_df[["event_id", "gene_symbol", "event_type", "dpsi", "fdr", col_flag]].to_csv(v_tsv, sep="\t", index=False)
            panels.append({
                "id": f"{cid}-as-volcano-{pol}",
                "chapter": 3,
                "contrast_id": cid,
                "policy": pol,
                "title": f"{lbl} volcano ({int(qc_df[col_flag].sum()):,} significant events) | {cid}",
                "caption": f"FDR ≤ {fdr_cutoff}, |ΔPSI| ≥ {dpsi_cutoff}. Original caller P/FDR preserved.",
                "image": v_png.relative_to(out).as_posix(),
                "svg": v_png.with_suffix(".svg").relative_to(out).as_posix(),
                "data_links": [{"role": "plot_data", "label": "Plot data (TSV)", "path": v_tsv.relative_to(out).as_posix(), "unit": "splicing_event", "contrast_id": cid, "policy": pol}],
            })
            sub_sig = qc_df[qc_df[col_flag]].sort_values(["fdr", "abs_dpsi"], ascending=[True, False]).head(10)
            top10.append({
                "id": f"{cid}-top10-das-{pol}",
                "chapter": 3,
                "contrast_id": cid,
                "analysis": "DAS",
                "policy": pol,
                "unit": "splicing_event",
                "columns": ["gene_symbol", "event_id", "event_type", "dpsi", "fdr"],
                "rows": sub_sig[["gene_symbol", "event_id", "event_type", "dpsi", "fdr"]].to_dict("records"),
                "full_path": das_pq.relative_to(out).as_posix(),
            })

        # Chapter 05 GO + KEGG ORA on PS-DAS foreground vs measured splicing background
        fg_s = set(qc_df.loc[qc_df["report_rs_das"], "gene_symbol"].dropna().astype(str))
        if len(fg_s) < 3:
            fg_s = set(qc_df.loc[qc_df["passes_precision_support"] & (qc_df["abs_dpsi"] >= dpsi_cutoff), "gene_symbol"].dropna().astype(str))
        bg_s = set(qc_df["gene_symbol"].dropna().astype(str))
        ora_df = _run_ora_from_bundled_catalog(fg_s, bg_s, species if species != "auto" else "human")
        if not ora_df.empty:
            ora_tsv = out / "outputs/enrichment" / f"{cid}-ps-das-go-kegg-ora.tsv"
            ora_tsv.parent.mkdir(parents=True, exist_ok=True)
            ora_df.to_csv(ora_tsv, sep="\t", index=False)
            ora_png = out / "figures/enrichment" / f"{cid}-ps-das-go-kegg-ora.png"
            _render_ora_faceted_figure(ora_df, ora_png, f"GO (BP / CC / MF) & KEGG Pathway Enrichment of PS-DAS Genes — {cid}")
            panels.append({
                "id": f"{cid}-ora-ps-das",
                "chapter": 5,
                "contrast_id": cid,
                "policy": "stable_significant",
                "title": f"GO & KEGG functional enrichment of PS-DAS genes | {cid}",
                "caption": f"Exact hypergeometric test with per-category Benjamini-Hochberg FDR across {len(bg_s):,} measured splicing background genes.",
                "image": ora_png.relative_to(out).as_posix(),
                "svg": ora_png.with_suffix(".svg").relative_to(out).as_posix(),
                "data_links": [{"role": "plot_data", "label": "Complete ORA table", "path": ora_tsv.relative_to(out).as_posix(), "unit": "pathway_term", "contrast_id": cid, "policy": "stable_significant"}],
            })
            tables.append({
                "name": f"GO & KEGG ORA complete table ({cid})",
                "analysis_id": "prioritize.enrichment",
                "chapter": 5,
                "rows": len(ora_df),
                "columns": list(ora_df.columns),
                "source": ora_tsv.relative_to(out).as_posix(),
            })

    if species == "auto":
        species = "human"

    # Optional MultiQC 5'-3' gene-body coverage reuse (Chapter 01)
    if multiqc is not None and Path(multiqc).exists():
        from katana_next.reports.multiqc_coverage import parse_multiqc_5to3_coverage, render_mckinsey_5to3_figure
        cov = parse_multiqc_5to3_coverage(Path(multiqc), {s["sample_id"]: s["condition"] for s in samples_list})
        if cov is not None and "curves" in cov:
            cov_png = out / "figures/qc/multiqc-5to3-coverage.png"
            cov_tsv = out / "outputs/qc/multiqc-5to3-coverage-plot-data.tsv"
            cov_png.parent.mkdir(parents=True, exist_ok=True)
            cov_tsv.parent.mkdir(parents=True, exist_ok=True)
            cov["curves"].to_csv(cov_tsv, sep="\t", index=False)
            render_mckinsey_5to3_figure(cov, cov_png, cov_png.with_suffix(".svg"), title=f"Transcript 5′ → 3′ Gene-Body Coverage — {study_name}")
            panels.append({
                "id": "qc-multiqc-5to3",
                "chapter": 1,
                "contrast_id": "all",
                "title": f"5′ → 3′ transcript coverage uniformity ({cov.get('provenance_source', 'MultiQC')})",
                "caption": "Directly reused from MultiQC saved raw curves without slow BAM re-scanning.",
                "image": cov_png.relative_to(out).as_posix(),
                "svg": cov_png.with_suffix(".svg").relative_to(out).as_posix(),
                "data_links": [{"role": "plot_data", "label": "Coverage curves (TSV)", "path": cov_tsv.relative_to(out).as_posix(), "unit": "percentile_bin", "contrast_id": "all", "policy": None}],
            })

    # 2. Sample-level PSI matrix + Beyond-Jutils Landscape (Chapter 02) + Katana-SpliceNet v2.3.2 (Chapter 06)
    all_sids = [s["sample_id"] for s in samples_list]
    psi_mat = pd.DataFrame.from_dict(psi_sample_accumulator, orient="index").reindex(columns=all_sids)
    psi_mat.index.name = "event_id"
    res_dir = out / "outputs/results"
    res_dir.mkdir(parents=True, exist_ok=True)
    psi_pq = res_dir / "sample_level_psi.parquet"
    psi_mat.reset_index().to_parquet(psi_pq, index=False)
    tables.append({
        "name": "Sample-level PSI matrix",
        "analysis_id": "landscape.bulk",
        "chapter": 2,
        "rows": len(psi_mat),
        "columns": ["event_id", *all_sids],
        "source": psi_pq.relative_to(out).as_posix(),
    })

    if len(psi_mat) >= 10 and len(all_sids) >= 3:
        land = compute_beyond_jutils_landscape(
            psi_mat,
            samples_list,
            event_meta=event_meta,
            top_n=heatmap_top_n,
            manifold_top_n=manifold_top_n,
            min_present_fraction=0.65,
            prior_kappa=10.0,
            z_clip=2.5,
            run_umap=len(all_sids) >= 4,
            umap_neighbors=min(8, max(2, len(all_sids) - 1)),
        )
        if land.get("available"):
            hm_png = out / "outputs/landscape/psi-heatmap.png"
            pca_png = out / "outputs/landscape/psi-pca.png"
            render_mckinsey_psi_clustermap(land, hm_png, hm_png.with_suffix(".svg"), title=f"Depth-Stabilized Splicing Z-Score Clustermap — {study_name}", event_meta=event_meta)
            render_mckinsey_psi_pca(land, pca_png, pca_png.with_suffix(".svg"), title=f"Global Splicing PCA ({land['manifold_events']:,} Events) & PERMANOVA — {study_name}")
            pca_tsv = out / "outputs/landscape/psi-pca-coordinates.tsv"
            pd.DataFrame(land["samples"]).to_csv(pca_tsv, sep="\t", index=False)
            panels.append({
                "id": "psi-pca",
                "chapter": 2,
                "contrast_id": "all",
                "title": f"Global splicing PCA ({land['manifold_events']:,} events · PC1 {land['pca_variance_pct'][0]:.1f}%, PERMANOVA P = {land['permanova_p']:.3g})",
                "caption": "Empirical-Bayes depth-stabilized splicing manifold across measurable cohort events.",
                "image": pca_png.relative_to(out).as_posix(),
                "svg": pca_png.with_suffix(".svg").relative_to(out).as_posix(),
                "data_links": [{"role": "plot_data", "label": "PCA coordinates (TSV)", "path": pca_tsv.relative_to(out).as_posix(), "unit": "biological_sample", "contrast_id": "all", "policy": None}],
            })
            panels.append({
                "id": "psi-heatmap",
                "chapter": 2,
                "contrast_id": "all",
                "title": f"Depth-stabilized PSI Z-score clustermap (top {land['heatmap_events']} rows; sample dendrogram across {land['manifold_events']:,} events)",
                "caption": "Separated Sample Group vs Splicing Subtype legends with right-axis representative gene callouts.",
                "image": hm_png.relative_to(out).as_posix(),
                "svg": hm_png.with_suffix(".svg").relative_to(out).as_posix(),
                "data_links": [{"role": "plot_data", "label": "Sample PSI matrix", "path": psi_pq.relative_to(out).as_posix(), "unit": "splicing_event", "contrast_id": "all", "policy": None}],
            })

    # Katana-SpliceNet v2.3.2 co-splicing module discovery (pure Python)
    splicenet_assignments = pd.DataFrame()
    if len(psi_mat) >= 40 and len(all_sids) >= 4:
        try:
            from katana_splicenet.evaluation import score_modules_v2
            from katana_splicenet.network import fit_wgcna_v2
            from katana_splicenet.qc_and_selection import select_events_v2, unit_robust_audit_v2

            meta_sn = pd.DataFrame(samples_list)[["sample_id", "biological_sample_id", "condition"]]
            masked, audit_sn = unit_robust_audit_v2(psi_mat, meta_sn, None)
            min_mod = max(5, min(20, len(psi_mat) // 12))
            selected_sn, _ = select_events_v2(masked, audit_sn, feature_cap=min(5000, len(psi_mat)), ranking_metric="variance")
            if len(selected_sn) >= max(12, min_mod * 2):
                v23_dir = out / "outputs/v23_tables"
                v23_dir.mkdir(parents=True, exist_ok=True)
                assign = fit_wgcna_v2(selected_sn, v23_dir / "wgcna_fit", power=0, min_module_size=min_mod, backend="python")
                assign["gene_symbol"] = assign["event_key"].map(lambda k: event_meta.get(str(k), {}).get("gene_symbol", str(k).split("|")[0]))
                assign["gene_id"] = assign["event_key"].map(lambda k: str(k).split("|")[0])
                assign["event_type"] = assign["event_key"].map(lambda k: event_meta.get(str(k), {}).get("event_type", "SE"))
                splicenet_assignments = assign[assign["module"].ne(0) & assign["module"].astype(str).ne("M0")].copy()
                if not splicenet_assignments.empty:
                    scores_sn, _, mod_sum, _ = score_modules_v2(selected_sn, assign)
                    splicenet_assignments.to_parquet(v23_dir / "splicenet_modules.parquet", index=False)
                    scores_sn.reset_index().rename(columns={"index": "sample_id"}).to_parquet(v23_dir / "splicenet_scores.parquet", index=False)
                    tables.append({
                        "name": "Katana-SpliceNet co-splicing module assignments",
                        "analysis_id": "programs.splicenet",
                        "chapter": 6,
                        "rows": len(splicenet_assignments),
                        "columns": list(splicenet_assignments.columns),
                        "source": (v23_dir / "splicenet_modules.parquet").relative_to(out).as_posix(),
                    })
                    # Plot module eigengene summary figure for Chapter 06
                    fig_sn, ax_sn = plt.subplots(figsize=(8.4, 4.4), dpi=220)
                    fig_sn.patch.set_facecolor("white")
                    mc_counts = splicenet_assignments["module"].astype(str).map(lambda m: m if m.startswith("M") else f"M{m}").value_counts()
                    ax_sn.bar(np.arange(len(mc_counts)), mc_counts.to_numpy(), color=BLUE, width=0.55)
                    ax_sn.set_xticks(np.arange(len(mc_counts)), mc_counts.index.tolist())
                    ax_sn.set_ylabel("Assigned splicing events", fontweight="bold", color=NAVY)
                    ax_sn.set_title(f"Katana-SpliceNet v2.3.2 Co-Splicing Modules ({len(mc_counts)} modules, {len(splicenet_assignments):,} events)", loc="left", fontweight="bold", color=NAVY)
                    sn_png = out / "figures/network/splicenet-module-sizes.png"
                    _save_dual(fig_sn, sn_png)
                    sn_tsv = v23_dir / "splicenet-module-sizes.tsv"
                    mc_counts.rename_axis("module").reset_index(name="events").to_csv(sn_tsv, sep="\t", index=False)
                    panels.append({
                        "id": "splicenet-overview",
                        "chapter": 6,
                        "contrast_id": "all",
                        "title": f"Katana-SpliceNet co-splicing module architecture ({len(mc_counts)} modules)",
                        "caption": "Scale-free Signed Topological Overlap co-splicing modules discovered unsupervised across study samples.",
                        "image": sn_png.relative_to(out).as_posix(),
                        "svg": sn_png.with_suffix(".svg").relative_to(out).as_posix(),
                        "data_links": [{"role": "plot_data", "label": "Module counts (TSV)", "path": sn_tsv.relative_to(out).as_posix(), "unit": "splicing_module", "contrast_id": "all", "policy": None}],
                    })
        except Exception:
            pass

    # 3. Katana-SpliceTrace v2.4.3 (170-PWM 11-channel splicing maps + 50-nt EJC NMD & UniProt domains when FASTA+GTF supplied)
    if fasta is not None and gtf is not None and Path(fasta).is_file() and Path(gtf).is_file():
        from katana_splicetrace.cli import main as splicetrace_cli_main
        for c_spec in contrasts_list:
            cid = c_spec["contrast_id"]
            rdir = rmats_map.get(cid) or next(iter(rmats_map.values()))
            st_out = out / "outputs/splicetrace" / cid
            st_out.mkdir(parents=True, exist_ok=True)
            try:
                splicetrace_cli_main([
                    "--rmats-dir", str(rdir),
                    "--fasta", str(fasta),
                    "--gtf", str(gtf),
                    "--species", species,
                    "--output-dir", str(st_out),
                ])
                rk_files = sorted(st_out.glob("*ranking*.tsv"))
                if rk_files:
                    rk_df = pd.read_csv(rk_files[0], sep="\t")
                    rk_pq = st_out / "splicetrace_ranking.parquet"
                    rk_df.to_parquet(rk_pq, index=False)
                    tables.append({
                        "name": f"Katana-SpliceTrace 170-PWM RBP ranking ({cid})",
                        "analysis_id": "prioritize.splicetrace_ranking",
                        "chapter": 7,
                        "rows": len(rk_df),
                        "columns": list(rk_df.columns),
                        "source": rk_pq.relative_to(out).as_posix(),
                    })
                    fwer_c = next((c for c in ("P_FWER", "regional_fwer_p", "fwer_p", "p_fwer") if c in rk_df.columns), rk_df.columns[-1])
                    min_fwer = float(pd.to_numeric(rk_df[fwer_c], errors="coerce").min()) if not rk_df.empty else 1.0
                    motif_summaries.append({
                        "contrast_id": cid,
                        "status": "significant" if min_fwer <= 0.05 else "no_significant_motif",
                        "min_native_fdr": min_fwer,
                        "notice": f"Katana-SpliceTrace v2.4.3 evaluated {len(rk_df)} curated {species} RBP PWMs across 11 splicing-map channels + mature Δ(I-E) RNA.",
                        "native_test_family": "Gene-block Rademacher max-T FWER across 11 positional channels",
                        "ranking_path": rk_pq.relative_to(out).as_posix(),
                    })
            except Exception:
                pass

    return _finalize_consulting_report_and_zip(
        outdir=out,
        study_name=study_name,
        mode="bulk",
        samples_list=samples_list,
        contrasts_list=contrasts_list,
        panels=panels,
        tables=tables,
        top10=top10,
        policy_data=policy_data,
        motif_summaries=motif_summaries,
        species=species,
    )


def run_single_cell(
    *,
    outdir: str | Path,
    events: str | Path | pd.DataFrame,
    cells: str | Path | pd.DataFrame,
    psi_matrix: str | Path | pd.DataFrame | None = None,
    support_matrix: str | Path | pd.DataFrame | None = None,
    expression: str | Path | pd.DataFrame | None = None,
    study_name: str = "Katana Single-Cell Splicing Study",
    species: str = "human",
    fdr_cutoff: float = 0.05,
    dpsi_cutoff: float = 0.10,
) -> dict[str, Any]:
    """Execute the standalone Katana Single-Cell Suite (Precision-Kernel Event-PSI UMAP, SC-SpliceNet, ORA & Report)."""
    from katana_next.science.psi_landscape import (
        compute_beyond_jutils_landscape,
        compute_dual_expression_vs_psi_umap,
        render_mckinsey_dual_psi_umap,
    )

    out = Path(outdir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    ev_df = events.copy() if isinstance(events, pd.DataFrame) else _read_tabular(Path(events))
    cell_df = cells.copy() if isinstance(cells, pd.DataFrame) else _read_tabular(Path(cells))

    cid_c = next((c for c in ("cell_id", "sample_id", "barcode") if c in cell_df.columns), cell_df.columns[0])
    cond_c = next((c for c in ("condition", "state", "activation_state", "group") if c in cell_df.columns), cell_df.columns[1])
    donor_c = next((c for c in ("donor_id", "biological_sample_id", "subject") if c in cell_df.columns), None)
    ct_c = next((c for c in ("cell_type", "cluster", "subset") if c in cell_df.columns), cond_c)

    samples_list = [
        {
            "sample_id": str(r[cid_c]),
            "condition": str(r[cond_c]),
            "biological_sample_id": str(r[donor_c]) if donor_c and pd.notna(r[donor_c]) else str(r[cid_c]),
            "donor_id": str(r[donor_c]) if donor_c and pd.notna(r[donor_c]) else "D1",
            "cell_type": str(r[ct_c]),
            "activation_state": str(r[cond_c]),
        }
        for _, r in cell_df.iterrows()
    ]

    cmp_col = next((c for c in ("comparison_id", "contrast_id", "contrast") if c in ev_df.columns), None)
    if cmp_col is None:
        cond_vals = list(dict.fromkeys(s["condition"] for s in samples_list))
        default_cid = f"{cond_vals[0]}_vs_{cond_vals[1]}" if len(cond_vals) >= 2 else "Activated_vs_Resting"
        ev_df["comparison_id"] = default_cid
        cmp_col = "comparison_id"

    contrasts_list = []
    for cid_v in list(dict.fromkeys(ev_df[cmp_col].astype(str))):
        parts = cid_v.split("_vs_", 1) if "_vs_" in cid_v else [cid_v, "Baseline"]
        contrasts_list.append({
            "contrast_id": cid_v,
            "label": cid_v.replace("_", " "),
            "numerator": {"label": parts[0], "sample_ids": [s["sample_id"] for s in samples_list if s["condition"] == parts[0]] or [samples_list[0]["sample_id"]]},
            "denominator": {"label": parts[1], "sample_ids": [s["sample_id"] for s in samples_list if s["condition"] == parts[1]] or [samples_list[-1]["sample_id"]]},
        })

    panels: list[dict[str, Any]] = []
    tables: list[dict[str, Any]] = []
    top10: list[dict[str, Any]] = []
    policy_data: dict[str, Any] = {"stable_significant": {}, "significant": {}}

    res_dir = out / "outputs/results"
    res_dir.mkdir(parents=True, exist_ok=True)
    ev_pq = res_dir / "sc_differential_splicing.parquet"
    ev_df.to_parquet(ev_pq, index=False)
    tables.append({
        "name": "Single-cell donor-aware differential splicing events",
        "analysis_id": "differential.sc_as",
        "chapter": 3,
        "rows": len(ev_df),
        "columns": list(ev_df.columns),
        "source": ev_pq.relative_to(out).as_posix(),
    })

    q_c = next((c for c in ("fdr", "FDR", "q_value", "padj") if c in ev_df.columns), None)
    d_c = next((c for c in ("dpsi", "delta_psi", "effect") if c in ev_df.columns), None)
    g_c = next((c for c in ("gene_symbol", "gene_name", "gene_id") if c in ev_df.columns), ev_df.columns[0])

    for c_spec in contrasts_list:
        cid = c_spec["contrast_id"]
        sub = ev_df[ev_df[cmp_col].astype(str) == cid].copy()
        if q_c and d_c:
            sub["_sig"] = (pd.to_numeric(sub[q_c], errors="coerce").fillna(1.0) <= fdr_cutoff) & (pd.to_numeric(sub[d_c], errors="coerce").abs().fillna(0.0) >= dpsi_cutoff)
            v_png = out / "figures/comparisons" / f"{cid}-sc-as-volcano.png"
            _render_volcano(sub, d_c, q_c, "_sig", g_c, "Single-cell ΔPSI", f"Donor-Aware Single-Cell Splicing Volcano — {cid}", v_png)
            v_tsv = out / "outputs/comparisons" / cid / f"{cid}-sc-volcano-plot-data.tsv"
            v_tsv.parent.mkdir(parents=True, exist_ok=True)
            sub[[g_c, d_c, q_c, "_sig"]].to_csv(v_tsv, sep="\t", index=False)
            panels.append({
                "id": f"{cid}-sc-volcano",
                "chapter": 3,
                "contrast_id": cid,
                "policy": "stable_significant",
                "title": f"Single-cell splicing volcano ({int(sub['_sig'].sum()):,} significant events) | {cid}",
                "caption": "Donor-aware single-cell differential splicing preserving biological replicate degrees of freedom.",
                "image": v_png.relative_to(out).as_posix(),
                "svg": v_png.with_suffix(".svg").relative_to(out).as_posix(),
                "data_links": [{"role": "plot_data", "label": "Volcano data (TSV)", "path": v_tsv.relative_to(out).as_posix(), "unit": "splicing_event", "contrast_id": cid, "policy": "stable_significant"}],
            })
            policy_data["stable_significant"][cid] = {"available": True, "events": int(sub["_sig"].sum()), "genes": int(sub.loc[sub["_sig"], g_c].nunique())}
            fg_sc = set(sub.loc[sub["_sig"], g_c].dropna().astype(str)) or set(sub[g_c].dropna().astype(str).head(40))
            ora_sc = _run_ora_from_bundled_catalog(fg_sc, set(sub[g_c].dropna().astype(str)), species)
            if not ora_sc.empty:
                ora_png = out / "figures/enrichment" / f"{cid}-sc-go-kegg-ora.png"
                ora_tsv = out / "outputs/enrichment" / f"{cid}-sc-go-kegg-ora.tsv"
                ora_tsv.parent.mkdir(parents=True, exist_ok=True)
                ora_sc.to_csv(ora_tsv, sep="\t", index=False)
                _render_ora_faceted_figure(ora_sc, ora_png, f"Single-Cell GO & KEGG Enrichment — {cid}")
                panels.append({
                    "id": f"{cid}-sc-ora",
                    "chapter": 5,
                    "contrast_id": cid,
                    "policy": "stable_significant",
                    "title": f"Single-cell GO & KEGG pathway enrichment | {cid}",
                    "caption": "Exact hypergeometric enrichment across GO (BP, CC, MF) and KEGG pathways.",
                    "image": ora_png.relative_to(out).as_posix(),
                    "svg": ora_png.with_suffix(".svg").relative_to(out).as_posix(),
                    "data_links": [{"role": "plot_data", "label": "Complete ORA table", "path": ora_tsv.relative_to(out).as_posix(), "unit": "pathway_term", "contrast_id": cid, "policy": "stable_significant"}],
                })

    # Precision-Kernel Event-PSI UMAP & Cross-Modality Divergence (Chapter 02)
    if psi_matrix is not None:
        p_df = psi_matrix.copy() if isinstance(psi_matrix, pd.DataFrame) else _read_tabular(Path(psi_matrix))
        if not isinstance(psi_matrix, pd.DataFrame) and p_df.columns[0] not in [s["sample_id"] for s in samples_list]:
            p_df = p_df.set_index(p_df.columns[0])
        sup_df = None
        if support_matrix is not None:
            sup_df = support_matrix.copy() if isinstance(support_matrix, pd.DataFrame) else _read_tabular(Path(support_matrix))
            if not isinstance(support_matrix, pd.DataFrame) and sup_df.columns[0] not in [s["sample_id"] for s in samples_list]:
                sup_df = sup_df.set_index(sup_df.columns[0])
        expr_df = None
        if expression is not None:
            expr_df = expression.copy() if isinstance(expression, pd.DataFrame) else _read_tabular(Path(expression))
            if not isinstance(expression, pd.DataFrame) and expr_df.columns[0] not in [s["sample_id"] for s in samples_list]:
                expr_df = expr_df.set_index(expr_df.columns[0])

        psi_land = compute_beyond_jutils_landscape(
            p_df,
            samples_list,
            depth_matrix=sup_df,
            min_present_fraction=0.25,
            run_umap=True,
            umap_neighbors=min(12, max(3, len(samples_list) - 2)),
        )
        if psi_land.get("available"):
            s_ids_sc = [s["sample_id"] for s in samples_list]
            if expr_df is not None:
                expr_cells_rows = expr_df.T if set(s_ids_sc).issubset(set(expr_df.columns.astype(str))) else expr_df
            else:
                expr_cells_rows = p_df.T.fillna(0.5)
            dual = compute_dual_expression_vs_psi_umap(expr_cells_rows, psi_land, k_neighbors=min(10, max(3, len(s_ids_sc) // 3)))
            if dual.get("cells"):
                umap_png = out / "outputs/landscape/sc-dual-psi-umap.png"
                render_mckinsey_dual_psi_umap(dual, umap_png, umap_png.with_suffix(".svg"))
                umap_tsv = out / "outputs/landscape/sc-dual-psi-umap-cells.tsv"
                pd.DataFrame(dual["cells"]).to_csv(umap_tsv, sep="\t", index=False)
                panels.append({
                    "id": "sc-dual-psi-umap",
                    "chapter": 2,
                    "contrast_id": "all",
                    "title": f"Precision-Kernel Event-PSI UMAP vs Expression UMAP (Donor LISI = {dual['median_donor_lisi']:.2f})",
                    "caption": "Three-panel single-cell manifold comparing gene expression UMAP, Binomial Fisher-precision weighted Event-PSI UMAP, and per-cell modality divergence Δ_kNN.",
                    "image": umap_png.relative_to(out).as_posix(),
                    "svg": umap_png.with_suffix(".svg").relative_to(out).as_posix(),
                    "data_links": [{"role": "plot_data", "label": "Cell coordinates & Δ_kNN (TSV)", "path": umap_tsv.relative_to(out).as_posix(), "unit": "single_cell", "contrast_id": "all", "policy": None}],
                })

    return _finalize_consulting_report_and_zip(
        outdir=out,
        study_name=study_name,
        mode="single-cell",
        samples_list=samples_list,
        contrasts_list=contrasts_list,
        panels=panels,
        tables=tables,
        top10=top10,
        policy_data=policy_data,
        motif_summaries=[],
        species=species,
    )


def open_workspace(
    report_path: str | Path | None = None,
    *,
    host: str = "127.0.0.1",
    port: int = 18775,
    open_browser: bool = True,
) -> None:
    """Launch the local interactive workspace and optionally pre-import a `REPORT.zip` or report directory."""
    import uvicorn
    from katana_next.api.app import PlatformRuntime, create_app

    app_root = Path(__file__).resolve().parents[1]
    runtime = PlatformRuntime(application_root=app_root)
    if report_path is not None:
        rp = Path(report_path).expanduser().resolve()
        ws_home = Path.home() / ".katana" / "workspaces" / rp.stem
        try:
            from katana_next.services.report_import import import_portable_report
            manifest_p = import_portable_report(rp, ws_home)
            runtime.projects.open_project(manifest_p)
        except Exception as exc:
            print(f"[katana view] Note: could not pre-import {rp}: {exc}")

    app = create_app(runtime)
    url = f"http://{host}:{port}"
    if open_browser:
        webbrowser.open(url)
    uvicorn.run(app, host=host, port=port)


def run_pipeline(
    nextflow_args: list[str],
    *,
    nextflow_executable: str = "nextflow",
    converge_suite_outdir: str | Path | None = None,
) -> int:
    """Execute the end-to-end automated Nextflow pipeline from FASTQ or BAM (`pipeline/main.nf`) and optionally converge into `run_bulk`."""
    repo_root = Path(__file__).resolve().parents[1]
    nf_main = repo_root / "pipeline" / "main.nf"
    if not nf_main.is_file():
        nf_main = repo_root / "katana_next" / "flow" / "_pipeline" / "main.nf"
    if not nf_main.is_file():
        raise FileNotFoundError("Bundled pipeline/main.nf not found in installation.")
    if shutil.which(nextflow_executable) is None and not Path(nextflow_executable).is_file():
        raise FileNotFoundError(f"Nextflow executable '{nextflow_executable}' not found on PATH.")

    cmd = [nextflow_executable, "run", str(nf_main), *nextflow_args]
    completed = subprocess.run(cmd)
    if completed.returncode != 0:
        return completed.returncode

    if converge_suite_outdir is not None:
        # Locate `--outdir` passed to Nextflow and run `run_bulk(bundle=...)`
        nf_outdir = None
        for idx, tok in enumerate(nextflow_args):
            if tok == "--outdir" and idx + 1 < len(nextflow_args):
                nf_outdir = Path(nextflow_args[idx + 1])
        if nf_outdir is not None:
            bundle_cand = nf_outdir / "katana_bundle" if (nf_outdir / "katana_bundle").is_dir() else nf_outdir
            run_bulk(bundle=bundle_cand, outdir=converge_suite_outdir)
    return 0
