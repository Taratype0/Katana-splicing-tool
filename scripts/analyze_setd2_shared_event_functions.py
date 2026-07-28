from __future__ import annotations

import argparse
import html
import json
import re
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.analyze_setd2_matched_background_mut_coverage import _load_rmats_comparison


EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
CORE_ORDER = (
    "CLK3-like",
    "PlaB-like",
    "MUT-CLK3-PlaB triple",
)
CORE_COLORS = {
    "CLK3-like": "#D55E00",
    "PlaB-like": "#0072B2",
    "MUT-CLK3-PlaB triple": "#6A51A3",
}
EVENT_COLORS = {
    "SE": "#4C78A8",
    "RI": "#F58518",
    "A3SS": "#54A24B",
    "A5SS": "#E45756",
    "MXE": "#9D6CC1",
}
GPROFILER_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
RC = {
    "font.family": "Arial",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.titlesize": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def _normalize_gene_id(value: object) -> str:
    return str(value).strip().strip('"').split(".", 1)[0]


def _parse_attributes(value: str) -> dict[str, str]:
    return {
        key: item
        for key, item in re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)"', value)
    }


def _save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_core_catalog(catalog: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for context in ("DP", "noDP"):
        subset = catalog.loc[catalog["context"].eq(context)].copy()
        clk = subset.loc[
            subset["CLK3_covers_MUT"].fillna(False)
            & subset["CLK3_same_direction"].fillna(False)
        ].copy()
        clk["core_set"] = "CLK3-like"
        clk["inhibitor"] = "CLK3"
        clk["inhibitor_FDR"] = clk["CLK3_FDR"]
        clk["inhibitor_dPSI"] = clk["CLK3_dPSI"]
        rows.append(clk)

        pla = subset.loc[
            subset["Pla_covers_MUT"].fillna(False)
            & subset["Pla_same_direction"].fillna(False)
        ].copy()
        pla["core_set"] = "PlaB-like"
        pla["inhibitor"] = "PlaB"
        pla["inhibitor_FDR"] = pla["Pla_FDR"]
        pla["inhibitor_dPSI"] = pla["Pla_dPSI"]
        rows.append(pla)

        triple = subset.loc[
            subset["CLK3_covers_MUT"].fillna(False)
            & subset["Pla_covers_MUT"].fillna(False)
            & subset["CLK3_same_direction"].fillna(False)
            & subset["Pla_same_direction"].fillna(False)
        ].copy()
        triple["core_set"] = "MUT-CLK3-PlaB triple"
        triple["inhibitor"] = "CLK3+PlaB"
        triple["inhibitor_FDR"] = triple[["CLK3_FDR", "Pla_FDR"]].max(axis=1)
        triple["inhibitor_dPSI"] = triple[["CLK3_dPSI", "Pla_dPSI"]].mean(axis=1)
        rows.append(triple)

    keep = [
        "context",
        "core_set",
        "event_key",
        "gene_id",
        "gene_symbol",
        "event_type",
        "chr",
        "strand",
        "event_coordinates",
        "MUT_FDR",
        "MUT_dPSI",
        "inhibitor",
        "inhibitor_FDR",
        "inhibitor_dPSI",
        "CLK3_FDR",
        "CLK3_dPSI",
        "Pla_FDR",
        "Pla_dPSI",
    ]
    output = pd.concat(rows, ignore_index=True)
    output["same_direction"] = (
        np.sign(output["MUT_dPSI"]) == np.sign(output["inhibitor_dPSI"])
    )
    output["abs_delta_dPSI"] = (
        output["MUT_dPSI"] - output["inhibitor_dPSI"]
    ).abs()
    return output.loc[:, [*keep, "same_direction", "abs_delta_dPSI"]]


def _event_regions(row: pd.Series) -> list[tuple[str, int, int]]:
    parts = str(row["event_key"]).split("|")
    event_type = str(row["event_type"])
    coordinates = [int(float(value)) for value in parts[4:]]
    if event_type == "SE":
        return [("skipped_exon", coordinates[0] + 1, coordinates[1])]
    if event_type == "RI":
        return [("retained_intron", coordinates[0] + 1, coordinates[1])]
    if event_type == "MXE":
        return [
            ("mutually_exclusive_exon_1", coordinates[0] + 1, coordinates[1]),
            ("mutually_exclusive_exon_2", coordinates[2] + 1, coordinates[3]),
        ]
    if event_type in {"A3SS", "A5SS"}:
        long_start, long_end, short_start, short_end = coordinates[:4]
        starts = sorted((long_start + 1, short_start + 1))
        ends = sorted((long_end, short_end))
        variable_start = starts[0]
        variable_end = starts[1] - 1 if starts[0] != starts[1] else ends[0] + 1
        if variable_end < variable_start:
            variable_start, variable_end = min(ends), max(ends)
        return [(f"{event_type}_variable_region", variable_start, variable_end)]
    return []


def annotate_events_with_gtf(events: pd.DataFrame, gtf_path: Path) -> pd.DataFrame:
    target_genes = {_normalize_gene_id(value) for value in events["gene_id"]}
    features: dict[str, list[dict[str, object]]] = defaultdict(list)
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"exon", "CDS", "UTR"}:
                continue
            attributes = _parse_attributes(fields[8])
            gene_id = _normalize_gene_id(attributes.get("gene_id", ""))
            if gene_id not in target_genes:
                continue
            features[gene_id].append(
                {
                    "chrom": fields[0],
                    "feature": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "transcript_id": attributes.get("transcript_id", ""),
                    "transcript_name": attributes.get("transcript_name", ""),
                    "exon_number": attributes.get("exon_number", ""),
                    "exon_id": attributes.get("exon_id", ""),
                    "gene_type": attributes.get("gene_type", ""),
                    "transcript_type": attributes.get("transcript_type", ""),
                }
            )

    rows: list[dict[str, object]] = []
    for row in events.drop_duplicates(["context", "core_set", "event_key"]).itertuples(index=False):
        series = pd.Series(row._asdict())
        regions = _event_regions(series)
        gene_features = features.get(_normalize_gene_id(row.gene_id), [])
        overlaps: list[dict[str, object]] = []
        for region_label, start, end in regions:
            for feature in gene_features:
                if feature["chrom"] != row.chr or feature["strand"] != row.strand:
                    continue
                if int(feature["end"]) < start or int(feature["start"]) > end:
                    continue
                overlaps.append({**feature, "region_label": region_label})
        feature_types = {str(item["feature"]) for item in overlaps}
        if row.event_type == "RI":
            consequence_class = "retained intron; isoform consequence unresolved"
        elif "CDS" in feature_types:
            consequence_class = "CDS-overlapping event"
        elif "UTR" in feature_types:
            consequence_class = "UTR-overlapping event"
        elif "exon" in feature_types:
            consequence_class = "annotated exon; coding consequence unresolved"
        else:
            consequence_class = "no overlapping GTF exon/CDS/UTR feature"
        rows.append(
            {
                "context": row.context,
                "core_set": row.core_set,
                "event_key": row.event_key,
                "gene_id": row.gene_id,
                "gene_symbol": row.gene_symbol,
                "event_type": row.event_type,
                "chr": row.chr,
                "strand": row.strand,
                "affected_regions": ";".join(
                    f"{label}:{row.chr}:{start}-{end}" for label, start, end in regions
                ),
                "affected_exon_numbers": ";".join(
                    sorted({str(item["exon_number"]) for item in overlaps if item["exon_number"]})
                ),
                "affected_exon_ids": ";".join(
                    sorted({str(item["exon_id"]) for item in overlaps if item["exon_id"]})
                ),
                "affected_transcript_ids": ";".join(
                    sorted({str(item["transcript_id"]) for item in overlaps if item["transcript_id"]})
                ),
                "affected_transcript_names": ";".join(
                    sorted({str(item["transcript_name"]) for item in overlaps if item["transcript_name"]})
                ),
                "n_affected_transcripts": len(
                    {str(item["transcript_id"]) for item in overlaps if item["transcript_id"]}
                ),
                "overlaps_CDS": "CDS" in feature_types,
                "overlaps_UTR": "UTR" in feature_types,
                "overlaps_annotated_exon": "exon" in feature_types,
                "gene_types": ";".join(
                    sorted({str(item["gene_type"]) for item in overlaps if item["gene_type"]})
                ),
                "predicted_consequence_class": consequence_class,
                "consequence_scope": (
                    "structural annotation only; reading frame, NMD and protein-domain effects "
                    "require explicit isoform reconstruction"
                ),
            }
        )
    return pd.DataFrame(rows)


def _gprofiler_query(
    name: str,
    genes: list[str],
    background: list[str],
    timeout_seconds: int = 120,
) -> tuple[pd.DataFrame, dict[str, object]]:
    query = sorted({str(gene) for gene in genes if str(gene).strip()})
    domain = sorted({str(gene) for gene in background if str(gene).strip()})
    payload = {
        "organism": "mmusculus",
        "query": query,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": domain,
        "no_evidences": False,
    }
    request = Request(
        GPROFILER_ENDPOINT,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urlopen(request, timeout=timeout_seconds) as response:  # noqa: S310
        result = json.loads(response.read().decode("utf-8"))
    rows = []
    for item in result.get("result", []):
        evidence = item.get("intersections", [])
        intersection_genes = [
            gene
            for gene, gene_evidence in zip(query, evidence, strict=False)
            if gene_evidence
        ]
        rows.append(
            {
                "analysis_set": name,
                "term_id": item.get("native", ""),
                "term_name": item.get("name", ""),
                "source": item.get("source", ""),
                "FDR": item.get("p_value"),
                "intersection_size": item.get("intersection_size"),
                "term_size": item.get("term_size"),
                "query_size": item.get("query_size"),
                "effective_domain_size": item.get("effective_domain_size"),
                "intersection_genes": ";".join(intersection_genes),
                "method": "g:Profiler g:GOSt FDR with a custom testable-gene background",
            }
        )
    provider = {
        "analysis_set": name,
        "provider": "g:Profiler g:GOSt",
        "endpoint": GPROFILER_ENDPOINT,
        "data_version": result.get("meta", {}).get("data_version", {}),
        "query_genes": len(query),
        "background_genes": len(domain),
    }
    return pd.DataFrame(rows), provider


def build_enrichment(
    core: pd.DataFrame,
    rmats: dict[str, pd.DataFrame],
    comparisons: dict[str, dict[str, str]],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames: list[pd.DataFrame] = []
    providers: list[dict[str, object]] = []
    for context, mapping in comparisons.items():
        mut = rmats[mapping["MUT"]]
        for core_set, inhibitor_key in (
            ("CLK3-like", "CLK3"),
            ("PlaB-like", "PlaB"),
            ("MUT-CLK3-PlaB triple", "triple"),
        ):
            subset = core.loc[
                core["context"].eq(context) & core["core_set"].eq(core_set)
            ]
            if inhibitor_key == "CLK3":
                tested = mut.merge(
                    rmats[mapping["CLK3"]][["event_key"]],
                    on="event_key",
                    how="inner",
                )
            elif inhibitor_key == "PlaB":
                tested = mut.merge(
                    rmats[mapping["PlaB"]][["event_key"]],
                    on="event_key",
                    how="inner",
                )
            else:
                tested = (
                    mut.merge(
                        rmats[mapping["CLK3"]][["event_key"]],
                        on="event_key",
                        how="inner",
                    )
                    .merge(
                        rmats[mapping["PlaB"]][["event_key"]],
                        on="event_key",
                        how="inner",
                    )
                )
            genes = sorted(set(subset["gene_symbol"].dropna().astype(str)))
            background = sorted(set(tested["gene_symbol"].dropna().astype(str)))
            analysis_set = f"{context} | {core_set}"
            try:
                result, provider = _gprofiler_query(analysis_set, genes, background)
                result.insert(0, "context", context)
                result.insert(1, "core_set", core_set)
                frames.append(result)
                providers.append({**provider, "status": "finished"})
            except Exception as exc:
                providers.append(
                    {
                        "analysis_set": analysis_set,
                        "provider": "g:Profiler g:GOSt",
                        "status": "failed",
                        "error": str(exc),
                        "query_genes": len(genes),
                        "background_genes": len(background),
                    }
                )
    enrichment = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return enrichment, providers


def _plot_core_composition(core: pd.DataFrame, path: Path) -> None:
    counts = (
        core.groupby(["context", "core_set", "event_type"], observed=True)["event_key"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(columns=EVENT_TYPES, fill_value=0)
    )
    labels = [f"{context}\n{core_set}" for context, core_set in counts.index]
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(8.4, 4.5))
        bottom = np.zeros(len(counts))
        for event_type in EVENT_TYPES:
            values = counts[event_type].to_numpy()
            ax.bar(
                np.arange(len(counts)),
                values,
                bottom=bottom,
                label=event_type,
                color=EVENT_COLORS[event_type],
                width=0.72,
            )
            bottom += values
        for index, total in enumerate(bottom):
            ax.text(index, total + max(bottom) * 0.015, f"n={int(total)}", ha="center", fontsize=7)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=25, ha="right")
        ax.set_ylabel("Direction-concordant shared significant AS events")
        ax.set_title("Composition of inhibitor-like Setd2 MUT splicing cores", loc="left", fontweight="bold")
        ax.legend(title="Event type", frameon=False, ncol=5, loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.25, top=0.88)
        _save_figure(fig, path)


def _plot_annotation_summary(annotation: pd.DataFrame, path: Path) -> None:
    summary = (
        annotation.groupby(["context", "core_set", "predicted_consequence_class"], observed=True)[
            "event_key"
        ]
        .nunique()
        .rename("n_events")
        .reset_index()
    )
    classes = summary["predicted_consequence_class"].drop_duplicates().tolist()
    labels = [
        f"{context}\n{core_set}"
        for context, core_set in summary[["context", "core_set"]].drop_duplicates().itertuples(index=False)
    ]
    keys = list(summary[["context", "core_set"]].drop_duplicates().itertuples(index=False, name=None))
    colors = ["#2C7FB8", "#7FCDBB", "#FDAE6B", "#BDBDBD", "#756BB1"]
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(8.6, 4.7))
        bottom = np.zeros(len(keys))
        for class_index, consequence in enumerate(classes):
            values = []
            for context, core_set in keys:
                match = summary.loc[
                    summary["context"].eq(context)
                    & summary["core_set"].eq(core_set)
                    & summary["predicted_consequence_class"].eq(consequence),
                    "n_events",
                ]
                values.append(int(match.iloc[0]) if not match.empty else 0)
            ax.bar(
                np.arange(len(keys)),
                values,
                bottom=bottom,
                label=consequence,
                color=colors[class_index % len(colors)],
                width=0.72,
            )
            bottom += np.asarray(values)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=25, ha="right")
        ax.set_ylabel("Annotated shared AS events")
        ax.set_title("GTF-based structural context of shared AS events", loc="left", fontweight="bold")
        ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
        fig.subplots_adjust(left=0.10, right=0.68, bottom=0.25, top=0.88)
        _save_figure(fig, path)


def _plot_enrichment(enrichment: pd.DataFrame, path: Path) -> None:
    if enrichment.empty:
        with plt.rc_context(RC):
            fig, ax = plt.subplots(figsize=(8, 2.5))
            ax.text(0.5, 0.5, "No enrichment result returned", ha="center", va="center")
            ax.axis("off")
            _save_figure(fig, path)
        return
    working = enrichment.copy()
    working["FDR"] = pd.to_numeric(working["FDR"], errors="coerce")
    working["score"] = -np.log10(working["FDR"].clip(lower=1e-300))
    significant = working.loc[working["FDR"].le(0.05)].nsmallest(12, "FDR").copy()
    best_by_core = (
        working.dropna(subset=["FDR"])
        .sort_values("FDR")
        .drop_duplicates("analysis_set")
        .sort_values(["context", "core_set"])
        .reset_index(drop=True)
    )
    with plt.rc_context(RC):
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(13.0, 4.8),
            gridspec_kw={"width_ratios": [1.45, 1.0]},
        )
        ax_terms, ax_core = axes

        if significant.empty:
            ax_terms.text(
                0.5,
                0.55,
                "No GO / KEGG / Reactome term passed FDR <= 0.05",
                ha="center",
                va="center",
                fontweight="bold",
            )
            ax_terms.text(
                0.5,
                0.42,
                "Shared-event cores do not show a broad functional enrichment signal.",
                ha="center",
                va="center",
                color="#555555",
            )
            ax_terms.axis("off")
        else:
            significant = significant.sort_values("FDR", ascending=False).reset_index(drop=True)
            y = np.arange(len(significant))
            colors = [
                CORE_COLORS.get(str(core_set), "#777777")
                for core_set in significant["core_set"]
            ]
            sizes = 45 + 16 * np.sqrt(
                pd.to_numeric(significant["intersection_size"], errors="coerce").fillna(0)
            )
            ax_terms.scatter(
                significant["score"],
                y,
                s=sizes,
                c=colors,
                edgecolor="#333333",
                linewidth=0.5,
            )
            ax_terms.set_yticks(y)
            ax_terms.set_yticklabels(
                [
                    f"{row.source} | {row.term_name}\n"
                    f"drivers: {row.intersection_genes}"
                    for row in significant.itertuples(index=False)
                ]
            )
            ax_terms.set_xlabel("-log10 enrichment FDR")
            ax_terms.set_title("FDR-significant functional terms", loc="left", fontweight="bold")
            ax_terms.spines[["top", "right"]].set_visible(False)

        threshold = -np.log10(0.05)
        y_core = np.arange(len(best_by_core))
        core_colors = [
            CORE_COLORS.get(str(row.core_set), "#777777")
            if row.FDR <= 0.05
            else "#C7C7C7"
            for row in best_by_core.itertuples(index=False)
        ]
        ax_core.barh(y_core, best_by_core["score"], color=core_colors, height=0.68)
        ax_core.axvline(threshold, color="#555555", linestyle="--", linewidth=0.9)
        ax_core.set_yticks(y_core)
        ax_core.set_yticklabels(
            [
                f"{row.context} | {row.core_set}"
                for row in best_by_core.itertuples(index=False)
            ]
        )
        ax_core.invert_yaxis()
        ax_core.set_xlabel("Best -log10 enrichment FDR")
        ax_core.set_title("Best result per shared-event core", loc="left", fontweight="bold")
        ax_core.text(
            threshold,
            len(best_by_core) - 0.15,
            "FDR 0.05",
            ha="left",
            va="bottom",
            fontsize=7,
            color="#555555",
        )
        ax_core.spines[["top", "right"]].set_visible(False)

        fig.suptitle(
            "Functional enrichment of direction-concordant shared AS-event genes",
            x=0.06,
            ha="left",
            fontweight="bold",
        )
        fig.subplots_adjust(left=0.29, right=0.98, bottom=0.16, top=0.82, wspace=0.48)
        _save_figure(fig, path)


def _summary_tables(
    core: pd.DataFrame,
    annotation: pd.DataFrame,
    enrichment: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    core_summary = (
        core.groupby(["context", "core_set"], observed=True)
        .agg(
            n_events=("event_key", "nunique"),
            n_genes=("gene_symbol", "nunique"),
            median_MUT_abs_dPSI=("MUT_dPSI", lambda values: values.abs().median()),
            median_inhibitor_abs_dPSI=("inhibitor_dPSI", lambda values: values.abs().median()),
            median_abs_delta_dPSI=("abs_delta_dPSI", "median"),
        )
        .reset_index()
    )
    annotation_summary = (
        annotation.groupby(["context", "core_set", "predicted_consequence_class"], observed=True)[
            "event_key"
        ]
        .nunique()
        .rename("n_events")
        .reset_index()
    )
    if not enrichment.empty:
        significant = (
            enrichment.loc[pd.to_numeric(enrichment["FDR"], errors="coerce").le(0.05)]
            .groupby(["context", "core_set"], observed=True)
            .agg(
                n_FDR_significant_terms=("term_id", "nunique"),
                best_enrichment_FDR=("FDR", "min"),
                best_term=("term_name", "first"),
            )
            .reset_index()
        )
        core_summary = core_summary.merge(
            significant,
            on=["context", "core_set"],
            how="left",
        )
    core_summary["n_FDR_significant_terms"] = (
        core_summary.get("n_FDR_significant_terms", 0).fillna(0).astype(int)
    )
    return core_summary, annotation_summary


def _build_report(
    core_summary: pd.DataFrame,
    enrichment: pd.DataFrame,
    output: Path,
    gtf_path: Path,
) -> None:
    lines = [
        "# Functional follow-up of shared significant AS events",
        "",
        "## Scope",
        "",
        "This report analyzes exact AS events that are significant in Setd2 MUT and an inhibitor "
        "contrast in the same DP background and have the same standardized dPSI sign. It is a "
        "functional-prioritization analysis, not evidence that a pathway or RBP caused the event.",
        "",
        "Significance filter: rMATS FDR <= 0.05 and abs(standardized dPSI) >= 0.10.",
        "",
        "## Core sizes",
        "",
        core_summary.to_markdown(index=False),
        "",
        "## Functional enrichment interpretation",
        "",
        "GO/KEGG/Reactome enrichment uses unique host genes. The background is restricted to genes "
        "with exact AS events testable in the relevant matched MUT and inhibitor contrasts. A GO "
        "term describes the host genes; it does not by itself prove the functional consequence of "
        "the splice event.",
        "",
    ]
    significant = enrichment.loc[
        pd.to_numeric(enrichment.get("FDR"), errors="coerce").le(0.05)
    ] if not enrichment.empty else pd.DataFrame()
    if significant.empty:
        lines.append(
            "No core produced an FDR-significant GO/KEGG/Reactome term. The shared events remain "
            "valid molecular overlap sets, but they do not support one broad pathway-level claim."
        )
    else:
        for (context, core_set), subset in significant.groupby(
            ["context", "core_set"], observed=True
        ):
            lines.extend(
                [
                    f"### {context} | {core_set}",
                    "",
                    subset.nsmallest(10, "FDR")[
                        [
                            "term_id",
                            "term_name",
                            "source",
                            "FDR",
                            "intersection_size",
                            "intersection_genes",
                        ]
                    ].to_markdown(index=False),
                    "",
                ]
            )
    lines.extend(
        [
            "## Structural annotation boundary",
            "",
            f"Events were intersected with `{gtf_path}`. CDS/UTR/exon overlap is descriptive. "
            "Reading-frame change, NMD and protein-domain gain/loss are not inferred without "
            "explicit transcript reconstruction and isoform-switch support.",
            "",
            "## Figures",
            "",
            "1. `plots/01_core_event_type_composition.png`",
            "2. `plots/02_functional_enrichment.png`",
            "3. `plots/03_gtf_structural_context.png`",
        ]
    )
    markdown = "\n".join(lines)
    (output / "SHARED_EVENT_FUNCTIONAL_REPORT.md").write_text(markdown, encoding="utf-8")
    body = markdown.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    body = re.sub(r"^# (.+)$", r"<h1>\1</h1>", body, flags=re.MULTILINE)
    body = re.sub(r"^## (.+)$", r"<h2>\1</h2>", body, flags=re.MULTILINE)
    body = re.sub(r"^### (.+)$", r"<h3>\1</h3>", body, flags=re.MULTILINE)
    body = body.replace("\n", "<br>\n")
    images = "".join(
        f'<figure><img src="plots/{name}" style="max-width:100%;height:auto"></figure>'
        for name in (
            "01_core_event_type_composition.png",
            "02_functional_enrichment.png",
            "03_gtf_structural_context.png",
        )
    )
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Shared AS-event functional follow-up</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1180px;margin:32px auto;color:#17202a;line-height:1.5}}
h1,h2,h3{{color:#102a43}} code{{background:#eef2f5;padding:2px 5px}} figure{{margin:28px 0}}
</style></head><body>{body}{images}</body></html>"""
    (output / "SHARED_EVENT_FUNCTIONAL_REPORT.html").write_text(document, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Functional follow-up of matched shared AS events.")
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path(
            r"E:\Project\RNASEQ\Result\Setd2_results\katana_output"
            r"\09_cosplice_network\06l_matched_background_MUT_coverage"
            r"\matched_background_MUT_event_catalog.tsv"
        ),
    )
    parser.add_argument(
        "--rmats-root",
        type=Path,
        default=Path(
            r"E:\Project\RNASEQ\Result\Setd2_results"
            r"\Setd2_Splicing_genomebam\star_salmon\rmats"
        ),
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=Path(r"E:\Project\RNASEQ\gencode.vM34.annotation.gtf"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            r"E:\Project\RNASEQ\Result\Setd2_results\katana_output"
            r"\09_cosplice_network\06m_shared_event_functional_followup"
        ),
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plot_dir = args.output / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    comparisons = {
        "DP": {
            "MUT": "MUT_DP-WT_DP",
            "CLK3": "WT_DP_CLK3-WT_DP",
            "PlaB": "WT_DP_Pla-WT_DP",
        },
        "noDP": {
            "MUT": "MUT_NO_DP-WT_NO_DP",
            "CLK3": "WT_NO_DP_CLK3-WT_NO_DP",
            "PlaB": "WT_NO_DP_Pla-WT_NO_DP",
        },
    }
    catalog = pd.read_csv(args.catalog, sep="\t", low_memory=False)
    core = build_core_catalog(catalog)
    annotation = annotate_events_with_gtf(core, args.gtf)

    rmats = {
        name: _load_rmats_comparison(args.rmats_root, name, "JCEC")
        for mapping in comparisons.values()
        for name in mapping.values()
    }
    enrichment, providers = build_enrichment(core, rmats, comparisons)
    core_summary, annotation_summary = _summary_tables(core, annotation, enrichment)

    core.to_csv(args.output / "shared_direction_concordant_event_catalog.tsv", sep="\t", index=False)
    annotation.to_csv(args.output / "shared_event_gtf_annotation.tsv", sep="\t", index=False)
    enrichment.to_csv(args.output / "shared_event_go_kegg_reactome_enrichment.tsv", sep="\t", index=False)
    core_summary.to_csv(args.output / "shared_event_core_summary.tsv", sep="\t", index=False)
    annotation_summary.to_csv(
        args.output / "shared_event_gtf_annotation_summary.tsv",
        sep="\t",
        index=False,
    )
    (args.output / "enrichment_provider.json").write_text(
        json.dumps(providers, indent=2),
        encoding="utf-8",
    )
    (args.output / "parameters.json").write_text(
        json.dumps(
            {
                "source_catalog": str(args.catalog),
                "rmats_root": str(args.rmats_root),
                "gtf": str(args.gtf),
                "rmats_mode": "JCEC",
                "rMATS_FDR_cutoff": 0.05,
                "abs_standardized_dPSI_cutoff": 0.10,
                "core_requirement": "exact shared significant event with the same standardized dPSI sign",
                "enrichment_background": "genes with exact AS events testable in the matched contrasts",
                "protein_consequence_boundary": (
                    "No frame/NMD/domain claim without explicit isoform reconstruction"
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    _plot_core_composition(core, plot_dir / "01_core_event_type_composition.png")
    _plot_enrichment(enrichment, plot_dir / "02_functional_enrichment.png")
    _plot_annotation_summary(annotation, plot_dir / "03_gtf_structural_context.png")
    _build_report(core_summary, enrichment, args.output, args.gtf)

    print(core_summary.to_string(index=False))
    if not enrichment.empty:
        significant = enrichment.loc[pd.to_numeric(enrichment["FDR"], errors="coerce").le(0.05)]
        print(f"\nFDR-significant enrichment terms: {len(significant)}")
        if not significant.empty:
            print(
                significant.nsmallest(20, "FDR")[
                    ["context", "core_set", "term_id", "term_name", "source", "FDR"]
                ].to_string(index=False)
            )
    print(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
