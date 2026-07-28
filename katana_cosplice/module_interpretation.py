from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom, pearsonr

from katana_cosplice.config import ModuleInterpretationConfig
from katana_cosplice.functional_enrichment import gprofiler_enrichment
from katana_cosplice.isoform_switch import enrich_isoform_switches, quantify_isoform_switches
from katana_cosplice.motif_scanner import scan_event_motifs
from katana_cosplice.module_annotation import _load_rbp_sets, _read_candidate_support
from katana_cosplice.utils import bh_adjust, normalize_gene_id


@dataclass
class ModuleInterpretationResult:
    hub_events: pd.DataFrame
    hub_genes: pd.DataFrame
    joint_evidence: pd.DataFrame
    rbp_motif_enrichment: pd.DataFrame
    rbp_motif_hits: pd.DataFrame
    rbp_expression_correlation: pd.DataFrame
    rbp_expression_matrix: pd.DataFrame
    gene_set_enrichment: pd.DataFrame
    isoform_switches: pd.DataFrame
    module_isoform_switch_enrichment: pd.DataFrame
    regulator_ranking: pd.DataFrame
    source_audit: pd.DataFrame


def _module_column(assignments: pd.DataFrame) -> str:
    for column in ("consensus_module", "module"):
        if column in assignments:
            return column
    raise ValueError("Module assignments require consensus_module or module.")


def _read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce")


def _gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def _quant_gene_expression(quant_root: Path, gtf: Path) -> pd.DataFrame:
    transcript_to_symbol: dict[str, str] = {}
    with gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attributes = _gtf_attributes(fields[8])
            transcript = attributes.get("transcript_id", "").split(".")[0]
            symbol = attributes.get("gene_name", "")
            if transcript and symbol:
                transcript_to_symbol[transcript] = symbol
    sample_series: dict[str, pd.Series] = {}
    for path in sorted(quant_root.glob("**/quant.sf")):
        sample = path.parent.name
        frame = pd.read_csv(path, sep="\t", usecols=["Name", "TPM"])
        frame["gene_symbol"] = frame["Name"].astype(str).str.split(".").str[0].map(transcript_to_symbol)
        frame["TPM"] = pd.to_numeric(frame["TPM"], errors="coerce")
        sample_series[sample] = frame.dropna(subset=["gene_symbol"]).groupby("gene_symbol", observed=True)["TPM"].sum()
    return pd.DataFrame(sample_series).fillna(0.0)


def _safe_pearson(left: pd.Series, right: pd.Series) -> tuple[float, float, int]:
    aligned = pd.concat([left, right], axis=1).dropna()
    if len(aligned) < 4 or aligned.iloc[:, 0].nunique() < 2 or aligned.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, len(aligned)
    correlation, pvalue = pearsonr(aligned.iloc[:, 0], aligned.iloc[:, 1])
    return float(correlation), float(pvalue), len(aligned)


def _hypergeom(universe: int, selected: int, module_size: int, overlap: int) -> float:
    if universe <= 0 or selected <= 0 or module_size <= 0 or overlap <= 0:
        return np.nan
    return float(hypergeom.sf(overlap - 1, universe, selected, module_size))


def _load_gmt(path: Path | None) -> tuple[dict[str, set[str]], dict[str, object]]:
    if path is None:
        return {}, {"source_type": "gene_sets_gmt", "path": "", "status": "not_configured", "rows": 0}
    if not path.exists():
        return {}, {"source_type": "gene_sets_gmt", "path": str(path), "status": "missing_file", "rows": 0}
    sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = {value.strip().casefold() for value in fields[2:] if value.strip()}
    return sets, {"source_type": "gene_sets_gmt", "path": str(path), "status": "loaded", "rows": len(sets)}


def interpret_modules(
    assignments: pd.DataFrame,
    psi: pd.DataFrame,
    program_scores: pd.DataFrame,
    event_metadata: pd.DataFrame,
    settings: ModuleInterpretationConfig,
    *,
    support_root: Path | None = None,
    rbp_gene_set: Path | None = None,
    quant_root: Path | None = None,
    gtf: Path | None = None,
    genome_fasta: Path | None = None,
    motif_hits_override: Path | None = None,
    isoform_switches_override: Path | None = None,
    rbp_expression_override: Path | None = None,
    sample_metadata: pd.DataFrame | None = None,
    condition_column: str = "condition",
) -> ModuleInterpretationResult:
    module_column = _module_column(assignments)
    membership = assignments[["event_key", module_column, *[
        name for name in ("consensus_membership_confidence", "module_membership", "network_connectivity")
        if name in assignments
    ]]].rename(columns={module_column: "module"})
    membership = membership.loc[membership["module"].astype(str).ne("grey")].drop_duplicates("event_key")
    events = membership.merge(event_metadata, on="event_key", how="left", validate="one_to_one")
    psi = psi.reindex(index=events["event_key"])
    scores = program_scores.copy()
    if not scores.index.isin(psi.columns).all() and "sample_id" in scores:
        scores = scores.set_index("sample_id")
    scores.index = scores.index.astype(str)

    hub_rows: list[dict[str, object]] = []
    for _, event in events.iterrows():
        module = str(event["module"])
        if module not in scores.columns or event["event_key"] not in psi.index:
            correlation, pvalue, n_samples = np.nan, np.nan, 0
        else:
            correlation, pvalue, n_samples = _safe_pearson(
                psi.loc[event["event_key"]].rename_axis("sample_id"),
                pd.to_numeric(scores[module], errors="coerce").rename_axis("sample_id"),
            )
        confidence = event.get("consensus_membership_confidence", event.get("module_membership", np.nan))
        connectivity = event.get("network_connectivity", np.nan)
        hub_rows.append(
            {
                "event_key": event["event_key"],
                "module": module,
                "gene_id": normalize_gene_id(event.get("gene_id")),
                "gene_symbol": str(event.get("gene_symbol", "") or ""),
                "event_type": event.get("event_type", ""),
                "module_membership_kME": correlation,
                "module_membership_pvalue": pvalue,
                "n_samples": n_samples,
                "consensus_confidence": confidence,
                "network_connectivity": connectivity,
            }
        )
    hub_events = pd.DataFrame(hub_rows)
    if not hub_events.empty:
        hub_events["module_membership_FDR"] = bh_adjust(hub_events["module_membership_pvalue"])
        for column in ("module_membership_kME", "consensus_confidence", "network_connectivity"):
            numeric = pd.to_numeric(hub_events[column], errors="coerce")
            if column == "module_membership_kME":
                numeric = numeric.abs()
            hub_events[f"{column}_percentile"] = numeric.groupby(hub_events["module"]).rank(pct=True)
        score_columns = [
            "module_membership_kME_percentile",
            "consensus_confidence_percentile",
            "network_connectivity_percentile",
        ]
        hub_events["hub_event_score"] = hub_events[score_columns].mean(axis=1, skipna=True)
        hub_events["hub_event_rank"] = hub_events.groupby("module")["hub_event_score"].rank(method="min", ascending=False)

    candidate_support, candidate_audit = _read_candidate_support(support_root)
    genes = hub_events.copy()
    genes["gene_key"] = np.where(
        genes["gene_id"].astype(str).ne(""),
        "id:" + genes["gene_id"].astype(str),
        "sym:" + genes["gene_symbol"].astype(str).str.casefold(),
    )
    hub_genes = genes.groupby(["module", "gene_key"], as_index=False).agg(
        gene_id=("gene_id", "first"),
        gene_symbol=("gene_symbol", "first"),
        n_module_events=("event_key", "nunique"),
        event_types=("event_type", lambda value: ";".join(sorted(set(map(str, value))))),
        max_abs_kME=("module_membership_kME", lambda value: pd.to_numeric(value, errors="coerce").abs().max()),
        mean_hub_event_score=("hub_event_score", "mean"),
        max_hub_event_score=("hub_event_score", "max"),
    )
    hub_genes["hub_gene_score"] = (
        hub_genes.groupby("module")["max_abs_kME"].rank(pct=True)
        + hub_genes.groupby("module")["max_hub_event_score"].rank(pct=True)
        + hub_genes.groupby("module")["n_module_events"].rank(pct=True)
    ) / 3.0
    hub_genes["hub_gene_rank"] = hub_genes.groupby("module")["hub_gene_score"].rank(method="min", ascending=False)

    joint = hub_genes.merge(candidate_support, on="gene_key", how="left", suffixes=("", "_support"))
    for flag in ("DEG_significant", "rMATS_significant", "DEXSeq_significant", "DTU_significant"):
        if flag not in joint:
            joint[flag] = False
        else:
            joint[flag] = joint[flag].map(lambda value: bool(value) if pd.notna(value) else False).astype(bool)
    joint["joint_support_count"] = joint[["DEG_significant", "rMATS_significant", "DEXSeq_significant", "DTU_significant"]].sum(axis=1)

    rbp_sets, rbp_audit = _load_rbp_sets(rbp_gene_set)
    rbp_symbols = set(rbp_sets["gene_symbol"].astype(str).str.strip().str.casefold())
    joint["RBP_gene"] = joint["gene_symbol"].astype(str).str.strip().str.casefold().isin(rbp_symbols)

    source_rows = [*candidate_audit.to_dict("records"), *rbp_audit.to_dict("records")]
    eligible_modules = set(
        hub_genes.groupby("module", observed=True)["gene_key"]
        .nunique()
        .loc[lambda values: values.ge(settings.min_module_genes)]
        .index.astype(str)
    )

    motif_rows: list[dict[str, object]] = []
    motif_path = motif_hits_override or settings.rbp_motif_hits
    generated_motif_hits = pd.DataFrame()
    if motif_path is None and settings.rbp_motif_database is not None and genome_fasta is not None:
        if settings.rbp_motif_database.exists() and genome_fasta.exists():
            generated_motif_hits = scan_event_motifs(
                events, genome_fasta, settings.rbp_motif_database,
                window_bp=settings.motif_window_bp,
                score_threshold=settings.motif_score_threshold,
                organism=settings.rbp_motif_organism,
            )
            hits = generated_motif_hits
            motif_path = settings.rbp_motif_database
            source_rows.append({
                "source_type": "rbp_motif_scan", "path": str(settings.rbp_motif_database),
                "status": "generated_from_genome_fasta", "rows": len(hits),
            })
        else:
            hits = pd.DataFrame()
            source_rows.append({
                "source_type": "rbp_motif_scan", "path": str(settings.rbp_motif_database),
                "status": "missing_motif_database_or_genome_fasta", "rows": 0,
            })
    elif motif_path is None:
        hits = pd.DataFrame()
        source_rows.append({"source_type": "rbp_motif_hits", "path": "", "status": "not_configured", "rows": 0})
    elif not motif_path.exists():
        hits = pd.DataFrame()
        source_rows.append({"source_type": "rbp_motif_hits", "path": str(motif_path), "status": "missing_file", "rows": 0})
    else:
        hits = pd.read_csv(motif_path, sep="," if motif_path.suffix.lower() == ".csv" else "\t", low_memory=False)
    all_motif_hits = hits.copy()
    if not hits.empty:
        motif_col = next((name for name in ("RBP", "rbp", "motif", "motif_id") if name in hits), None)
        if "event_key" not in hits or motif_col is None:
            raise ValueError("RBP motif-hit table requires event_key and RBP/rbp/motif/motif_id.")
        hits = hits.drop_duplicates(["event_key", motif_col])
        universe = set(events["event_key"])
        for motif, subset in hits.groupby(motif_col, observed=True):
            selected = set(subset["event_key"]) & universe
            for module, module_events in events.groupby("module", observed=True):
                if str(module) not in eligible_modules:
                    continue
                members = set(module_events["event_key"])
                overlap = members & selected
                motif_rows.append({
                    "module": module,
                    "RBP_or_motif": motif,
                    "event_universe": len(universe),
                    "motif_events_in_universe": len(selected),
                    "module_events": len(members),
                    "overlap_events": len(overlap),
                    "pvalue": _hypergeom(len(universe), len(selected), len(members), len(overlap)),
                    "method": "rMAPS2-inspired event-level positional motif-hit enrichment",
                })
        if generated_motif_hits.empty:
            source_rows.append({"source_type": "rbp_motif_hits", "path": str(motif_path), "status": "loaded", "rows": len(hits)})
    motif_enrichment = pd.DataFrame(
        motif_rows,
        columns=[
            "module", "RBP_or_motif", "event_universe", "motif_events_in_universe",
            "module_events", "overlap_events", "pvalue", "method",
            "FDR", "significant_at_configured_FDR",
        ],
    )
    if not motif_enrichment.empty:
        motif_enrichment["FDR"] = bh_adjust(motif_enrichment["pvalue"])
        motif_enrichment["significant_at_configured_FDR"] = motif_enrichment["FDR"].le(settings.fdr_cutoff)
        motif_enrichment = motif_enrichment.sort_values(["FDR", "pvalue"], na_position="last")

    expression_rows: list[dict[str, object]] = []
    expression_path = settings.rbp_expression_matrix
    expression: pd.DataFrame | None = None
    expression_source = ""
    if rbp_expression_override is not None and rbp_expression_override.exists():
        expression = _read_matrix(rbp_expression_override)
        expression_source = str(rbp_expression_override)
    elif expression_path is not None and expression_path.exists():
        expression = _read_matrix(expression_path)
        expression_source = str(expression_path)
    elif quant_root is not None and gtf is not None and quant_root.exists() and gtf.exists():
        expression = _quant_gene_expression(quant_root, gtf)
        expression_source = f"quant.sf under {quant_root}; transcript-to-gene mapping from {gtf}"
    if expression is None:
        if expression_path is None:
            source_rows.append({"source_type": "rbp_expression_matrix", "path": "", "status": "not_configured_and_quant_gtf_unavailable", "rows": 0})
        else:
            source_rows.append({"source_type": "rbp_expression_matrix", "path": str(expression_path), "status": "missing_file_and_quant_gtf_unavailable", "rows": 0})
    else:
        expression.index = expression.index.astype(str).str.strip().str.casefold()
        expression = expression.groupby(level=0, observed=True).sum()
        common_samples = expression.columns.intersection(scores.index)
        rbp_names = {
            str(value).strip().casefold(): str(value).strip()
            for value in rbp_sets["gene_symbol"].astype(str)
            if str(value).strip()
        }
        for rbp_key in sorted(set(expression.index) & set(rbp_names)):
            rbp = rbp_names[rbp_key]
            for module in scores.columns:
                if str(module) not in eligible_modules:
                    continue
                correlation, pvalue, n_samples = _safe_pearson(
                    expression.loc[rbp_key, common_samples], scores.loc[common_samples, module]
                )
                expression_rows.append({"module": module, "RBP": rbp, "correlation": correlation, "pvalue": pvalue, "n_samples": n_samples})
        source_rows.append({"source_type": "rbp_expression_matrix", "path": expression_source, "status": "loaded", "rows": len(expression)})
    expression_correlation = pd.DataFrame(
        expression_rows,
        columns=[
            "module", "RBP", "correlation", "pvalue", "n_samples", "FDR",
            "significant_at_configured_FDR",
        ],
    )
    if not expression_correlation.empty:
        expression_correlation["FDR"] = bh_adjust(expression_correlation["pvalue"])
        expression_correlation["significant_at_configured_FDR"] = expression_correlation["FDR"].le(settings.fdr_cutoff)
        expression_correlation = expression_correlation.sort_values(["FDR", "pvalue"], na_position="last")

    gene_sets, gmt_audit = _load_gmt(settings.gene_sets_gmt)
    source_rows.append(gmt_audit)
    enrichment_rows: list[dict[str, object]] = []
    universe_symbols = set(hub_genes["gene_symbol"].astype(str).str.casefold()) - {""}
    for term, members in gene_sets.items():
        selected = members & universe_symbols
        for module, module_genes in hub_genes.groupby("module", observed=True):
            if str(module) not in eligible_modules:
                continue
            module_symbols = set(module_genes["gene_symbol"].astype(str).str.casefold()) - {""}
            overlap = module_symbols & selected
            enrichment_rows.append({
                "module": module, "term": term, "universe_genes": len(universe_symbols),
                "term_genes_in_universe": len(selected), "module_genes": len(module_symbols),
                "overlap_count": len(overlap), "overlap_genes": ";".join(sorted(overlap)),
                "pvalue": _hypergeom(len(universe_symbols), len(selected), len(module_symbols), len(overlap)),
            })
    gene_set_enrichment = pd.DataFrame(
        enrichment_rows,
        columns=[
            "module", "term", "universe_genes", "term_genes_in_universe", "module_genes",
            "overlap_count", "overlap_genes", "pvalue",
            "FDR", "significant_at_configured_FDR",
        ],
    )
    if not gene_set_enrichment.empty:
        gene_set_enrichment["FDR"] = bh_adjust(gene_set_enrichment["pvalue"])
        gene_set_enrichment["significant_at_configured_FDR"] = gene_set_enrichment["FDR"].le(settings.fdr_cutoff)
        gene_set_enrichment = gene_set_enrichment.sort_values(["FDR", "pvalue"], na_position="last")
    elif settings.gprofiler_enabled:
        try:
            module_symbols = {
                str(module): part["gene_symbol"].astype(str).tolist()
                for module, part in hub_genes.loc[hub_genes["module"].astype(str).isin(eligible_modules)].groupby("module", observed=True)
            }
            gene_set_enrichment, provider = gprofiler_enrichment(
                module_symbols, organism=settings.gprofiler_organism
            )
            if not gene_set_enrichment.empty:
                gene_set_enrichment["significant_at_configured_FDR"] = gene_set_enrichment["FDR"].le(settings.fdr_cutoff)
            source_rows.append({
                "source_type": "gene_sets_gprofiler", "path": provider["endpoint"],
                "status": "loaded", "rows": len(gene_set_enrichment),
                "data_version": str(provider.get("data_version", "")),
            })
        except Exception as exc:
            source_rows.append({
                "source_type": "gene_sets_gprofiler", "path": "g:Profiler g:GOSt",
                "status": "request_failed", "rows": 0, "message": str(exc),
            })

    isoform_switches = pd.DataFrame()
    isoform_enrichment = pd.DataFrame()
    if isoform_switches_override is not None and isoform_switches_override.exists():
        isoform_switches = pd.read_csv(isoform_switches_override, sep="\t", low_memory=False)
        isoform_enrichment = enrich_isoform_switches(isoform_switches, joint)
        source_rows.append({
            "source_type": "isoform_switch", "path": str(isoform_switches_override),
            "status": "cache_hit", "rows": len(isoform_switches),
        })
    elif quant_root is not None and gtf is not None and sample_metadata is not None and quant_root.exists() and gtf.exists():
        isoform_switches, isoform_enrichment = quantify_isoform_switches(
            quant_root, gtf, sample_metadata, joint,
            condition_column=condition_column,
            min_gene_tpm=settings.isoform_min_gene_tpm,
            min_abs_delta_if=settings.isoform_min_abs_delta_if,
        )
        source_rows.append({
            "source_type": "isoform_switch", "path": f"quant.sf under {quant_root}; {gtf}",
            "status": "generated_from_quant_gtf", "rows": len(isoform_switches),
        })
    else:
        source_rows.append({
            "source_type": "isoform_switch", "path": "", "status": "not_configured", "rows": 0,
        })

    regulator_rows: list[dict[str, object]] = []
    rbp_joint = joint.loc[joint["RBP_gene"]].copy()
    for _, row in rbp_joint.iterrows():
        module = row["module"]
        if str(module) not in eligible_modules:
            continue
        symbol = row["gene_symbol"]
        corr = expression_correlation.loc[
            expression_correlation["module"].eq(module)
            & expression_correlation["RBP"].astype(str).str.casefold().eq(str(symbol).casefold())
        ] if not expression_correlation.empty else pd.DataFrame()
        motif = motif_enrichment.loc[
            motif_enrichment["module"].eq(module) & motif_enrichment["RBP_or_motif"].astype(str).str.casefold().eq(str(symbol).casefold())
        ] if not motif_enrichment.empty else pd.DataFrame()
        corr_fdr = float(corr["FDR"].min()) if not corr.empty else np.nan
        motif_fdr = float(motif["FDR"].min()) if not motif.empty else np.nan
        evidence = [float(row["hub_gene_score"]), min(float(row["joint_support_count"]) / 4.0, 1.0)]
        if np.isfinite(corr_fdr):
            evidence.append(min(-np.log10(max(corr_fdr, 1e-300)) / 10.0, 1.0))
        if np.isfinite(motif_fdr):
            evidence.append(min(-np.log10(max(motif_fdr, 1e-300)) / 10.0, 1.0))
        regulator_rows.append({
            "module": module, "RBP": symbol, "hub_gene_score": row["hub_gene_score"],
            "joint_support_count": row["joint_support_count"], "expression_correlation_FDR": corr_fdr,
            "motif_enrichment_FDR": motif_fdr, "regulator_priority_score": float(np.mean(evidence)),
            "interpretation": "regulator prioritization only; not causal inference",
        })
    regulator_ranking = pd.DataFrame(
        regulator_rows,
        columns=[
            "module", "RBP", "hub_gene_score", "joint_support_count",
            "expression_correlation_FDR", "motif_enrichment_FDR",
            "regulator_priority_score", "interpretation",
            "regulator_rank",
        ],
    )
    if not regulator_ranking.empty:
        regulator_ranking["regulator_rank"] = regulator_ranking.groupby("module")["regulator_priority_score"].rank(method="min", ascending=False)
        regulator_ranking = regulator_ranking.sort_values(["module", "regulator_rank"])

    return ModuleInterpretationResult(
        hub_events=hub_events.sort_values(["module", "hub_event_rank"]).reset_index(drop=True),
        hub_genes=hub_genes.sort_values(["module", "hub_gene_rank"]).reset_index(drop=True),
        joint_evidence=joint.sort_values(["module", "joint_support_count", "hub_gene_score"], ascending=[True, False, False]).reset_index(drop=True),
        rbp_motif_enrichment=motif_enrichment.reset_index(drop=True),
        rbp_motif_hits=all_motif_hits.reset_index(drop=True),
        rbp_expression_correlation=expression_correlation.reset_index(drop=True),
        rbp_expression_matrix=(expression if expression is not None else pd.DataFrame()),
        gene_set_enrichment=gene_set_enrichment.reset_index(drop=True),
        isoform_switches=isoform_switches.reset_index(drop=True),
        module_isoform_switch_enrichment=isoform_enrichment.reset_index(drop=True),
        regulator_ranking=regulator_ranking.reset_index(drop=True),
        source_audit=pd.DataFrame(source_rows),
    )


def render_module_interpretation_plots(result: ModuleInterpretationResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    if not result.hub_events.empty:
        top = result.hub_events.loc[result.hub_events["hub_event_rank"].le(10)].copy()
        fig, ax = plt.subplots(figsize=(max(9, top["module"].nunique() * 1.5), max(5, len(top) * 0.18)))
        sns.scatterplot(data=top, x="module", y="hub_event_rank", size="hub_event_score", hue="module_membership_kME", palette="vlag", ax=ax, legend="brief")
        ax.invert_yaxis()
        ax.set_title("Top hub events per CoSplice module")
        fig.tight_layout()
        path = output_dir / "module_hub_events.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    if not result.rbp_expression_correlation.empty:
        significant = result.rbp_expression_correlation.nsmallest(100, "FDR")
        matrix = significant.pivot_table(index="RBP", columns="module", values="correlation", aggfunc="first")
        if not matrix.empty:
            fig, ax = plt.subplots(figsize=(max(8, matrix.shape[1] * 0.8), max(5, matrix.shape[0] * 0.3)))
            sns.heatmap(matrix, cmap="vlag", center=0, ax=ax, cbar_kws={"label": "RBP expression vs module score (Pearson r)"})
            ax.set_title("RBP expression support for CoSplice modules")
            fig.tight_layout()
            path = output_dir / "rbp_module_expression_correlation.png"
            fig.savefig(path, dpi=220, bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
    if not result.rbp_motif_enrichment.empty:
        frame = result.rbp_motif_enrichment.nsmallest(40, "FDR").copy()
        frame["label"] = frame["RBP_or_motif"].astype(str) + " | " + frame["module"].astype(str)
        frame["-log10 FDR"] = -np.log10(pd.to_numeric(frame["FDR"], errors="coerce").clip(lower=1e-300))
        fig, ax = plt.subplots(figsize=(9, max(5, len(frame) * 0.25)))
        sns.barplot(data=frame, x="-log10 FDR", y="label", color="#0f766e", ax=ax)
        ax.set_title("RBP motif support around CoSplice events")
        fig.tight_layout()
        path = output_dir / "rbp_motif_module_enrichment.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    if not result.module_isoform_switch_enrichment.empty:
        frame = result.module_isoform_switch_enrichment.copy()
        if "switch_scope" in frame:
            frame = frame.loc[frame["switch_scope"].eq("DTU_supported_switch")]
        if not frame.empty:
            frame["-log10 FDR"] = -np.log10(pd.to_numeric(frame["FDR"], errors="coerce").clip(lower=1e-300))
            fig, ax = plt.subplots(figsize=(max(8, len(frame) * 0.7), 5))
            sns.barplot(data=frame, x="module", y="-log10 FDR", color="#d97706", ax=ax)
            ax.set_title("DTU-supported dominant-isoform switch enrichment by CoSplice module")
            fig.tight_layout()
            path = output_dir / "module_isoform_switch_enrichment.png"
            fig.savefig(path, dpi=220, bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
    return paths
