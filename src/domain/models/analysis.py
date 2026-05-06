from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass
class ThresholdConfig:
    splicing_fdr: float = 0.05
    splicing_dpsi: float = 0.10
    program_delta_dpsi: float = 0.10
    deg_padj: float = 0.05
    deg_log2fc: float = 1.0
    dexseq_qvalue: float = 0.05
    dtu_qvalue: float = 0.05


@dataclass
class AnalysisRequest:
    comparison_ids: list[str] = field(default_factory=list)
    program_comparison_ids: list[str] = field(default_factory=list)
    analysis_modules: list[str] = field(default_factory=list)
    rmats_mode: str = "JC"
    reverse_directions: dict[str, bool] = field(default_factory=dict)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)


@dataclass
class AnalysisResults:
    splicing_landscape: pd.DataFrame = field(default_factory=pd.DataFrame)
    program_events: pd.DataFrame = field(default_factory=pd.DataFrame)
    program_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    mechanism_support: pd.DataFrame = field(default_factory=pd.DataFrame)
    tx_splicing_gene_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    tx_splicing_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    candidate_gene_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    cross_comparison_candidate_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    ranked_candidates: pd.DataFrame = field(default_factory=pd.DataFrame)
    shortlist_dp: pd.DataFrame = field(default_factory=pd.DataFrame)
    shortlist_ko: pd.DataFrame = field(default_factory=pd.DataFrame)
    cards_shortlist: pd.DataFrame = field(default_factory=pd.DataFrame)
    cards_expression_support: pd.DataFrame = field(default_factory=pd.DataFrame)
    sashimi_manifest: pd.DataFrame = field(default_factory=pd.DataFrame)
    isoform_manifest: pd.DataFrame = field(default_factory=pd.DataFrame)
    jutils_manifest: pd.DataFrame = field(default_factory=pd.DataFrame)


@dataclass
class AnalysisRunState:
    request: AnalysisRequest | None = None
    status: str = "idle"
    message: str = "No analysis has been run yet."
    results: AnalysisResults = field(default_factory=AnalysisResults)
