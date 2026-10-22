// Display names never replace durable workflow IDs, artifact names or cache keys.
export const SPLICE_MODULE_NAME = "Katana-Splice-Network";
export const MOTIFTRACE_MODULE_NAME = "Katana MotifTrace";

export function displayModuleName(value: string): string {
  return value.replace(/(?<![\w./\\-])(?:Katana-Splice(?:-Network)?|katana-splice-network|(?:Katana-)?CoSplice|cosplice)(?![\w./\\-])/g, SPLICE_MODULE_NAME);
}

export function displayArtifactName(value: string): string {
  return ({
    consensus_modules__modules: "Module sizes and stability",
    consensus_modules__assignments: "Event membership and confidence",
    consensus_modules__bootstrap: "50 resampling runs: coverage and confidence",
    consensus_modules__scores: "Fixed-module activity across samples",
    network_observed__scores: "Discovery-module activity across samples",
    network_observed__edges: "Splice-event association edges",
    stable_qc__matrix: "Measured event PSI matrix",
    stable_qc__metrics: "Event measurement and support checks",
    matrix__events: "Event coordinates and gene names",
    matrix__imputation_means: "Discovery-only imputation audit",
    module_functions__enrichment: "Module functions and contributing genes",
    module_functions__background: "Measured background genes",
  }[value] ?? displayModuleName(value));
}

export function displayParameterOption(key: string, value: string): string {
  if (key === "event_set") return ({stable_significant: "RS-DAS: significant + read support (default)", significant: "Ordinary DAS: significant events"}[value] ?? value);
  if (key === "pairing") return ({auto:"Detect from sample identities", paired:"Paired samples", unpaired:"Independent samples"}[value] ?? value);
  return ({
    direct_experiment: "Experimental motifs only",
    all_evidence_separate: "Experimental, inferred and predicted (separate)",
    boundary_association: "Upstream motif / CLIP regulatory clues",
    se_sequence_impact: "Legacy SE fixed-site threshold analysis",
    completed_sequence_maps: "Katana MotifTrace: completed SE RNA and motifs",
    event_consequences: "Katana MotifTrace: paired RNA / five event types / region statistics",
    completed_parallel_models: "MotifTrace: completed MAPP / usage results (separate P values)",
    preview_only: "User library: sequence positions only (no new P/FDR)",
    fixed_enrichment: "User library: explicit enrichment in a separate FDR family",
    selected_comparison: "Selected comparison only (original FDR retained)",
    all_frozen_comparisons: "All comparisons in the frozen result (original FDR retained)",
    if_configured: "Use available CLIP data",
    required: "Require CLIP data",
    diagnostic_only: "Report uncertainty separately",
    raw_reliability_weighted: "Weight observed PSI (experimental)",
  }[value] ?? value);
}
