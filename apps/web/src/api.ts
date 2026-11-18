import type {
  Artifact,
  Job,
  ProjectPayload,
  Selection,
  TablePage,
  WorkflowDefinition,
  WorkflowCache,
  Workspace,
} from "./types";
import { displayModuleName } from "./branding";

const API_ROOT = import.meta.env.VITE_API_ROOT ?? "/api";

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${API_ROOT}${path}`, {
    headers: { "Content-Type": "application/json", ...init?.headers },
    ...init,
  });
  if (!response.ok) {
    const payload = await response.json().catch(() => ({ detail: response.statusText }));
    throw new Error(payload.detail ?? `Request failed (${response.status})`);
  }
  return response.json() as Promise<T>;
}

export const api = {
  browseLocalFiles: (path = "") => request<import("./components/PathPicker").LocalDirectory>("/local-files/browse", {method:"POST",body:JSON.stringify({path})}),
  prepareTargetPipeline: (project:string,payload:unknown) => request<{task_id:string;status:string;command:string;download_url:string}>(`/projects/${encodeURIComponent(project)}/targets/nextflow`, {method:"POST",body:JSON.stringify(payload)}),
  importReport: (payload: {report_path:string;output_root:string}) => request<ProjectPayload>("/projects/import-report", {method:"POST",body:JSON.stringify(payload)}),
  browseReport: (project:string, artifact:string, payload:unknown) => request<TablePage>(`/projects/${encodeURIComponent(project)}/report-tables/${encodeURIComponent(artifact)}`, {method:"POST",body:JSON.stringify(payload)}),
  plotReport: (project:string, artifact:string, payload:unknown) => request<import("./components/PopulationPlot").Population>(`/projects/${encodeURIComponent(project)}/report-tables/${encodeURIComponent(artifact)}/plot`, {method:"POST",body:JSON.stringify(payload)}),
  selectReport: (project:string, artifact:string, payload:unknown) => request<{rows:Record<string,unknown>[];total:number}>(`/projects/${encodeURIComponent(project)}/report-tables/${encodeURIComponent(artifact)}/selection`, {method:"POST",body:JSON.stringify(payload)}),
  enrichReport: (project:string, artifact:string, payload:unknown) => request<import("./components/LocalEnrichment").EnrichmentResult>(`/projects/${encodeURIComponent(project)}/report-tables/${encodeURIComponent(artifact)}/enrichment`, {method:"POST",body:JSON.stringify(payload)}),
  compactReport: (project:string) => request<{archive_url:string;directory:string}>(`/projects/${encodeURIComponent(project)}/report-compact`, {method:"POST"}),
  exportReportRows: async (project:string, artifact:string, payload:unknown) => {
    const response=await fetch(`${API_ROOT}/projects/${encodeURIComponent(project)}/report-tables/${encodeURIComponent(artifact)}/export`, {method:"POST",headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});
    if(!response.ok) throw new Error((await response.json()).detail);
    const url=URL.createObjectURL(await response.blob());const a=document.createElement('a');a.href=url;a.download='filtered-results.tsv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  },
  motifTraceModes: () => request<MotifTraceCatalog>("/motiftrace/modes"),
  configureGeneFilter: (projectId: string, settings: import("./components/GeneFilterPanel").GeneFilterSettings) =>
    request<ProjectPayload>(`/projects/${encodeURIComponent(projectId)}/gene-filter`, {method:"POST",body:JSON.stringify(settings)}),
  moduleDashboard: (projectId: string, artifactId: string, module: string, condition = "all", selected_event?: string, hide_same_gene = false) =>
    request<ModuleDashboard>(`/projects/${encodeURIComponent(projectId)}/artifacts/${encodeURIComponent(artifactId)}/module-dashboard?${new URLSearchParams({module, condition, hide_same_gene: String(hide_same_gene), ...(selected_event ? {selected_event} : {})})}`),
  lineage: (project: string, artifact: string) => request<{scope:string; nodes:Array<{artifact_id:string;name:string;workflow:string;entity_type:string;rows:number|null;state:string;parents?:string[]}>}>(`/projects/${encodeURIComponent(project)}/artifacts/${encodeURIComponent(artifact)}/lineage`),
  editableDesign: (id: string) => request<{revision: string; samples: import("./components/DesignEditor").DesignSample[]}>(`/projects/${encodeURIComponent(id)}/design/edit`),
  editDesign: (id: string, edit: unknown, apply: boolean) => request<{applied: boolean; revision?: string; review: DesignReview; project?: ProjectPayload}>(`/projects/${encodeURIComponent(id)}/design/edit?apply=${apply}`, {method: "POST", body: JSON.stringify(edit)}),
  moduleEvidence: (projectId: string, artifactId: string, module: string) =>
    request<ModuleEvidence>(`/projects/${encodeURIComponent(projectId)}/artifacts/${encodeURIComponent(artifactId)}/module-evidence?${new URLSearchParams({module})}`),
  workflowPlan: (projectId: string, workflowId: string, parameters: Record<string, unknown>) =>
    request<RunPlan>(`/projects/${encodeURIComponent(projectId)}/workflows/${encodeURIComponent(workflowId)}/plan`, {
      method: "POST", body: JSON.stringify({ parameters, force: false }),
    }),
  projectDesign: (projectId: string) => request<DesignReview>(`/projects/${encodeURIComponent(projectId)}/design`),
  workspaces: () => request<Workspace[]>("/workspaces"),
  workflows: () => request<WorkflowDefinition[]>("/workflows").then((workflows) => {
    const primary = ["programs.cosplice", "programs.psi_network", "programs.hub_events", "programs.djexpress"];
    const priority = (id: string) => primary.includes(id) ? primary.indexOf(id) : primary.length;
    return workflows.sort((a,b)=>priority(a.workflow_id)-priority(b.workflow_id)).map((workflow) => ({
    ...workflow,
    label: displayModuleName(workflow.label),
    question: displayModuleName(workflow.question),
    method: displayModuleName(workflow.method),
    output_summary: displayModuleName(workflow.output_summary),
    limitations: displayModuleName(workflow.limitations),
    parameters: workflow.parameters.map((parameter) => ({...parameter, label: displayModuleName(parameter.label), help: displayModuleName(parameter.help)})),
  })); }),
  singleCellConditions: (projectId: string) => request<string[]>(`/projects/${encodeURIComponent(projectId)}/single-cell/conditions`),
  openProject: (manifestPath: string) =>
    request<ProjectPayload>("/projects/open", {
      method: "POST",
      body: JSON.stringify({ manifest_path: manifestPath }),
    }),
  migrateProject: (payload: { legacy_config_path: string; output_root: string; manifest_path?: string }) =>
    request<ProjectPayload>("/projects/migrate", {
      method: "POST",
      body: JSON.stringify(payload),
    }),
  artifacts: (projectId: string, analysisId?: string) => {
    const query = analysisId ? `?analysis_id=${encodeURIComponent(analysisId)}` : "";
    return request<Artifact[]>(`/projects/${encodeURIComponent(projectId)}/artifacts${query}`);
  },
  configureResources: (payload: { manifest_path: string; resource_bundle_path: string; output_root: string }) =>
    request<ProjectPayload>("/projects/configure-resources", { method: "POST", body: JSON.stringify(payload) }),
  workflowCache: (
    projectId: string,
    workflowId: string,
    parameters: Record<string, unknown>,
  ) =>
    request<WorkflowCache>(
      `/projects/${encodeURIComponent(projectId)}/workflows/${encodeURIComponent(workflowId)}/cache`,
      {
        method: "POST",
        body: JSON.stringify({ parameters, force: false }),
      },
    ),
  artifactRows: (
    projectId: string,
    artifactId: string,
    options: { offset: number; limit: number; search: string; sortBy?: string; descending?: boolean },
  ) => {
    const query = new URLSearchParams({
      offset: String(options.offset),
      limit: String(options.limit),
      search: options.search,
    });
    if (options.sortBy) query.set("sort_by", options.sortBy);
    if (options.descending) query.set("descending", "true");
    return request<TablePage>(
      `/projects/${encodeURIComponent(projectId)}/artifacts/${encodeURIComponent(artifactId)}/rows?${query}`,
    );
  },
  artifactContent: (projectId: string, artifactId: string) =>
    request<Record<string, unknown>>(
      `/projects/${encodeURIComponent(projectId)}/artifacts/${encodeURIComponent(artifactId)}/content`,
    ),
  artifactFileUrl: (projectId: string, artifactId: string, download = false) =>
    `${API_ROOT}/projects/${encodeURIComponent(projectId)}/artifacts/${encodeURIComponent(artifactId)}/file${
      download ? "?download=true" : ""
    }`,
  artifactExportUrl: (projectId: string, artifactId: string, format: string) =>
    `${API_ROOT}/projects/${encodeURIComponent(projectId)}/artifacts/${encodeURIComponent(artifactId)}/export/${format}`,
  selections: (projectId: string) =>
    request<Selection[]>(`/projects/${encodeURIComponent(projectId)}/selections`),
  createSelection: (
    projectId: string,
    payload: {
      name: string;
      entity_type: string;
      source_artifact_id: string;
      entity_ids: string[];
      filters: Record<string, unknown>;
    },
  ) =>
    request<Selection>(`/projects/${encodeURIComponent(projectId)}/selections`, {
      method: "POST",
      body: JSON.stringify(payload),
    }),
  runWorkflow: (projectId: string, workflowId: string, parameters: Record<string, unknown>, force = false) =>
    request<Job>(
      `/projects/${encodeURIComponent(projectId)}/workflows/${encodeURIComponent(workflowId)}/run`,
      {
        method: "POST",
        body: JSON.stringify({ parameters, force }),
      },
    ),
  openOutput: (projectId: string) =>
    request<{ opened: boolean; path: string }>(`/projects/${encodeURIComponent(projectId)}/open-output`, {
      method: "POST",
    }),
  jobs: (projectId: string) => request<Job[]>(`/jobs?project_id=${encodeURIComponent(projectId)}`),
  retryJob: (jobId: string) => request<Job>(`/jobs/${encodeURIComponent(jobId)}/retry`, {method: "POST"}),
  cancelJob: (jobId: string) => request<Job>(`/jobs/${encodeURIComponent(jobId)}/cancel`, { method: "POST" }),
  exportTargets: (projectId: string, payload: import("./types").TargetBatch) =>
    request<import("./types").TargetExport>(`/projects/${encodeURIComponent(projectId)}/targets/export`, {method: "POST", body: JSON.stringify(payload)}),
  planTargets: (projectId: string, payload: import("./types").TargetBatch) =>
    request<import("./types").TargetPlan>(`/projects/${encodeURIComponent(projectId)}/targets/plan`, {method: "POST", body: JSON.stringify(payload)}),
  runTargets: (projectId: string, payload: import("./types").TargetBatch) =>
    request<import("./types").TargetRun>(`/projects/${encodeURIComponent(projectId)}/targets/run`, {method: "POST", body: JSON.stringify(payload)}),
  targetRun: (projectId: string, runId: string) =>
    request<import("./types").TargetRun>(`/projects/${encodeURIComponent(projectId)}/targets/runs/${encodeURIComponent(runId)}`),
  importTargetResults: (projectId: string, manifestPath: string) =>
    request<import("./types").TargetRun>(`/projects/${encodeURIComponent(projectId)}/targets/results/import`, {method: "POST", body: JSON.stringify({manifest_path: manifestPath})}),
  targetFileUrl: (projectId: string, runId: string, path: string) =>
    `${API_ROOT}/projects/${encodeURIComponent(projectId)}/targets/runs/${encodeURIComponent(runId)}/files/${path.split("/").map(encodeURIComponent).join("/")}`,
  splicetraceWorkbench: (projectId: string, comparison = "") =>
    request<SpliceTracePayload>(`/projects/${encodeURIComponent(projectId)}/splicetrace-workbench?${new URLSearchParams(comparison ? {comparison} : {})}`),
  quadrantExplorer: (projectId: string, comparison = "") =>
    request<QuadrantPayload>(`/projects/${encodeURIComponent(projectId)}/quadrant-explorer?${new URLSearchParams(comparison ? {comparison} : {})}`),
  contrastRescue: (projectId: string, payload: {contrast_a: string; contrast_b: string; min_abs_dpsi: number; max_fdr: number; relation: "reversed" | "concordant" | "specific_to_a"; rs_only: boolean}) =>
    request<ContrastRescueResult>(`/projects/${encodeURIComponent(projectId)}/contrast-rescue`, {method: "POST", body: JSON.stringify(payload)}),
  singleCellManifold: (projectId: string, feature = "") =>
    request<SingleCellManifoldPayload>(`/projects/${encodeURIComponent(projectId)}/single-cell-manifold?${new URLSearchParams(feature ? {feature} : {})}`),
  genecardPreview: (projectId: string, payload: {gene_query?: string; event_id?: string; comparison?: string; custom_sequence?: string; custom_motif?: string}) =>
    request<GeneCardPayload>(`/projects/${encodeURIComponent(projectId)}/genecard-preview`, {method: "POST", body: JSON.stringify(payload)}),
  computePsiLandscape: (projectId: string, payload: {comparison?: string; mode?: "unsupervised" | "supervised"; top_n?: number; prior_kappa?: number; z_clip?: number; metric?: "euclidean" | "cosine" | "correlation"; event_types?: string[]; umap_neighbors?: number; umap_min_dist?: number}) =>
    request<LivePsiLandscapePayload>(`/projects/${encodeURIComponent(projectId)}/compute/psi-landscape`, {method: "POST", body: JSON.stringify(payload)}),
  computeLiveOra: (projectId: string, payload: {comparison?: string; preset?: string; custom_genes?: string[]; fdr_cutoff?: number; effect_cutoff?: number; min_term_size?: number; max_term_size?: number}) =>
    request<LiveOraResultPayload>(`/projects/${encodeURIComponent(projectId)}/compute/live-ora`, {method: "POST", body: JSON.stringify(payload)}),
  computeMotifScan: (projectId: string, payload: {comparison?: string; query_motif_or_rbp?: string; custom_sequence?: string}) =>
    request<LiveMotifScanResult>(`/projects/${encodeURIComponent(projectId)}/compute/motif-scan`, {method: "POST", body: JSON.stringify(payload)}),
  computeTargetRanker: (projectId: string, payload: {comparison?: string; w_dpsi?: number; w_isoform?: number; w_splicenet?: number; w_motif?: number; w_proteome?: number; min_abs_dpsi?: number; max_fdr?: number; deg_independent_only?: boolean; require_domain_or_nmd?: boolean}) =>
    request<LiveTargetRankerResult>(`/projects/${encodeURIComponent(projectId)}/compute/target-ranker`, {method: "POST", body: JSON.stringify(payload)}),
};

export type SpliceTracePayload = {
  available: boolean;
  comparison_id: string;
  rbps: Array<{
    rbp: string;
    best_region: string;
    stat_z: number | null;
    nominal_p: number | null;
    fwer_p: number | null;
    host_deg_log2fc: number | null;
    host_deg_fdr: number | null;
    pwm?: {motif_id: string; matrix: number[][]} | null;
  }>;
  regions: Record<string, Array<{region: string; stat_z: number | null; nominal_p: number | null}>>;
  targets: Array<{
    event_id: string;
    gene_id: string;
    gene_symbol: string;
    event_type: string;
    dpsi: number | null;
    frame_status: string;
    nmd_switch_class: string;
    delta_aa: number | null;
    disrupted_protein_features: string;
    cis_rbps: string;
  }>;
  bridge: Array<Record<string, unknown>>;
};

export type QuadrantPayload = {
  comparison_id: string;
  das_artifact_id?: string;
  deg_artifact_id?: string;
  total_co_measured_genes?: number;
  rows: Array<{
    gene_id: string;
    gene_symbol: string;
    event_id: string;
    das_row_id: string;
    log2fc_val: number;
    deg_fdr_val: number;
    dpsi_val: number;
    das_fdr_val: number;
    rs: boolean;
    nmd_switch_class?: string | null;
    disrupted_protein_features?: string | null;
  }>;
};

export type ContrastRescueResult = {
  contrast_a: string;
  contrast_b: string;
  relation: string;
  source_artifact_id?: string;
  total: number;
  rows: Array<{
    event_id: string;
    gene_id?: string;
    gene_symbol: string;
    event_type?: string;
    row_id_a?: string;
    dpsi_a: number;
    fdr_a: number;
    dpsi_b: number | null;
    fdr_b: number | null;
    nmd_switch_class?: string | null;
    disrupted_protein_features?: string | null;
  }>;
};

export type SingleCellManifoldPayload = {
  available: boolean;
  feature_label: string;
  module_columns: string[];
  cell_types: string[];
  donors: string[];
  psi_events: string[];
  expr_genes: string[];
  cells: Array<{
    cell_id: string;
    donor_id: string;
    condition: string;
    activation_arm: string;
    cell_type: string;
    expr_PC1: number;
    expr_PC2: number;
    psi_PC1: number;
    psi_PC2: number;
    feature_value?: number | null;
    [key: string]: unknown;
  }>;
};

export type GeneCardPayload = {
  found: boolean;
  reason?: string;
  comparison_id?: string;
  gene_id?: string;
  gene_symbol?: string;
  event_id?: string;
  event_type?: string;
  dpsi?: number | null;
  fdr?: number | null;
  das_artifact_id?: string;
  das_row_id?: string;
  host_events?: Array<{event_id: string; event_type: string; dpsi: number | null; fdr: number | null; row_id: string}>;
  track1_quantification?: {
    deg: {log2fc: number | null; fdr: number | null; baseMean: number | null};
    expression_samples: Array<{sample_id: string; condition: string; arm: string; value: number; unit: string}>;
    psi_samples: Array<{sample_id: string; condition: string; psi: number; ijc: number | null; sjc: number | null}>;
    dtu_isoforms: Array<{transcript_id: string; effect: number | null; fdr: number | null}>;
  };
  track2_sashimi?: {
    chrom: string;
    strand: string;
    event_type: string;
    exons: Array<{label: string; role: string; start: number; end: number}>;
    junctions: Array<{label: string; role: string; from_idx: number; to_idx: number}>;
    arms: Array<{label: string; mean_psi: number | null; mean_ijc: number | null; mean_sjc: number | null; n_samples: number}>;
    bam_bound_samples: number;
  };
  track3_motifs?: {
    status: "present" | "scanned_zero" | "unavailable";
    mode?: string;
    hits?: Array<{rbp: string; region: string; delta_score: number | null; consensus?: string; slot: number}>;
    note?: string;
  };
  track4_proteome?: {
    status: "coding_architecture" | "noncoding" | "scanned_unmapped" | "unavailable";
    computed_mode?: string;
    transcript_id?: string;
    uniprot_acc?: string;
    coding_region_location?: string;
    frame_status?: string;
    nmd_switch_class?: string;
    aa_len_inclusion?: number | null;
    aa_len_exclusion?: number | null;
    delta_aa?: number | null;
    ptc_to_last_ejc_nt?: number | null;
    disrupted_protein_features?: string;
    summary?: string;
    note?: string;
    architecture?: Record<string, unknown>;
  };
};

export type LivePsiLandscapePayload = {
  available: boolean;
  reason?: string;
  run_id?: string;
  assay_mode?: "bulk" | "single-cell";
  selection_mode: string;
  metric: string;
  prior_kappa: number;
  z_clip: number;
  umap_method: string;
  total_candidate_events: number;
  selected_events: number;
  manifold_events?: number;
  heatmap_events?: number;
  permanova_r2: number;
  permanova_p: number;
  pca_variance_pct: number[];
  dual_manifold_metrics?: {
    mean_divergence: number;
    median_donor_lisi: number;
    cell_type_silhouette_psi: number;
    cell_type_silhouette_expr: number;
  };
  samples: Array<{
    sample_id: string;
    condition: string;
    donor_id: string;
    cell_type: string;
    observed_event_fraction: number;
    PC1: number;
    PC2: number;
    PC3: number;
    psi_UMAP1: number;
    psi_UMAP2: number;
    expr_UMAP1?: number;
    expr_UMAP2?: number;
    modality_divergence?: number;
    donor_lisi?: number;
  }>;
  ordered_samples: string[];
  ordered_conditions: string[];
  heatmap_rows: Array<{
    event_id: string;
    gene_symbol: string;
    event_type: string;
    excess_biovar: number;
    z_values: number[];
    psi_values: number[];
  }>;
  files?: Record<string, string>;
};

export type LiveOraRow = {
  ontology_category: string;
  category_code: string;
  term: string;
  overlap_genes: number;
  term_size_in_bg: number;
  selected_genes: number;
  background_genes: number;
  gene_ratio: number;
  fold_enrichment: number;
  pvalue: number;
  fdr: number;
  hit_genes: string;
};

export type LiveOraResultPayload = {
  available: boolean;
  reason?: string;
  run_id?: string;
  species?: string;
  comparison_id?: string;
  preset_label?: string;
  foreground_count?: number;
  background_count?: number;
  significant_fdr05_count?: number;
  facets?: Record<string, LiveOraRow[]>;
  top_rows?: LiveOraRow[];
  files?: Record<string, string>;
};

export type LiveMotifScanResult = {
  query: string;
  comparison_id: string;
  motif_id: string;
  consensus: string;
  pwm_matrix: number[][] | null;
  catalog_size: number;
  available_rbps: string[];
  homologous_rbps: Array<{rbp: string; motif_id: string; consensus: string}>;
  inclusion_hits: number;
  skipping_hits: number;
  direction_binomial_p: number;
  matched_events: Array<{
    event_id: string;
    gene_symbol: string;
    event_type: string;
    dpsi: number | null;
    region: string;
    nmd_switch_class: string;
    disrupted_protein_features: string;
  }>;
  sequence_scan_hits: Array<{rbp: string; region: string; delta_score: number; consensus: string; count: number}>;
};

export type LiveTargetRankerResult = {
  available: boolean;
  comparison_id: string;
  total_qualified: number;
  rows: Array<{
    comparison_id: string;
    event_id: string;
    gene_id: string;
    gene_symbol: string;
    event_type: string;
    dpsi: number;
    fdr: number;
    splicenet_module?: string;
    cis_rbp_motifs?: string;
    nmd_switch_class?: string;
    disrupted_uniprot_domains?: string;
    composite_score: number;
  }>;
};

export type MotifTraceMode = {
  mode_id: string; label: string; workflow_id: string;
  fixed_parameters: Record<string, unknown>; description: string;
  operation: "new_analysis" | "completed_result_import" | "existing_result_integration"; requires_selection: boolean;
};
export type MotifTraceCatalog = {
  package_version: string; scientific_release: string; modes: MotifTraceMode[];
};

export type RunPlan = { ready: boolean; errors: string[]; warnings: string[]; scope: string; input_bytes: number;
  memory_estimate: string; parent_artifacts: {artifact_id: string; name: string; workflow: string; run_fingerprint: string}[] };
export type DesignReview = { factor_policy: string; batch_design: {reason: string; estimable: boolean | null};
  contrasts: {contrast_id: string; label: string; numerator: string; denominator: string; numerator_biological_units: number;
    denominator_biological_units: number; pairing: string; effect_definition: string}[];
  comparison_pairs: {contrast_A: string; contrast_B: string; label_A: string; label_B: string; matched_denominator: boolean;
    shared_samples: string[]; shared_biological_units: string[]; interpretation: string}[] };
export type ModuleDashboard = {
  module: string; condition: string; conditions: string[]; status: string; method: string;
  biological_samples: number; minimum_biological_samples: number; total_events: number; total_edges: number; same_gene_edges: number;
  nodes: {event_key: string; gene: string; connectivity: number; x: number; y: number}[];
  edges: {event_a: string; event_b: string; weight: number; correlation: number | null; pair_n: number | null; same_gene: boolean}[];
  hub_events: HubEvent[];
  hub_status?: string; hub_method?: string; hub_analysis_run?: string; hub_comparison?: HubEvent[];
  hub_label_status?: string;
  network_label_policy?: {policy_id?: string; apply_threshold_labels?: boolean; hard_filter?: boolean}|null;
  hub_label_policy?: {policy_id?: string; apply_threshold_labels?: boolean; hard_filter?: boolean}|null;
  event_options: {event_key: string; gene: string}[];
  heatmap: {events: string[]; samples: string[]; conditions: string[]; z: (number|null)[][]; scale: string};
  selected_event: string; event_profile: {sample_id: string; biological_sample_id: string; condition: string; psi: number|null}[];
  scores: {sample_id: string; condition: string; value: number|null}[]; score_definition: string;
  display_rule: string; interpretation: string; discovery_run: string; connectivity_definition: string;
  enrichment: {term: string; overlap_count: number; overlap_genes: string; pvalue: number|null; FDR: number|null}[];
  enrichment_total: number; enrichment_scope: string;
  sources: {artifact_id: string; name: string; run_fingerprint: string}[];
};

export type HubEvent = {event_key: string; gene: string; connectivity: number;
  kWithin?: number|null; kWithin_no_same_gene?: number|null; kME?: number|null;
  rank_all_edges?: number|null; rank_no_same_gene?: number|null; rank_change_no_same_gene?: number|null;
  bootstrap_frequency_all_edges?: number|null; bootstrap_frequency_no_same_gene?: number|null;
  hub_candidate?: boolean|null; candidate_all_edges?: boolean|null; candidate_no_same_gene?: boolean|null; selection_reason?: string;
  preset_selection_frequency_all_edges?: number|null; preset_selection_frequency_no_same_gene?: number|null;
  label_status?: string; hub_label_status?: string; estimation_status?: string; estimation_reason?: string};

export type ModuleEvidence = {module: string; discovery_run: string; scope: string; tables: {
  artifact_id: string; name: string; entity_type: string; row_count: number; preview_limit: number; rows: Record<string,unknown>[]; run_fingerprint: string;
}[]};
