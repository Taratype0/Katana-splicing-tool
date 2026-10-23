export type Workspace = {
  workspace_id: string;
  label: string;
  purpose: string;
  icon: string;
  enabled: boolean;
  phase: number;
};

export type WorkflowParameter = {
  key: string;
  label: string;
  kind: "number" | "select" | "boolean" | "text" | "selection" | "condition" | "artifact";
  default: unknown;
  help: string;
  options: string[];
  minimum?: number;
  maximum?: number;
  step?: number;
};

export type WorkflowDefinition = {
  workflow_id: string;
  workspace_id: string;
  label: string;
  question: string;
  method: string;
  output_summary: string;
  analysis_unit: string;
  expected_artifacts: number;
  parameters: WorkflowParameter[];
  required_asset_kinds: string[];
  optional_asset_kinds: string[];
  supports_contrast_selection: boolean;
  enabled: boolean;
  limitations: string;
  presentation_role?: "primary" | "advanced" | "auxiliary";
  das_selection?: {
    mode: "switchable" | "rs_required" | "measurement_population" | "not_applicable";
    label: string;
    description: string;
    default_event_set: string | null;
    policy_version?: string;
  };
};

export type FileFingerprint = {
  path: string;
  exists: boolean;
  size?: number;
  digest?: string;
  hash_method: string;
};

export type ManifestSample = {
  sample_id: string;
  biological_sample_id: string;
  condition: string;
  batch?: string;
  attributes: Record<string, string>;
  assets: Record<string, string>;
};

export type ManifestContrast = {
  contrast_id: string;
  label: string;
  numerator: { label: string; sample_ids: string[] };
  denominator: { label: string; sample_ids: string[] };
};

export type ProjectPayload = {
  snapshot: {
    project_id: string;
    manifest_path: string;
    name: string;
    input_root: string;
    output_root: string;
    genome_build: string;
    manifest_fingerprint: string;
    input_fingerprints: FileFingerprint[];
    sample_count: number;
    contrast_count: number;
    warnings: string[];
  };
  manifest: {
    name: string;
    contrasts: ManifestContrast[];
    samples: ManifestSample[];
    metadata: Record<string, unknown>;
  };
  artifact_count: number;
  selection_count: number;
};

export type Artifact = {
  artifact_id: string;
  project_id: string;
  analysis_id: string;
  name: string;
  kind: "table" | "plot_spec" | "selection" | "report" | "file";
  entity_type: string;
  status: string;
  fingerprint: string;
  run_fingerprint: string;
  path: string;
  row_count?: number | null;
  parent_artifact_ids: string[];
  contrast_ids: string[];
  parameters: Record<string, unknown>;
  metadata: Record<string, unknown>;
  warnings: string[];
  created_at: string;
};

export type WorkflowCache = {
  run_fingerprint: string;
  cache_hit: boolean;
  artifacts: Artifact[];
};

export type TablePage = {
  columns: string[];
  rows: Array<Record<string, unknown>>;
  offset: number;
  limit: number;
  total: number;
};

export type JobEvent = {
  sequence: number;
  timestamp: string;
  progress: number;
  message: string;
};

export type Job = {
  job_id: string;
  project_id: string;
  analysis_id: string;
  status: "queued" | "running" | "complete" | "failed" | "cancelled" | "cache_hit";
  run_fingerprint: string;
  progress: number;
  message: string;
  cache_status: "hit" | "rebuild";
  artifact_ids: string[];
  events: JobEvent[];
  error?: string;
  created_at: string;
};

export type Selection = {
  compatible_workflows?: string[];
  selection_id: string;
  project_id: string;
  name: string;
  entity_type: string;
  source_artifact_id: string;
  entity_ids: string[];
  created_at: string;
};

export type TargetKind = "gene" | "event";
export type TargetSelection = {
  source_artifact_id: string;
  row_ids: string[];
  target_kind: TargetKind;
  filters: Record<string, unknown>;
};
export type TargetBasketEntry = TargetSelection & {
  basket_id: string;
  source_name: string;
  labels: string[];
};
export type TargetResource = {resource_id: string; kind: string; sample_id?: string; [key: string]: unknown};
export type TargetExport = {
  request_id: string; request_path: string; archive_path: string; bundle_directory: string;
  target_count: number; targets: Array<Record<string, unknown>>; resources: TargetResource[];
  download_url: string;
};
export type TargetBatch = {
  selections?: TargetSelection[]; request_id?: string; modules?: string[];
  contrast_ids?: string[]; resources?: Record<string, string>; workers?: number;
};
export type TargetPlanMessage = string | {message?: string; target_id?: string; resource_id?: string; [key: string]: unknown};
export type TargetPlan = {ready: boolean; errors: TargetPlanMessage[]; warnings: TargetPlanMessage[]; tasks: unknown[]; request_id: string; [key: string]: unknown};
export type TargetRun = {
  run_id: string; request_id: string; status: string; execution?: string; error?: string;
  progress?: Record<string, unknown>; report_url?: string; archive_url?: string; manifest_url?: string;
  target_statuses?: Record<string, {type?: string; status?: string; [key: string]: unknown}>;
  result?: {status: string; targets: Array<{target_id: string; target_kind?: string; gene_id?: string; gene_symbol?: string; event_id?: string; status: string; report?: string; error?: string; reason?: string;
    sashimi?: Array<{status: string; reason?: string}>; [key: string]: unknown}>; [key: string]: unknown} | null;
  imported?: boolean;
};
