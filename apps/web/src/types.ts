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
  kind: "number" | "select" | "boolean" | "text" | "selection" | "condition";
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
  selection_id: string;
  project_id: string;
  name: string;
  entity_type: string;
  source_artifact_id: string;
  entity_ids: string[];
  created_at: string;
};
