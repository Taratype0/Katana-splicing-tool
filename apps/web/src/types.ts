export type Workspace = {
  workspace_id: string;
  label: string;
  purpose: string;
  icon: string;
  enabled: boolean;
  phase: number;
};

export type FileFingerprint = {
  path: string;
  exists: boolean;
  size?: number;
  digest?: string;
  hash_method: string;
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
    contrasts: Array<{ contrast_id: string; label: string }>;
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
  row_count?: number;
  metadata: Record<string, unknown>;
  created_at: string;
};

export type TablePage = {
  columns: string[];
  rows: Array<Record<string, unknown>>;
  offset: number;
  limit: number;
  total: number;
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
