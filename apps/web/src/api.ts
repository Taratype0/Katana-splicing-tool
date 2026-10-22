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
  hub_candidate?: boolean; candidate_all_edges?: boolean; candidate_no_same_gene?: boolean; selection_reason?: string};

export type ModuleEvidence = {module: string; discovery_run: string; scope: string; tables: {
  artifact_id: string; name: string; entity_type: string; row_count: number; preview_limit: number; rows: Record<string,unknown>[]; run_fingerprint: string;
}[]};
