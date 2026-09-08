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
  workspaces: () => request<Workspace[]>("/workspaces"),
  workflows: () => request<WorkflowDefinition[]>("/workflows"),
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
  cancelJob: (jobId: string) => request<Job>(`/jobs/${encodeURIComponent(jobId)}/cancel`, { method: "POST" }),
};
