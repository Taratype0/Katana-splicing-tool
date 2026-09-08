import type { Artifact, Job, ProjectPayload, Selection, TablePage, Workspace } from "./types";

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
  openProject: (manifestPath: string) =>
    request<ProjectPayload>("/projects/open", {
      method: "POST",
      body: JSON.stringify({ manifest_path: manifestPath }),
    }),
  artifacts: (projectId: string) => request<Artifact[]>(`/projects/${projectId}/artifacts`),
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
    return request<TablePage>(`/projects/${projectId}/artifacts/${artifactId}/rows?${query}`);
  },
  artifactContent: (projectId: string, artifactId: string) =>
    request<Record<string, unknown>>(`/projects/${projectId}/artifacts/${artifactId}/content`),
  selections: (projectId: string) => request<Selection[]>(`/projects/${projectId}/selections`),
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
    request<Selection>(`/projects/${projectId}/selections`, {
      method: "POST",
      body: JSON.stringify(payload),
    }),
  runPhase1: (projectId: string, force = false) =>
    request<Job>(`/projects/${projectId}/workflows/phase1-profile/run`, {
      method: "POST",
      body: JSON.stringify({ parameters: {}, force }),
    }),
  jobs: (projectId: string) => request<Job[]>(`/jobs?project_id=${encodeURIComponent(projectId)}`),
};
