import { startTransition, useDeferredValue, useEffect, useMemo, useRef, useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import {
  Archive,
  Braces,
  DatabaseZap,
  FolderOpen,
  FolderX,
  Play,
  RefreshCw,
  Rows3,
} from "lucide-react";
import { api } from "./api";
import { ArtifactViewer } from "./components/ArtifactViewer";
import { ContextBar } from "./components/ContextBar";
import { Inspector } from "./components/Inspector";
import { JobDrawer } from "./components/JobDrawer";
import { ProjectGate } from "./components/ProjectGate";
import { SelectionBasket } from "./components/SelectionBasket";
import { WorkflowControls } from "./components/WorkflowControls";
import { WorkspaceNav } from "./components/WorkspaceNav";
import type { Artifact, ProjectPayload, WorkflowDefinition } from "./types";

const RECENT_MANIFEST_KEY = "katana-next:recent-manifest";
const COMPARISON_WORKSPACES = new Set(["quality", "landscape", "differential", "prioritize", "validate", "programs"]);

function defaultsFor(workflow?: WorkflowDefinition) {
  return Object.fromEntries((workflow?.parameters ?? []).map((parameter) => [parameter.key, parameter.default]));
}

function artifactPriority(artifact: Artifact) {
  const extension = artifact.path.split(".").pop()?.toLowerCase();
  if (artifact.kind === "plot_spec") return 0;
  if (["png", "svg", "jpg", "jpeg", "webp", "html", "pdf"].includes(extension ?? "")) return 1;
  if (artifact.kind === "table") return 2;
  return 3;
}

export default function App() {
  const queryClient = useQueryClient();
  const [project, setProject] = useState<ProjectPayload>();
  const [activeWorkspace, setActiveWorkspace] = useState("project");
  const [activeWorkflowId, setActiveWorkflowId] = useState("");
  const [selectedContrasts, setSelectedContrasts] = useState<Set<string>>(new Set());
  const [parameterStore, setParameterStore] = useState<Record<string, Record<string, unknown>>>({});
  const [activeArtifactId, setActiveArtifactId] = useState("");
  const [focusedRow, setFocusedRow] = useState<Record<string, unknown>>();
  const [selectedIds, setSelectedIds] = useState<Set<string>>(new Set());
  const completedJobs = useRef(new Set<string>());
  const automaticOpenAttempted = useRef(false);
  const recentManifest = localStorage.getItem(RECENT_MANIFEST_KEY) ?? "";

  const workspaces = useQuery({ queryKey: ["workspaces"], queryFn: api.workspaces });
  const workflows = useQuery({ queryKey: ["workflows"], queryFn: api.workflows });

  function acceptProject(payload: ProjectPayload) {
    localStorage.setItem(RECENT_MANIFEST_KEY, payload.snapshot.manifest_path);
    startTransition(() => {
      setProject(payload);
      setSelectedContrasts(new Set(payload.manifest.contrasts.map((contrast) => contrast.contrast_id)));
      setSelectedIds(new Set());
      setFocusedRow(undefined);
      setActiveWorkspace("project");
    });
  }

  const openProject = useMutation({ mutationFn: api.openProject, onSuccess: acceptProject });
  const migrateProject = useMutation({ mutationFn: api.migrateProject, onSuccess: acceptProject });

  useEffect(() => {
    if (!recentManifest || automaticOpenAttempted.current) return;
    automaticOpenAttempted.current = true;
    openProject.mutate(recentManifest);
  }, [recentManifest]);

  const workspaceWorkflows = useMemo(
    () => (workflows.data ?? []).filter((workflow) => workflow.workspace_id === activeWorkspace && workflow.enabled),
    [activeWorkspace, workflows.data],
  );

  useEffect(() => {
    if (!workspaceWorkflows.length) return;
    if (!workspaceWorkflows.some((workflow) => workflow.workflow_id === activeWorkflowId)) {
      setActiveWorkflowId(workspaceWorkflows[0].workflow_id);
    }
  }, [activeWorkflowId, workspaceWorkflows]);

  const activeWorkflow = (workflows.data ?? []).find((workflow) => workflow.workflow_id === activeWorkflowId);
  const projectId = project?.snapshot.project_id ?? "";
  const activeParameters = activeWorkflow ? (parameterStore[activeWorkflow.workflow_id] ?? defaultsFor(activeWorkflow)) : {};
  const requiresComparisons = Boolean(activeWorkflow && project?.manifest.contrasts.length && COMPARISON_WORKSPACES.has(activeWorkflow.workspace_id));
  const resolvedParameters = useMemo(
    () => ({
      ...activeParameters,
      ...(requiresComparisons ? { contrast_ids: [...selectedContrasts].sort() } : {}),
    }),
    [activeParameters, requiresComparisons, selectedContrasts],
  );
  const deferredParameterSignature = useDeferredValue(JSON.stringify(resolvedParameters));
  const deferredParameters = useMemo(
    () => JSON.parse(deferredParameterSignature) as Record<string, unknown>,
    [deferredParameterSignature],
  );

  const artifactQuery = useQuery({
    queryKey: ["workflow-cache", projectId, activeWorkflowId, deferredParameterSignature],
    queryFn: () => api.workflowCache(projectId, activeWorkflowId, deferredParameters),
    enabled: Boolean(projectId && activeWorkflowId),
  });
  const jobsQuery = useQuery({
    queryKey: ["jobs", projectId],
    queryFn: () => api.jobs(projectId),
    enabled: Boolean(projectId),
    refetchInterval: (query) => query.state.data?.some((job) => ["queued", "running"].includes(job.status)) ? 700 : false,
  });
  const selectionsQuery = useQuery({
    queryKey: ["selections", projectId],
    queryFn: () => api.selections(projectId),
    enabled: Boolean(projectId),
  });

  useEffect(() => {
    for (const job of jobsQuery.data ?? []) {
      if (!["complete", "failed", "cancelled", "cache_hit"].includes(job.status) || completedJobs.current.has(job.job_id)) continue;
      completedJobs.current.add(job.job_id);
      void queryClient.invalidateQueries({ queryKey: ["workflow-cache", projectId, job.analysis_id] });
    }
  }, [jobsQuery.data, projectId, queryClient]);

  const latestArtifacts = useMemo(() => {
    const compatible = artifactQuery.data?.artifacts ?? [];
    return [...compatible].sort((a, b) => artifactPriority(a) - artifactPriority(b));
  }, [artifactQuery.data]);

  useEffect(() => {
    if (!latestArtifacts.length) {
      setActiveArtifactId("");
      return;
    }
    if (!latestArtifacts.some((artifact) => artifact.artifact_id === activeArtifactId)) {
      setActiveArtifactId(latestArtifacts[0].artifact_id);
    }
    setFocusedRow(undefined);
    setSelectedIds(new Set());
  }, [activeWorkflowId, latestArtifacts]);

  const activeArtifact = latestArtifacts.find((artifact) => artifact.artifact_id === activeArtifactId) ?? latestArtifacts[0];
  const activeJob = (jobsQuery.data ?? []).find((job) => job.analysis_id === activeWorkflowId && ["queued", "running"].includes(job.status));
  const lastWorkflowJob = (jobsQuery.data ?? []).find((job) => job.analysis_id === activeWorkflowId);
  const canRun = Boolean(activeWorkflow && !activeJob && (!requiresComparisons || selectedContrasts.size > 0));

  const runWorkflow = useMutation({
    mutationFn: (force: boolean) => {
      if (!activeWorkflow) throw new Error("Choose a workflow first.");
      return api.runWorkflow(projectId, activeWorkflow.workflow_id, resolvedParameters, force);
    },
    onSuccess: () => void queryClient.invalidateQueries({ queryKey: ["jobs", projectId] }),
  });

  const openOutput = useMutation({ mutationFn: () => api.openOutput(projectId) });
  const saveSelection = useMutation({
    mutationFn: () => {
      if (!activeArtifact || activeArtifact.kind !== "table") throw new Error("Choose a table first.");
      return api.createSelection(projectId, {
        name: `${activeArtifact.entity_type} selection · ${new Date().toLocaleDateString()}`,
        entity_type: activeArtifact.entity_type,
        source_artifact_id: activeArtifact.artifact_id,
        entity_ids: [...selectedIds],
        filters: {},
      });
    },
    onSuccess: () => {
      setSelectedIds(new Set());
      void queryClient.invalidateQueries({ queryKey: ["selections", projectId] });
    },
  });

  function closeProject() {
    localStorage.removeItem(RECENT_MANIFEST_KEY);
    startTransition(() => {
      setProject(undefined);
      setFocusedRow(undefined);
      setSelectedIds(new Set());
      setParameterStore({});
    });
  }

  function updateParameter(key: string, value: unknown) {
    if (!activeWorkflow) return;
    setParameterStore((current) => ({
      ...current,
      [activeWorkflow.workflow_id]: { ...(current[activeWorkflow.workflow_id] ?? defaultsFor(activeWorkflow)), [key]: value },
    }));
  }

  function resetParameters() {
    if (!activeWorkflow) return;
    setParameterStore((current) => ({ ...current, [activeWorkflow.workflow_id]: defaultsFor(activeWorkflow) }));
  }

  function toggleContrast(contrastId: string) {
    setSelectedContrasts((current) => {
      const next = new Set(current);
      if (next.has(contrastId)) next.delete(contrastId); else next.add(contrastId);
      return next;
    });
  }

  function toggleSelection(id: string) {
    setSelectedIds((current) => {
      const next = new Set(current);
      if (next.has(id)) next.delete(id); else next.add(id);
      return next;
    });
  }

  if (!project) {
    return (
      <ProjectGate
        defaultPath={recentManifest}
        loading={openProject.isPending || migrateProject.isPending}
        error={openProject.error?.message ?? migrateProject.error?.message}
        onOpen={(path) => openProject.mutate(path)}
        onMigrate={(request) => migrateProject.mutate(request)}
      />
    );
  }

  const conditions = [...new Set(project.manifest.samples.map((sample) => sample.condition))];
  return (
    <div className="app-frame">
      <WorkspaceNav workspaces={workspaces.data ?? []} activeId={activeWorkspace} onChange={setActiveWorkspace} />
      <section className="app-main">
        <header className="topbar">
          <div><span className="eyebrow">Evidence-led RNA analysis</span><h1>{workspaces.data?.find((item) => item.workspace_id === activeWorkspace)?.label ?? "Katana"}</h1></div>
          <div className="topbar-actions">
            <span className="version-badge"><Braces size={15} />2.0</span>
            <button type="button" className="ghost-action" onClick={() => openOutput.mutate()}><FolderOpen size={16} />Output folder</button>
            <button type="button" className="ghost-action" onClick={closeProject}><FolderX size={16} />Close</button>
          </div>
        </header>

        <ContextBar project={project} selectedContrastCount={selectedContrasts.size} workflowLabel={activeWorkflow?.label ?? "Choose workflow"} />

        <div className="workspace-grid">
          <main className="workbench">
            <div className="workflow-switcher" role="tablist" aria-label="Workspace analyses">
              {workspaceWorkflows.map((workflow) => (
                <button
                  type="button"
                  role="tab"
                  aria-selected={activeWorkflowId === workflow.workflow_id}
                  className={activeWorkflowId === workflow.workflow_id ? "active" : ""}
                  key={workflow.workflow_id}
                  onClick={() => setActiveWorkflowId(workflow.workflow_id)}
                >
                  <span>{workflow.label}</span><small>{workflow.analysis_unit}</small>
                </button>
              ))}
            </div>

            {activeWorkflow && (
              <>
                <div className="workbench-heading">
                  <div><span className="eyebrow">{activeWorkflow.workflow_id}</span><h2>{activeWorkflow.question}</h2><p>{activeWorkflow.output_summary}</p></div>
                  <div className="run-actions">
                    <button type="button" className="secondary-action" disabled={!canRun || runWorkflow.isPending} onClick={() => runWorkflow.mutate(true)}><RefreshCw size={16} />Rebuild</button>
                    <button type="button" className="primary-action" disabled={!canRun || runWorkflow.isPending} onClick={() => runWorkflow.mutate(false)}>
                      {activeJob ? <RefreshCw className="spin" size={16} /> : <Play size={16} />}{activeJob ? "Running" : "Run analysis"}
                    </button>
                  </div>
                </div>

                {requiresComparisons && selectedContrasts.size === 0 && <div className="warning-strip">Select at least one comparison. An empty selection never silently means all.</div>}
                {(runWorkflow.error || lastWorkflowJob?.error) && <div className="warning-strip error-strip">{runWorkflow.error?.message ?? lastWorkflowJob?.error}</div>}

                <WorkflowControls
                  workflow={activeWorkflow}
                  contrasts={project.manifest.contrasts}
                  selectedContrasts={selectedContrasts}
                  parameters={activeParameters}
                  selections={selectionsQuery.data ?? []}
                  conditions={conditions}
                  onToggleContrast={toggleContrast}
                  onSelectAllContrasts={() => setSelectedContrasts(new Set(project.manifest.contrasts.map((item) => item.contrast_id)))}
                  onClearContrasts={() => setSelectedContrasts(new Set())}
                  onParameterChange={updateParameter}
                  onReset={resetParameters}
                />

                <div className="metric-strip">
                  <article><Archive size={18} /><div><small>Cached outputs</small><strong>{latestArtifacts.length}</strong></div></article>
                  <article><DatabaseZap size={18} /><div><small>Last access</small><strong>{lastWorkflowJob?.cache_status === "hit" ? "Compatible cache hit" : lastWorkflowJob?.status ?? "Not run this session"}</strong></div></article>
                  <article><Rows3 size={18} /><div><small>Analysis unit</small><strong>{activeWorkflow.analysis_unit}</strong></div></article>
                </div>

                <section className="result-card">
                  <div className="artifact-nav">
                    <div><span className="eyebrow">Latest compatible run</span><strong>Results</strong></div>
                    <div className="artifact-tabs">
                      {latestArtifacts.map((artifact) => (
                        <button type="button" className={activeArtifact?.artifact_id === artifact.artifact_id ? "active" : ""} key={artifact.artifact_id} onClick={() => { setActiveArtifactId(artifact.artifact_id); setFocusedRow(undefined); setSelectedIds(new Set()); }}>
                          {artifact.name}<small>{artifact.kind.replace("_", " ")}{artifact.row_count != null ? ` · ${artifact.row_count.toLocaleString()}` : ""}</small>
                        </button>
                      ))}
                    </div>
                    {activeArtifact && <span className="result-provenance">run {activeArtifact.run_fingerprint.slice(0, 10)}</span>}
                  </div>
                  <div className="result-view">
                    <ArtifactViewer projectId={projectId} artifact={activeArtifact} selectedIds={selectedIds} onToggle={toggleSelection} onFocus={setFocusedRow} />
                  </div>
                </section>
              </>
            )}
          </main>
          <Inspector workflow={activeWorkflow} artifact={activeArtifact} row={focusedRow} />
        </div>

        <SelectionBasket count={selectedIds.size} saving={saveSelection.isPending} onSave={() => saveSelection.mutate()} onClear={() => setSelectedIds(new Set())} />
        <JobDrawer jobs={jobsQuery.data ?? []} />
      </section>
    </div>
  );
}
