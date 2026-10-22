import { startTransition, useDeferredValue, useEffect, useMemo, useRef, useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import {
  Archive,
  DatabaseZap,
  FolderOpen,
  FolderX,
  Play,
  RefreshCw,
  Rows3,
  Info,
} from "lucide-react";
import { api, type MotifTraceMode } from "./api";
import { displayModuleName, displayArtifactName } from "./branding";
import { DEFAULT_DAS_EVENT_SET, isDasEventSet, isDasSwitchable, parametersWithDasSelection, type DasEventSet } from "./dasSelection";
import { ArtifactViewer } from "./components/ArtifactViewer";
import { ContextBar } from "./components/ContextBar";
import { Inspector } from "./components/Inspector";
import { JobDrawer } from "./components/JobDrawer";
import { ProjectGate } from "./components/ProjectGate";
import { GeneFilterPanel } from "./components/GeneFilterPanel";
import { SelectionBasket } from "./components/SelectionBasket";
import { WorkflowControls } from "./components/WorkflowControls";
import { WorkspaceNav } from "./components/WorkspaceNav";
import { RunPlanPanel } from "./components/RunPlanPanel";
import { RbpEvidenceGuide } from "./components/RbpEvidenceGuide";
import { MotifTraceWorkbench, motifTraceModeParameters } from "./components/MotifTraceWorkbench";
import type { Artifact, ProjectPayload, WorkflowDefinition } from "./types";

const RECENT_MANIFEST_KEY = "katana-next:recent-manifest";

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
  const [showAdvancedAnalyses, setShowAdvancedAnalyses] = useState(false);
  const [selectedContrasts, setSelectedContrasts] = useState<Set<string>>(new Set());
  const [parameterStore, setParameterStore] = useState<Record<string, Record<string, unknown>>>({});
  const [dasEventSet, setDasEventSet] = useState<DasEventSet>(DEFAULT_DAS_EVENT_SET);
  const [activeArtifactId, setActiveArtifactId] = useState("");
  const [focusedRow, setFocusedRow] = useState<Record<string, unknown>>();
  const [detailOpen, setDetailOpen] = useState(false);
  const [aboutOpen, setAboutOpen] = useState(false);
  const [selectedIds, setSelectedIds] = useState<Set<string>>(new Set());
  const completedJobs = useRef(new Set<string>());
  const workspaceScroll = useRef<HTMLElement>(null);
  const automaticOpenAttempted = useRef(false);
  const recentManifest = localStorage.getItem(RECENT_MANIFEST_KEY) ?? "";

  useEffect(() => {
    workspaceScroll.current?.scrollTo({ top: 0, left: 0 });
  }, [activeWorkspace, activeWorkflowId]);

  const workspaces = useQuery({ queryKey: ["workspaces"], queryFn: api.workspaces });
  const workflows = useQuery({ queryKey: ["workflows"], queryFn: api.workflows });
  const motifTraceCatalog = useQuery({ queryKey: ["motiftrace-modes"], queryFn: api.motifTraceModes,
    enabled: activeWorkspace === "prioritize" });

  function acceptProject(payload: ProjectPayload) {
    localStorage.setItem(RECENT_MANIFEST_KEY, payload.snapshot.manifest_path);
    startTransition(() => {
      if (project?.snapshot.manifest_path !== payload.snapshot.manifest_path) {
        setDasEventSet(DEFAULT_DAS_EVENT_SET);
        setParameterStore({});
      }
      setProject(payload);
      setSelectedContrasts(new Set(payload.manifest.contrasts.map((contrast) => contrast.contrast_id)));
      setSelectedIds(new Set());
      setFocusedRow(undefined);
      setActiveWorkspace("project");
    });
  }

  const openProject = useMutation({ mutationFn: api.openProject, onSuccess: acceptProject });
  const migrateProject = useMutation({ mutationFn: api.migrateProject, onSuccess: acceptProject });
  const configureResources = useMutation({ mutationFn: api.configureResources, onSuccess: acceptProject });
  const configureGeneFilter = useMutation({
    mutationFn: (settings: import("./components/GeneFilterPanel").GeneFilterSettings) => api.configureGeneFilter(project!.snapshot.project_id, settings),
    onSuccess: acceptProject,
  });

  useEffect(() => {
    if (!recentManifest || automaticOpenAttempted.current) return;
    automaticOpenAttempted.current = true;
    openProject.mutate(recentManifest);
  }, [recentManifest]);

  const workspaceWorkflows = useMemo(
    () => (workflows.data ?? []).filter((workflow) => workflow.workspace_id === activeWorkspace && workflow.enabled && (showAdvancedAnalyses || workflow.workflow_id === activeWorkflowId || !workflow.presentation_role || workflow.presentation_role === "primary")),
    [activeWorkspace, activeWorkflowId, workflows.data, showAdvancedAnalyses],
  );

  useEffect(() => {
    if (!workspaceWorkflows.length) return;
    if (!workspaceWorkflows.some((workflow) => workflow.workflow_id === activeWorkflowId)) {
      setActiveWorkflowId(workspaceWorkflows[0].workflow_id);
    }
  }, [activeWorkflowId, workspaceWorkflows]);

  const activeWorkflow = (workflows.data ?? []).find((workflow) => workflow.workflow_id === activeWorkflowId);
  const projectId = project?.snapshot.project_id ?? "";
  const sourceArtifacts = useQuery({ queryKey: ["source-artifacts", projectId], queryFn: () => api.artifacts(projectId), enabled: Boolean(projectId) });
  const scConditions = useQuery({ queryKey: ["single-cell-conditions", projectId], queryFn: () => api.singleCellConditions(projectId), enabled: Boolean(projectId) && activeWorkspace === "single-cell" });
  const activeParameters = parametersWithDasSelection(activeWorkflow,
    activeWorkflow ? (parameterStore[activeWorkflow.workflow_id] ?? defaultsFor(activeWorkflow)) : {}, dasEventSet);
  const requiresComparisons = Boolean(
    activeWorkflow?.supports_contrast_selection && project?.manifest.contrasts.length,
  );
  const allComparisonsSelected = Boolean(
    project?.manifest.contrasts.length
      && selectedContrasts.size === project.manifest.contrasts.length
      && project.manifest.contrasts.every((contrast) => selectedContrasts.has(contrast.contrast_id)),
  );
  const resolvedParameters = useMemo(
    () => ({
      ...activeParameters,
      ...(requiresComparisons && !allComparisonsSelected
        ? { contrast_ids: [...selectedContrasts].sort() }
        : {}),
    }),
    [activeParameters, allComparisonsSelected, requiresComparisons, selectedContrasts],
  );
  const deferredParameterSignature = useDeferredValue(JSON.stringify(resolvedParameters));

  const artifactQuery = useQuery({
    queryKey: ["workflow-cache", projectId, activeWorkflowId, deferredParameterSignature],
    queryFn: ({ queryKey }) => {
      const [, cacheProjectId, cacheWorkflowId, parameterSignature] = queryKey as [string, string, string, string];
      return api.workflowCache(
        cacheProjectId,
        cacheWorkflowId,
        JSON.parse(parameterSignature) as Record<string, unknown>,
      );
    },
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
      void queryClient.invalidateQueries({ queryKey: ["source-artifacts", projectId] });
      void queryClient.invalidateQueries({ queryKey: ["module-evidence", projectId] });
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
  const savedSettingsArtifact = (sourceArtifacts.data ?? []).find((artifact) =>
    artifact.analysis_id === activeWorkflowId
    && artifact.metadata.project_snapshot_fingerprint === project?.snapshot.manifest_fingerprint,
  );
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
      setDasEventSet(DEFAULT_DAS_EVENT_SET);
    });
  }

  function updateParameter(key: string, value: unknown) {
    if (!activeWorkflow) return;
    if (key === "event_set" && isDasSwitchable(activeWorkflow) && isDasEventSet(value)) setDasEventSet(value);
    setParameterStore((current) => ({
      ...current,
      [activeWorkflow.workflow_id]: { ...(current[activeWorkflow.workflow_id] ?? defaultsFor(activeWorkflow)), [key]: value },
    }));
  }

  function selectMotifTraceMode(mode: MotifTraceMode) {
    const target=(workflows.data ?? []).find(workflow=>workflow.workflow_id===mode.workflow_id);
    if (!target) return;
    setParameterStore(current=>({...current,[mode.workflow_id]:motifTraceModeParameters(mode,target,activeParameters,current[mode.workflow_id])}));
    setActiveWorkflowId(mode.workflow_id);
  }

  function resetParameters() {
    if (!activeWorkflow) return;
    if (isDasSwitchable(activeWorkflow)) setDasEventSet(DEFAULT_DAS_EVENT_SET);
    setParameterStore((current) => ({ ...current, [activeWorkflow.workflow_id]: defaultsFor(activeWorkflow) }));
  }

  function restoreSavedParameters() {
    if (!activeWorkflow || !savedSettingsArtifact) return;
    if (isDasSwitchable(activeWorkflow) && isDasEventSet(savedSettingsArtifact.parameters.event_set)) {
      setDasEventSet(savedSettingsArtifact.parameters.event_set);
    }
    setParameterStore((current) => ({
      ...current,
      [activeWorkflow.workflow_id]: { ...defaultsFor(activeWorkflow), ...savedSettingsArtifact.parameters },
    }));
    if (requiresComparisons) {
      setSelectedContrasts(new Set(savedSettingsArtifact.contrast_ids.length
        ? savedSettingsArtifact.contrast_ids
        : project?.manifest.contrasts.map((contrast) => contrast.contrast_id) ?? []));
    }
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
        loading={openProject.isPending || migrateProject.isPending || configureResources.isPending}
        error={openProject.error?.message ?? migrateProject.error?.message ?? configureResources.error?.message}
        onOpen={(path) => openProject.mutate(path)}
        onMigrate={(request) => migrateProject.mutate(request)}
        onConfigureResources={(request) => configureResources.mutate(request)}
      />
    );
  }

  const conditions = activeWorkspace === "single-cell" ? (scConditions.data ?? []) : [...new Set(project.manifest.samples.map((sample) => sample.condition))];
  return (
    <div className="app-frame">
      <WorkspaceNav workspaces={workspaces.data ?? []} activeId={activeWorkspace} onChange={setActiveWorkspace} />
      <section ref={workspaceScroll} className="app-main" tabIndex={0} aria-label="Analysis workspace">
        <header className="topbar">
          <div><span className="eyebrow">RNA splicing analysis</span><h1>{workspaces.data?.find((item) => item.workspace_id === activeWorkspace)?.label ?? "Katana"}</h1></div>
          <div className="topbar-actions">
            <button type="button" className="ghost-action" aria-expanded={detailOpen} aria-controls="scientific-detail" onClick={() => setDetailOpen(!detailOpen)}>{detailOpen ? "Hide methods" : "Methods & details"}</button>
            <button type="button" className="ghost-action" aria-expanded={aboutOpen} onClick={() => setAboutOpen(!aboutOpen)}><Info size={16} />About</button>
            <button type="button" className="ghost-action" onClick={() => openOutput.mutate()}><FolderOpen size={16} />Output folder</button>
            <button type="button" className="ghost-action" onClick={closeProject}><FolderX size={16} />Close</button>
          </div>
        </header>
        {aboutOpen && <section className="software-info" aria-label="About Katana"><strong>Katana</strong><p>Bulk and single-cell RNA splicing analysis, visualization and reporting.</p><small>Software version: 2.0.0-beta.2. Browser and desktop use the same local analysis engine. Input files remain on this computer.</small></section>}

        <ContextBar project={project} selectedContrastCount={selectedContrasts.size} workflowLabel={activeWorkflow?.label ?? "Choose workflow"} />
        {activeWorkspace === "project" && <GeneFilterPanel key={project.snapshot.manifest_fingerprint} project={project}
          busy={configureGeneFilter.isPending} error={configureGeneFilter.error?.message}
          onApply={settings=>configureGeneFilter.mutate(settings)}/>}
        {activeWorkspace === "project" && project.snapshot.warnings.length > 0 && (
          <details className="warning-strip" open>
            <summary>Project configuration needs attention</summary>
            {project.snapshot.warnings.map((warning) => <p key={warning}>{warning}</p>)}
            <small>Close this project and use GO / RBP resources to attach a prepared bundle to a new project copy.</small>
          </details>
        )}

        <div className={`workspace-grid${detailOpen ? "" : " detail-closed"}`}>
          <main className="workbench">
            {activeWorkspace === "prioritize" && <MotifTraceWorkbench projectId={projectId}
              catalog={motifTraceCatalog.data} loading={motifTraceCatalog.isLoading} error={motifTraceCatalog.error?.message}
              workflows={workflows.data ?? []}
              artifacts={(sourceArtifacts.data ?? []).filter(a=>a.metadata.project_snapshot_fingerprint===project.snapshot.manifest_fingerprint)}
              activeWorkflowId={activeWorkflowId} parameters={activeParameters} onSelect={selectMotifTraceMode} />}
            {activeWorkspace === "prioritize" && (
              <RbpEvidenceGuide workflows={workflows.data ?? []}
                artifacts={(sourceArtifacts.data ?? []).filter(a => a.metadata.project_snapshot_fingerprint === project.snapshot.manifest_fingerprint)}
                activeWorkflowId={activeWorkflowId} parameters={activeParameters}
                onSelect={(id, patch) => {
                  setActiveWorkflowId(id);
                  if (patch) setParameterStore(current => ({...current, [id]: {
                    ...defaultsFor((workflows.data ?? []).find(w => w.workflow_id === id)), ...current[id], ...patch,
                  }}));
                }} />
            )}
            {activeWorkspace === "prioritize" && (
              <label className="advanced-analyses-toggle">
                <input type="checkbox" checked={showAdvancedAnalyses} onChange={(event) => setShowAdvancedAnalyses(event.target.checked)} />
                Show additional RBP method steps and legacy analyses
              </label>
            )}
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
                  <span>{workflow.label}</span>
                </button>
              ))}
            </div>

            {activeWorkflow && (
              <>
                <div className="workbench-heading">
                  <div><span className="eyebrow">Analysis</span><h2>{activeWorkflow.label}</h2><p>{activeWorkflow.question}</p><p><strong>Outputs: </strong>{activeWorkflow.output_summary}</p></div>
                  <div className="run-actions">
                    <button type="button" className="secondary-action" disabled={!canRun || runWorkflow.isPending} onClick={() => runWorkflow.mutate(true)}><RefreshCw size={16} />Re-run analysis</button>
                    <button type="button" className="primary-action" disabled={!canRun || runWorkflow.isPending} onClick={() => runWorkflow.mutate(false)}>
                      {activeJob ? <RefreshCw className="spin" size={16} /> : <Play size={16} />}{activeJob ? "Running" : "Run analysis"}
                    </button>
                  </div>
                </div>

                {requiresComparisons && selectedContrasts.size === 0 && <div className="warning-strip">Select at least one comparison. An empty selection never silently means all.</div>}
                {(runWorkflow.error || lastWorkflowJob?.error) && <div className="warning-strip error-strip">{runWorkflow.error?.message ?? lastWorkflowJob?.error}</div>}
                {activeWorkflow.workflow_id === "differential.as" && <div className="run-plan-panel"><p>This module uses the completed AS read-support result selected below. Set FDR and |dPSI| in Quality Control. Read support requires IJC + SJC of at least 10 in every selected sample.</p><button type="button" className="text-action" onClick={() => {setActiveWorkspace("quality"); setActiveWorkflowId("quality.stable_as");}}>Edit DAS thresholds</button></div>}
                {artifactQuery.error && <div className="warning-strip error-strip">Could not read compatible cache: {artifactQuery.error.message}</div>}
                {!artifactQuery.isFetching && !latestArtifacts.length && savedSettingsArtifact && (
                  <div className="warning-strip">
                    <span>Saved outputs exist with different settings or an older workflow version. Restore their settings to check compatible cache. This does not run an analysis.</span>
                    <button className="secondary-button" onClick={restoreSavedParameters}>Restore saved settings</button>
                  </div>
                )}

                <WorkflowControls
                  workflow={activeWorkflow}
                  contrasts={project.manifest.contrasts}
                  selectedContrasts={selectedContrasts}
                  parameters={activeParameters}
                  selections={selectionsQuery.data ?? []}
                  conditions={conditions}
                  artifacts={(sourceArtifacts.data ?? []).filter((a) => a.metadata.project_snapshot_fingerprint === project.snapshot.manifest_fingerprint)}
                  onToggleContrast={toggleContrast}
                  onSelectAllContrasts={() => setSelectedContrasts(new Set(project.manifest.contrasts.map((item) => item.contrast_id)))}
                  onClearContrasts={() => setSelectedContrasts(new Set())}
                  onParameterChange={updateParameter}
                  onReset={resetParameters}
                />
                <RunPlanPanel key={`${projectId}:${activeWorkflow.workflow_id}:${JSON.stringify(resolvedParameters)}`}
                  projectId={projectId} workflowId={activeWorkflow.workflow_id} parameters={resolvedParameters}
                  onDesignSaved={payload => {setProject(payload); setSelectedIds(new Set()); setFocusedRow(undefined); setActiveArtifactId("");}} />

                <div className="metric-strip">
                  <article><Archive size={18} /><div><small>Cached outputs</small><strong>{artifactQuery.isFetching ? "..." : latestArtifacts.length}</strong></div></article>
                  <article><DatabaseZap size={18} /><div><small>Last access</small><strong>{lastWorkflowJob?.cache_status === "hit" ? "Compatible cache hit" : lastWorkflowJob?.status ?? "Not run this session"}</strong></div></article>
                  <article><Rows3 size={18} /><div><small>Analysis unit</small><strong>{activeWorkflow.analysis_unit}</strong></div></article>
                </div>

                <section className="result-card">
                  <div className="artifact-nav">
                    <div><span className="eyebrow">Latest compatible run</span><strong>Results</strong></div>
                    <div className="artifact-tabs">
                      {latestArtifacts.map((artifact) => (
                        <button type="button" className={activeArtifact?.artifact_id === artifact.artifact_id ? "active" : ""} key={artifact.artifact_id} onClick={() => { setActiveArtifactId(artifact.artifact_id); setFocusedRow(undefined); setSelectedIds(new Set()); }}>
                          {displayArtifactName(artifact.name)}<small>{artifact.kind.replace("_", " ")}{artifact.row_count != null ? ` · ${artifact.row_count.toLocaleString()}` : ""}</small>
                        </button>
                      ))}
                    </div>
                    {activeArtifact && <span className="result-provenance">run {activeArtifact.run_fingerprint.slice(0, 10)}</span>}
                  </div>
                  <div className="result-view">
                    <ArtifactViewer projectId={projectId} artifact={activeArtifact} selectedIds={selectedIds} onToggle={toggleSelection} onFocus={(row) => { setFocusedRow(row); setDetailOpen(true); }} />
                  </div>
                </section>
              </>
            )}
          </main>
          {detailOpen && <Inspector workflow={activeWorkflow} artifact={activeArtifact} row={focusedRow} />}
        </div>

        <SelectionBasket count={selectedIds.size} saving={saveSelection.isPending} onSave={() => saveSelection.mutate()} onClear={() => setSelectedIds(new Set())} />
        <JobDrawer jobs={jobsQuery.data ?? []} />
      </section>
    </div>
  );
}
