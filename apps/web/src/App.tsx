import { startTransition, useEffect, useRef, useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import {
  Archive,
  ArrowRight,
  Braces,
  ChartNoAxesColumnIncreasing,
  DatabaseZap,
  FolderX,
  Play,
  RefreshCw,
  Table2,
} from "lucide-react";
import { api } from "./api";
import { ArtifactTable } from "./components/ArtifactTable";
import { ContextBar } from "./components/ContextBar";
import { Inspector } from "./components/Inspector";
import { JobDrawer } from "./components/JobDrawer";
import { PlotViewer } from "./components/PlotViewer";
import { ProjectGate } from "./components/ProjectGate";
import { SelectionBasket } from "./components/SelectionBasket";
import { WorkspaceNav } from "./components/WorkspaceNav";
import type { Artifact, ProjectPayload } from "./types";

const RECENT_MANIFEST_KEY = "katana-next:recent-manifest";

export default function App() {
  const queryClient = useQueryClient();
  const [project, setProject] = useState<ProjectPayload>();
  const [activeWorkspace, setActiveWorkspace] = useState("project");
  const [contrastId, setContrastId] = useState("");
  const [view, setView] = useState<"plot" | "table">("plot");
  const [focusedRow, setFocusedRow] = useState<Record<string, unknown>>();
  const [selectedIds, setSelectedIds] = useState<Set<string>>(new Set());
  const completedJobs = useRef(new Set<string>());
  const automaticOpenAttempted = useRef(false);
  const recentManifest = localStorage.getItem(RECENT_MANIFEST_KEY) ?? "";

  const workspaces = useQuery({ queryKey: ["workspaces"], queryFn: api.workspaces });
  const openProject = useMutation({
    mutationFn: api.openProject,
    onSuccess: (payload, manifestPath) => {
      localStorage.setItem(RECENT_MANIFEST_KEY, manifestPath);
      startTransition(() => {
        setProject(payload);
        setSelectedIds(new Set());
        setFocusedRow(undefined);
      });
    },
  });

  useEffect(() => {
    if (!recentManifest || automaticOpenAttempted.current) return;
    automaticOpenAttempted.current = true;
    openProject.mutate(recentManifest);
  }, [recentManifest]);

  const projectId = project?.snapshot.project_id ?? "";
  const artifactQuery = useQuery({
    queryKey: ["artifacts", projectId],
    queryFn: () => api.artifacts(projectId),
    enabled: Boolean(projectId),
  });
  const jobsQuery = useQuery({
    queryKey: ["jobs", projectId],
    queryFn: () => api.jobs(projectId),
    enabled: Boolean(projectId),
    refetchInterval: (query) => {
      const jobs = query.state.data;
      return jobs?.some((job) => job.status === "queued" || job.status === "running")
        ? 750
        : false;
    },
  });
  const selectionsQuery = useQuery({
    queryKey: ["selections", projectId],
    queryFn: () => api.selections(projectId),
    enabled: Boolean(projectId),
  });

  useEffect(() => {
    for (const job of jobsQuery.data ?? []) {
      if (!["complete", "failed", "cancelled", "cache_hit"].includes(job.status)) continue;
      if (completedJobs.current.has(job.job_id)) continue;
      completedJobs.current.add(job.job_id);
      void queryClient.invalidateQueries({ queryKey: ["artifacts", projectId] });
    }
  }, [jobsQuery.data, projectId, queryClient]);

  const runWorkflow = useMutation({
    mutationFn: (force: boolean) => api.runPhase1(projectId, force),
    onSuccess: () => void queryClient.invalidateQueries({ queryKey: ["jobs", projectId] }),
  });

  const artifacts = artifactQuery.data ?? [];
  const tableArtifact = artifacts.find((artifact) => artifact.kind === "table");
  const plotArtifact = artifacts.find((artifact) => artifact.kind === "plot_spec");
  const saveSelection = useMutation({
    mutationFn: () => {
      if (!tableArtifact) throw new Error("No source table is available.");
      return api.createSelection(projectId, {
        name: `Review events ${new Date().toLocaleDateString()}`,
        entity_type: tableArtifact.entity_type,
        source_artifact_id: tableArtifact.artifact_id,
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
    });
  }

  function toggleSelection(id: string) {
    setSelectedIds((current) => {
      const next = new Set(current);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  }

  if (!project) {
    return (
      <ProjectGate
        defaultPath={recentManifest}
        loading={openProject.isPending}
        error={openProject.error?.message}
        onOpen={(path) => openProject.mutate(path)}
      />
    );
  }

  const running = (jobsQuery.data ?? []).some((job) => job.status === "queued" || job.status === "running");
  const activeArtifact: Artifact | undefined = view === "plot" ? plotArtifact : tableArtifact;

  return (
    <div className="app-frame">
      <WorkspaceNav
        workspaces={workspaces.data ?? []}
        activeId={activeWorkspace}
        onChange={setActiveWorkspace}
      />
      <section className="app-main">
        <header className="topbar">
          <div>
            <span className="eyebrow">Platform foundation</span>
            <h1>Project evidence workspace</h1>
          </div>
          <div className="topbar-actions">
            <span className="alpha-badge"><Braces size={15} />2.0 alpha</span>
            <button type="button" className="ghost-action" onClick={closeProject}><FolderX size={16} />Close project</button>
          </div>
        </header>

        <ContextBar project={project} contrastId={contrastId} onContrastChange={setContrastId} />

        <div className="workspace-grid">
          <main className="workbench">
            <div className="workbench-heading">
              <div>
                <span className="eyebrow">Project / Artifact contract</span>
                <h2>Input and cache verification</h2>
                <p>Opening this page reads persisted artifacts only. Run is explicit and isolated.</p>
              </div>
              <div className="run-actions">
                <button
                  type="button"
                  className="primary-action"
                  disabled={running || runWorkflow.isPending}
                  onClick={() => runWorkflow.mutate(false)}
                >
                  {running ? <RefreshCw className="spin" size={17} /> : <Play size={17} />}
                  {running ? "Running check" : artifacts.length ? "Run / use cache" : "Run platform check"}
                </button>
                {artifacts.length > 0 && (
                  <button
                    type="button"
                    className="icon-action"
                    title="Force rebuild"
                    aria-label="Force rebuild"
                    disabled={running}
                    onClick={() => runWorkflow.mutate(true)}
                  ><RefreshCw size={17} /></button>
                )}
              </div>
            </div>

            {project.snapshot.warnings.length > 0 && (
              <div className="warning-strip">{project.snapshot.warnings.join(" ")}</div>
            )}

            <section className="metric-strip">
              <article><DatabaseZap size={18} /><span><small>Cache state</small><strong>{artifacts.length ? `${artifacts.length} artifacts ready` : "No cached result"}</strong></span></article>
              <article><Archive size={18} /><span><small>Saved selections</small><strong>{selectionsQuery.data?.length ?? project.selection_count}</strong></span></article>
              <article><Table2 size={18} /><span><small>Samples / contrasts</small><strong>{project.snapshot.sample_count} / {project.snapshot.contrast_count}</strong></span></article>
            </section>

            <section className="result-card">
              <div className="result-tabs" role="tablist">
                <button role="tab" aria-selected={view === "plot"} className={view === "plot" ? "active" : ""} onClick={() => setView("plot")}>
                  <ChartNoAxesColumnIncreasing size={16} />Plot
                </button>
                <button role="tab" aria-selected={view === "table"} className={view === "table" ? "active" : ""} onClick={() => setView("table")}>
                  <Table2 size={16} />Table
                </button>
                <span className="result-provenance">
                  {activeArtifact ? `Cache ${activeArtifact.run_fingerprint.slice(0, 10)}` : "Awaiting explicit run"}
                </span>
              </div>
              <div className="result-view">
                {view === "plot" ? (
                  <PlotViewer projectId={projectId} artifact={plotArtifact} />
                ) : (
                  <ArtifactTable
                    projectId={projectId}
                    artifact={tableArtifact}
                    selectedIds={selectedIds}
                    onToggle={toggleSelection}
                    onFocus={setFocusedRow}
                  />
                )}
              </div>
              {!artifacts.length && !running && (
                <div className="first-run-callout">
                  <span>Platform acceptance fixture</span>
                  <p>This run creates one Parquet table and one responsive Vega-Lite plot. It does not launch DEG, DAS, CoSplice, or any legacy module.</p>
                  <button type="button" onClick={() => runWorkflow.mutate(false)}>Run isolated check<ArrowRight size={16} /></button>
                </div>
              )}
            </section>
          </main>

          <Inspector artifact={activeArtifact} row={focusedRow} />
        </div>

        <SelectionBasket
          count={selectedIds.size}
          saving={saveSelection.isPending}
          onSave={() => saveSelection.mutate()}
          onClear={() => setSelectedIds(new Set())}
        />
        <JobDrawer jobs={jobsQuery.data ?? []} />
      </section>
    </div>
  );
}
