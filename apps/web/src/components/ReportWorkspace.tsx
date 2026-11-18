import {useState} from "react";
import {useQuery} from "@tanstack/react-query";
import type {Artifact, ProjectPayload, TargetBasketEntry} from "../types";
import type {ApplicationRoute} from "./ApplicationTree";
import {ReportExplorer} from "./ReportExplorer";
import {ReportFigures, reportGet, reportUrl, type ReportPanel} from "./ReportFigures";
import {ReportNetwork} from "./ReportNetwork";
import {SpliceTraceExplorer} from "./SpliceTraceExplorer";
import {QuadrantExplorer, ContrastRescueExplorer} from "./QuadrantAndRescueExplorer";
import {SingleCellManifoldExplorer} from "./SingleCellManifoldExplorer";
import {GeneCardDossier} from "./GeneCardDossier";
import {PsiLandscapeLab} from "./PsiLandscapeLab";
import {LiveOraWorkbench} from "./LiveOraWorkbench";

export function ReportWorkspace({
  project,
  artifacts,
  route,
  onNavigate,
  onSetup,
  onAdd,
  onOpenGeneCard,
}: {
  project: ProjectPayload;
  artifacts: Artifact[];
  route: ApplicationRoute;
  onNavigate: (id: string) => void;
  onSetup: () => void;
  onAdd: (entry: TargetBasketEntry) => void;
  onOpenGeneCard?: (gene: string, eventId?: string, comparison?: string) => void;
}) {
  const pid = project.snapshot.project_id;
  const [comparison, setComparison] = useState(() => {
    try {
      return localStorage.getItem(`katana:comparison:${pid}`) || project.manifest.contrasts[0]?.contrast_id || "";
    } catch {
      return project.manifest.contrasts[0]?.contrast_id || "";
    }
  });
  const [gene, setGene] = useState("");
  const [analysis, setAnalysis] = useState("");

  const presentation = useQuery({
    queryKey: ["report-presentation", pid],
    queryFn: () =>
      reportGet<{
        panels: ReportPanel[];
        complete_results: number;
        splicetrace_available?: boolean;
        motif_summary?: {contrast_id: string; sequence_report?: Record<string, unknown>}[];
      }>(pid, "report-presentation"),
  });

  const rx = route.match ? new RegExp(route.match, "i") : null;
  const candidates = artifacts.filter(a =>
    route.id === "reports" || route.id === "thresholds"
      ? true
      : rx
      ? rx.test(String(a.metadata.source_analysis_id))
      : route.chapters?.includes(Number(a.metadata.chapter)) || String(a.metadata.source_analysis_id).startsWith("project.")
  );
  const groups = [...new Set(candidates.map(a => String(a.metadata.source_analysis_id || "Other results")))];
  const active = groups.includes(analysis) ? analysis : groups[0];
  const tables = candidates.filter(
    a =>
      String(a.metadata.source_analysis_id || "Other results") === active &&
      (!a.contrast_ids?.length || a.contrast_ids.includes("all") || a.contrast_ids.includes(comparison))
  );
  const panels = (presentation.data?.panels || []).filter(
    p =>
      (!p.contrast_id || p.contrast_id === "all" || p.contrast_id === comparison) &&
      (route.id === "reports" ||
        (route.figure ? new RegExp(route.figure, "i").test(p.title) : route.chapters?.includes(Number(p.chapter))))
  );
  const tableFirst = ["deg", "das", "thresholds"].includes(route.id) || !!route.eventType;
  const projectPage = ["project", "comparisons", "samples"].includes(route.id);
  const inspectGene = (g: string, ev?: string) => {
    if (onOpenGeneCard) onOpenGeneCard(g, ev, comparison);
    else onNavigate("targets");
  };
  const isSingleCell = project.assay?.mode === "single-cell";

  return (
    <main className="workbench report-workbench">
      <div className="workbench-heading">
        <div>
          <span className="eyebrow">{project.assay?.label || "Transcriptomic Splicing Suite"} · {project.snapshot.genome_build}</span>
          <h2>{route.label}</h2>
        </div>
        <div style={{display: "flex", gap: 8, alignItems: "center"}}>
          <a className="secondary-action compact" href={reportUrl(pid, "REPORT.html")} target="_blank" rel="noreferrer">
            Standalone Study Report ↗
          </a>
          <button className="secondary-action compact" onClick={onSetup}>
            Analysis setup &amp; inputs
          </button>
        </div>
      </div>

      {route.id === "network" ? (
        <>
          <ReportNetwork projectId={pid} onGene={setGene} onInspectGeneCard={inspectGene} />
          {gene && (
            <section className="result-card mck-panel">
              <h3>Gene follow-up · {gene}</h3>
              <ReportExplorer
                key={gene}
                projectId={pid}
                artifacts={artifacts.filter(a => /^(differential.deg|differential.as)$/.test(String(a.metadata.source_analysis_id)))}
                initialSearch={gene}
                onAdd={onAdd}
              />
            </section>
          )}
        </>
      ) : (
        <>
          <div className="workflow-controls report-comparisons">
            <label>
              Active Comparison
              <select
                aria-label="Analysis comparison"
                value={comparison}
                onChange={e => {
                  setComparison(e.target.value);
                  try {
                    localStorage.setItem(`katana:comparison:${pid}`, e.target.value);
                  } catch {}
                }}
              >
                {project.manifest.contrasts.map(c => (
                  <option key={c.contrast_id} value={c.contrast_id}>
                    {c.label || c.contrast_id}
                  </option>
                ))}
              </select>
            </label>
          </div>

          {/* 01 · Executive Study Design & Sample Stratification */}
          {projectPage && (
            <section className="result-card report-design mck-panel">
              <div className="mck-panel-header">
                <div>
                  <span className="eyebrow">CHAPTER 01 · EXECUTIVE STUDY SCOPE &amp; REPLICATE ARCHITECTURE</span>
                  <h3>Study Stratification, Sample Ledger &amp; Interactive Lab Shortcuts</h3>
                </div>
              </div>
              <div className="mck-kpi-grid" style={{gridTemplateColumns: "repeat(4, minmax(0, 1fr))"}}>
                <div className="mck-kpi-box">
                  <small>{isSingleCell ? "QC-Passing Cells" : "Biological Samples"}</small>
                  <strong>{project.manifest.samples.length}</strong>
                  <span>{isSingleCell ? "Donor-stratified Smart-seq3" : "Bulk RNA-seq libraries"}</span>
                </div>
                <div className="mck-kpi-box">
                  <small>Declared Comparisons</small>
                  <strong>{project.manifest.contrasts.length}</strong>
                  <span>Exact numerator vs denominator</span>
                </div>
                <div className="mck-kpi-box">
                  <small>Reference Build</small>
                  <strong>{project.snapshot.genome_build}</strong>
                  <span>0-based half-open splice coords</span>
                </div>
                <div className="mck-kpi-box">
                  <small>Delivered Tables</small>
                  <strong>{presentation.data?.complete_results ?? artifacts.length}</strong>
                  <span>SHA256-verified Parquet / TSV</span>
                </div>
              </div>
              <div className="report-table-wrap" style={{maxHeight: 280}}>
                <table>
                  <thead>
                    <tr>
                      <th>{isSingleCell ? "Cell ID" : "Sample ID"}</th>
                      <th>Condition / Arm</th>
                      <th>Biological Replicate / Donor</th>
                    </tr>
                  </thead>
                  <tbody>
                    {project.manifest.samples.map(s => (
                      <tr key={s.sample_id}>
                        <td><strong>{s.sample_id}</strong></td>
                        <td>{s.condition}</td>
                        <td>{String(s.biological_sample_id || "Recorded")}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <div className="mck-filter-pills" style={{marginTop: 14}}>
                <button type="button" className="primary-action compact" onClick={() => onNavigate("patterns")}>
                  02 · Live PSI Z-Map, PCA &amp; UMAP Lab →
                </button>
                <button type="button" className="secondary-action compact" onClick={() => onNavigate("das")}>
                  03 · Differential PS-DAS &amp; DEG →
                </button>
                <button type="button" className="secondary-action compact" onClick={() => onNavigate("integrated")}>
                  04 · DEG × DAS &amp; Rescue Filter →
                </button>
                <button type="button" className="secondary-action compact" onClick={() => onNavigate("enrichment")}>
                  05 · Live GO &amp; KEGG ORA Lab →
                </button>
                <button type="button" className="secondary-action compact" onClick={() => onNavigate("network")}>
                  06 · Katana-SpliceNet →
                </button>
                <button type="button" className="secondary-action compact" onClick={() => onNavigate("motif")}>
                  07 · SpliceTrace &amp; PWM Scanner →
                </button>
                <button type="button" className="secondary-action compact" onClick={() => onNavigate("targets")}>
                  08 · 4-Track GeneCard &amp; Sashimi →
                </button>
              </div>
            </section>
          )}

          {/* 02 · Live PSI Landscape, Beyond-Jutils Z-Score Clustermap, PCA & Precision-Kernel UMAP Lab */}
          {route.id === "patterns" && (
            <PsiLandscapeLab
              projectId={pid}
              comparison={comparison}
              isSingleCell={isSingleCell}
              onInspectGeneCard={inspectGene}
            />
          )}

          {/* 02.1 · Single-Cell Dual Manifold Feature Overlay Explorer */}
          {(route.id === "sc-manifold" || (isSingleCell && route.id === "patterns")) && (
            <SingleCellManifoldExplorer projectId={pid} />
          )}

          {/* 04 · Cross-Layer DEG × DAS Quadrant & Multi-Contrast Rescue Filter */}
          {(route.id === "integrated" || route.id === "candidates") && (
            <QuadrantExplorer projectId={pid} comparison={comparison} onInspectGeneCard={inspectGene} />
          )}
          {(route.id === "integrated" || route.id === "cross") && project.manifest.contrasts.length > 1 && (
            <ContrastRescueExplorer projectId={pid} contrasts={project.manifest.contrasts} onInspectGeneCard={inspectGene} />
          )}

          {/* 05 · Live GO (BP/CC/MF) + KEGG Hypergeometric ORA Lab */}
          {route.id === "enrichment" && (
            <LiveOraWorkbench projectId={pid} comparison={comparison} onInspectGeneCard={g => inspectGene(g)} />
          )}

          {/* 07 · Katana-SpliceTrace Positional RBP, Live IUPAC / 170-PWM Scanner & Proteome Lab */}
          {(route.id === "motif" || route.parent === "motif") && (
            <SpliceTraceExplorer projectId={pid} comparison={comparison} onInspectGeneCard={inspectGene} />
          )}

          {/* 08.x · Sashimi / Scorecard quick preview */}
          {(route.id === "sashimi" || route.id === "candidates") && (
            <GeneCardDossier projectId={pid} contrasts={project.manifest.contrasts} entries={[]} initialComparison={comparison} />
          )}

          {/* Delivered High-Resolution Consulting Report Exhibits */}
          {!tableFirst && <ReportFigures key={comparison} projectId={pid} panels={panels} />}

          {/* Complete DuckDB-Backed Searchable & Filterable Evidence Tables */}
          {!!groups.length && (
            <>
              <div className="workflow-switcher" role="tablist" aria-label="Analysis methods">
                {groups.map(g => (
                  <button
                    role="tab"
                    aria-selected={g === active}
                    className={g === active ? "active" : ""}
                    key={g}
                    onClick={() => setAnalysis(g)}
                  >
                    <span>{g.replaceAll(".", " / ")}</span>
                  </button>
                ))}
              </div>
              {tableFirst ? (
                <ReportExplorer
                  key={`${active}:${comparison}`}
                  projectId={pid}
                  artifacts={tables}
                  eventType={route.eventType}
                  onAdd={onAdd}
                />
              ) : (
                <details className="result-card mck-panel" open={route.id === "reports"}>
                  <summary style={{cursor: "pointer", fontWeight: 700, color: "#051C2C"}}>
                    Complete Delivered Statistical Tables, Threshold Filtering &amp; Target Basket Selection
                  </summary>
                  <ReportExplorer key={`${active}:${comparison}`} projectId={pid} artifacts={tables} onAdd={onAdd} />
                </details>
              )}
            </>
          )}
        </>
      )}
      {presentation.error && <p role="alert">{presentation.error.message}</p>}
    </main>
  );
}
