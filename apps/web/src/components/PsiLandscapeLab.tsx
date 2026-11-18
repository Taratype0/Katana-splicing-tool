import {useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {api, type LivePsiLandscapePayload} from "../api";

const PALETTE = ["#051C2C", "#2251FF", "#C83E4D", "#0085A1", "#D97706", "#4F46E5", "#0F766E", "#475569"];
// Distinct jewel/slate palette for Splicing Subtypes — never shares hex colors with Sample Group PALETTE!
const EVENT_COLORS: Record<string, string> = {
  SE: "#64748B",
  RI: "#8B5CF6",
  A3SS: "#10B981",
  A5SS: "#F59E0B",
  MXE: "#EC4899",
};

function zScoreColor(z: number, clip: number): string {
  const frac = Math.max(-1, Math.min(1, z / Math.max(0.5, clip)));
  if (frac >= 0) {
    const r = Math.round(255 - frac * (255 - 184));
    const g = Math.round(255 - frac * (255 - 58));
    const b = Math.round(255 - frac * (255 - 75));
    return `rgb(${r},${g},${b})`;
  }
  const a = -frac;
  const r = Math.round(255 - a * 255);
  const g = Math.round(255 - a * (255 - 119));
  const b = Math.round(255 - a * (255 - 200));
  return `rgb(${r},${g},${b})`;
}

function ManifoldScatterSvg({
  title,
  subtitle,
  xLabel,
  yLabel,
  points,
  colorMode,
}: {
  title: string;
  subtitle: string;
  xLabel: string;
  yLabel: string;
  points: Array<{id: string; x: number; y: number; group: string; divergence?: number; meta: string}>;
  colorMode: "categorical" | "divergence";
}) {
  if (!points.length) return null;
  const W = 490, H = 310, padL = 50, padR = 20, padT = 20, padB = 42;
  const plotW = W - padL - padR, plotH = H - padT - padB;
  const xs = points.map(p => p.x);
  const ys = points.map(p => p.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  const spanX = Math.max(1e-4, maxX - minX) * 1.16;
  const spanY = Math.max(1e-4, maxY - minY) * 1.16;
  const midX = (minX + maxX) / 2, midY = (minY + maxY) / 2;
  const mapX = (v: number) => padL + ((v - (midX - spanX / 2)) / spanX) * plotW;
  const mapY = (v: number) => padT + (1 - (v - (midY - spanY / 2)) / spanY) * plotH;

  const groups = Array.from(new Set(points.map(p => p.group)));
  const groupCol = Object.fromEntries(groups.map((g, idx) => [g, PALETTE[idx % PALETTE.length]]));

  return (
    <div className="mck-subcard">
      <div className="mck-track-head">
        <strong>{title}</strong>
        <span>{subtitle}</span>
      </div>
      <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label={title} style={{width: "100%", height: "auto"}}>
        <rect x={padL} y={padT} width={plotW} height={plotH} fill="#F8FAFC" stroke="#E2E8F0" />
        <line x1={padL} y1={padT + plotH / 2} x2={padL + plotW} y2={padT + plotH / 2} stroke="#E2E8F0" strokeDasharray="3,3" />
        <line x1={padL + plotW / 2} y1={padT} x2={padL + plotW / 2} y2={padT + plotH} stroke="#E2E8F0" strokeDasharray="3,3" />

        {points.map((p, idx) => {
          const cx = mapX(p.x), cy = mapY(p.y);
          let fill = groupCol[p.group] || "#051C2C";
          if (colorMode === "divergence") {
            const d = Math.max(0, Math.min(1, p.divergence ?? 0));
            fill = d > 0.65 ? "#B83A4B" : d > 0.35 ? "#2251FF" : "#94A3B8";
          }
          return (
            <g key={p.id}>
              <circle
                cx={cx}
                cy={cy}
                r={points.length <= 36 ? 6.2 : 4.2}
                fill={fill}
                fillOpacity={0.88}
                stroke="#FFFFFF"
                strokeWidth={0.9}
              >
                <title>{p.meta}</title>
              </circle>
              {points.length <= 24 && (
                <text x={cx + 7} y={cy + (idx % 2 === 0 ? -4 : 9)} fontSize={8.2} fontWeight={600} fill="#051C2C">
                  {p.id}
                </text>
              )}
            </g>
          );
        })}

        <text x={padL + plotW / 2} y={H - 10} textAnchor="middle" fontSize={9.5} fontWeight={600} fill="#051C2C">
          {xLabel}
        </text>
        <text x={14} y={padT + plotH / 2} textAnchor="middle" fontSize={9.5} fontWeight={600} fill="#051C2C" transform={`rotate(-90 14 ${padT + plotH / 2})`}>
          {yLabel}
        </text>
      </svg>

      {colorMode === "categorical" ? (
        <div className="mck-legend-inline">
          <strong style={{fontSize: 11, color: "#051C2C", marginRight: 4}}>Sample Group:</strong>
          {groups.slice(0, 8).map(g => (
            <span key={g}>
              <i style={{background: groupCol[g]}} />
              {g}
            </span>
          ))}
        </div>
      ) : (
        <div className="mck-legend-inline">
          <span><i style={{background: "#94A3B8"}} />Conserved (Δ_kNN &lt; 0.35)</span>
          <span><i style={{background: "#2251FF"}} />Moderate divergence</span>
          <span><i style={{background: "#B83A4B"}} />Splicing-divergent state (Δ_kNN &gt; 0.65)</span>
        </div>
      )}
    </div>
  );
}

export function PsiLandscapeLab({
  projectId,
  comparison,
  isSingleCell,
  onInspectGeneCard,
}: {
  projectId: string;
  comparison: string;
  isSingleCell?: boolean;
  onInspectGeneCard?: (gene: string, eventId?: string) => void;
}) {
  const [mode, setMode] = useState<"unsupervised" | "supervised">("unsupervised");
  const [manifoldTopN, setManifoldTopN] = useState(2000);
  const [topN, setTopN] = useState(150);
  const [priorKappa, setPriorKappa] = useState(10);
  const [metric, setMetric] = useState<"euclidean" | "cosine" | "correlation">("euclidean");
  const [eventTypeFilter, setEventTypeFilter] = useState<string>("ALL");
  const [umapNeighbors, setUmapNeighbors] = useState(isSingleCell ? 14 : 8);
  const [colorBy, setColorBy] = useState<"condition" | "cell_type" | "donor" | "divergence">(isSingleCell ? "cell_type" : "condition");

  const [committed, setCommitted] = useState({
    mode: "unsupervised" as "unsupervised" | "supervised",
    manifold_top_n: 2000,
    top_n: 150,
    prior_kappa: 10,
    metric: "euclidean" as "euclidean" | "cosine" | "correlation",
    event_types: [] as string[],
    umap_neighbors: isSingleCell ? 14 : 8,
  });

  const query = useQuery<LivePsiLandscapePayload>({
    queryKey: ["live-psi-landscape", projectId, comparison, committed],
    queryFn: () =>
      typeof api.computePsiLandscape === "function"
        ? api.computePsiLandscape(projectId, {comparison, ...committed})
        : Promise.resolve({available: false} as LivePsiLandscapePayload),
  });

  function applyRecompute() {
    setCommitted({
      mode,
      manifold_top_n: manifoldTopN,
      top_n: topN,
      prior_kappa: priorKappa,
      metric,
      event_types: eventTypeFilter === "ALL" ? [] : [eventTypeFilter],
      umap_neighbors: umapNeighbors,
    });
  }

  const data = query.data;
  if (typeof api.computePsiLandscape !== "function") return null;

  const groupGetter = (s: LivePsiLandscapePayload["samples"][number]) => {
    if (colorBy === "cell_type") return s.cell_type || s.condition || "Cohort";
    if (colorBy === "donor") return s.donor_id || s.condition || "Sample";
    return s.condition || "Sample";
  };

  const manifoldUsed = data?.manifold_events ?? data?.selected_events ?? 0;
  const heatmapUsed = data?.heatmap_events ?? data?.heatmap_rows?.length ?? 0;

  return (
    <section className="result-card mck-panel" aria-label="Live PSI Z-Score Clustermap, PCA & Precision-Kernel UMAP Lab">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">GLOBAL SPLICING MANIFOLD · DEPTH-STABILIZED PSI, PERMANOVA PCA &amp; FISHER-KERNEL UMAP</span>
          <h3>Interactive Splicing PCA, Event-PSI UMAP &amp; Z-Score Clustermap Lab</h3>
        </div>
        {data?.files && (
          <div className="mck-filter-pills">
            <a className="secondary-action compact" href={data.files.clustermap_png} target="_blank" rel="noreferrer">Clustermap PNG ↗</a>
            <a className="secondary-action compact" href={data.files.pca_png} target="_blank" rel="noreferrer">PCA PNG ↗</a>
            {data.files.dual_umap_png && <a className="secondary-action compact" href={data.files.dual_umap_png} target="_blank" rel="noreferrer">Dual UMAP PNG ↗</a>}
            <a className="secondary-action compact" href={data.files.coordinates_tsv} download>Coordinates TSV ↓</a>
          </div>
        )}
      </div>

      {/* Plain-Language Methodological Summary */}
      <div className="run-plan-panel" style={{marginBottom: 12, fontSize: 12.5, lineHeight: 1.55}}>
        <strong>Methodological Upgrades Over Standard PSI Heatmaps (e.g. Jutils):</strong>{" "}
        <strong>(1) Depth-Stabilized PSI:</strong> shrinks noisy low-coverage ratios (such as 1-read <code>0% / 100%</code> extremes) toward the event mean (κ={priorKappa} reads) without deleting incomplete rows;{" "}
        <strong>(2) Binomial-Noise Deconvolution:</strong> subtracts expected count-sampling variance before selecting the top <strong>{manifoldUsed.toLocaleString()}</strong> biologically variable events for global PCA/UMAP;{" "}
        <strong>(3) Exact Statistical Separation:</strong> tests group separation via permutation PERMANOVA R² and joint read-precision UMAP kernel.
      </div>

      {/* Live Parameter Control Bar */}
      <div className="mck-controls-bar">
        <label>
          Feature Scope
          <select aria-label="Landscape feature selection mode" value={mode} onChange={e => setMode(e.target.value as "unsupervised" | "supervised")}>
            <option value="unsupervised">Unsupervised Biological Variance (Zero Leakage)</option>
            <option value="supervised">Supervised PS-DAS Signature ({comparison || "Contrast"})</option>
          </select>
        </label>
        <label>
          PCA / UMAP Events
          <select aria-label="Manifold events count" value={manifoldTopN} onChange={e => setManifoldTopN(Number(e.target.value))}>
            <option value={500}>Top 500 events</option>
            <option value={1000}>Top 1,000 events</option>
            <option value={2000}>Top 2,000 events (standard)</option>
            <option value={5000}>Top 5,000 events</option>
            <option value={10000}>Up to 10,000 / All PS-DAS</option>
          </select>
        </label>
        <label>
          Heatmap Display Rows
          <select aria-label="Top splicing events count" value={topN} onChange={e => setTopN(Number(e.target.value))}>
            {[60, 100, 150, 250].map(n => <option key={n} value={n}>Top {n} rows</option>)}
          </select>
        </label>
        <label>
          Depth Prior (κ)
          <select aria-label="Empirical Bayes prior kappa" value={priorKappa} onChange={e => setPriorKappa(Number(e.target.value))}>
            <option value={5}>κ = 5 reads (light)</option>
            <option value={10}>κ = 10 reads (balanced)</option>
            <option value={20}>κ = 20 reads (conservative)</option>
          </select>
        </label>
        <label>
          Linkage Metric
          <select aria-label="Clustering distance metric" value={metric} onChange={e => setMetric(e.target.value as "euclidean" | "cosine" | "correlation")}>
            <option value="euclidean">Euclidean · Ward</option>
            <option value="cosine">Cosine · Average</option>
            <option value="correlation">Pearson Correlation</option>
          </select>
        </label>
        <label>
          Event Subtype
          <select aria-label="Filter splicing event subtype" value={eventTypeFilter} onChange={e => setEventTypeFilter(e.target.value)}>
            <option value="ALL">All Types (SE/RI/A3/A5/MXE)</option>
            {["SE", "RI", "A3SS", "A5SS", "MXE"].map(et => <option key={et} value={et}>{et} only</option>)}
          </select>
        </label>
        <label>
          UMAP Neighbors
          <input
            aria-label="UMAP neighbor count"
            type="number"
            min={3}
            max={40}
            value={umapNeighbors}
            onChange={e => setUmapNeighbors(Math.max(3, Math.min(40, Number(e.target.value) || 10)))}
            style={{width: 68}}
          />
        </label>
        <label>
          Color Overlay
          <select aria-label="Manifold color overlay" value={colorBy} onChange={e => setColorBy(e.target.value as typeof colorBy)}>
            <option value="condition">Condition / Arm</option>
            <option value="cell_type">Cell Type / Subtype</option>
            <option value="donor">Biological Donor / Replicate</option>
            {data?.dual_manifold_metrics && <option value="divergence">Modality Divergence (Δ_kNN)</option>}
          </select>
        </label>
        <button type="button" className="primary-action compact" onClick={applyRecompute} disabled={query.isFetching}>
          {query.isFetching ? "Computing…" : "Recompute Manifold & Z-Map"}
        </button>
      </div>

      {query.isLoading && <p className="mck-hint">Fitting depth-stabilized PSI Z-scores, Ward dendrogram, PERMANOVA PCA &amp; Precision-Kernel UMAP…</p>}
      {query.error && <p className="form-error" role="alert">{query.error instanceof Error ? query.error.message : String(query.error)}</p>}
      {data && !data.available && <p className="mck-compact-note">{data.reason || "PSI matrix not available for live refitting."}</p>}

      {data?.available && (
        <>
          <div className="mck-kpi-grid">
            <div className="mck-kpi-box">
              <small>PERMANOVA Group Separation R²</small>
              <strong>{(data.permanova_r2 * 100).toFixed(1)}%</strong>
              <span>Permutation P = {data.permanova_p.toPrecision(3)}</span>
            </div>
            <div className="mck-kpi-box">
              <small>PCA Variance (PC1 + PC2)</small>
              <strong>{((data.pca_variance_pct[0] || 0) + (data.pca_variance_pct[1] || 0)).toFixed(1)}%</strong>
              <span>PC1 {data.pca_variance_pct[0]}% · PC2 {data.pca_variance_pct[1]}%</span>
            </div>
            <div className="mck-kpi-box">
              <small>PCA / UMAP Manifold Events</small>
              <strong>{manifoldUsed.toLocaleString()} / {data.total_candidate_events.toLocaleString()}</strong>
              <span>{data.selection_mode.startsWith("supervised") ? `All matching PS-DAS · Heatmap top ${heatmapUsed}` : `Top ${manifoldUsed.toLocaleString()} bio-variable · Heatmap top ${heatmapUsed}`}</span>
            </div>
            {data.dual_manifold_metrics ? (
              <div className="mck-kpi-box">
                <small>Cell-Type Silhouette &amp; Donor LISI</small>
                <strong>PSI {data.dual_manifold_metrics.cell_type_silhouette_psi.toFixed(2)} · LISI {data.dual_manifold_metrics.median_donor_lisi.toFixed(2)}</strong>
                <span>Mean modality Δ_kNN = {data.dual_manifold_metrics.mean_divergence.toFixed(2)}</span>
              </div>
            ) : (
              <div className="mck-kpi-box">
                <small>Splicing UMAP Distance Kernel</small>
                <strong>Fisher Precision</strong>
                <span>{data.umap_method}</span>
              </div>
            )}
          </div>

          <div className="mck-split-grid">
            <ManifoldScatterSvg
              title={`A. Depth-Stabilized Splicing PCA (${manifoldUsed.toLocaleString()} events)`}
              subtitle={`PC1 (${data.pca_variance_pct[0]}%) vs PC2 (${data.pca_variance_pct[1]}%)`}
              xLabel={`PC1 (${data.pca_variance_pct[0]}% variance)`}
              yLabel={`PC2 (${data.pca_variance_pct[1]}% variance)`}
              colorMode={colorBy === "divergence" ? "divergence" : "categorical"}
              points={data.samples.map(s => ({
                id: s.sample_id,
                x: s.PC1,
                y: s.PC2,
                group: groupGetter(s),
                divergence: s.modality_divergence,
                meta: `${s.sample_id} | ${s.condition}${s.cell_type && s.cell_type !== s.condition ? ` · ${s.cell_type}` : ""} | Detected ${(s.observed_event_fraction * 100).toFixed(0)}%`,
              }))}
            />
            <ManifoldScatterSvg
              title={`B. Precision-Kernel Event-PSI UMAP (${manifoldUsed.toLocaleString()} events)`}
              subtitle={data.dual_manifold_metrics ? "Procrustes-aligned against Gene Expression UMAP" : "Pairwise Binomial Fisher-precision kernel"}
              xLabel="Precision Event-PSI UMAP 1"
              yLabel="Precision Event-PSI UMAP 2"
              colorMode={colorBy === "divergence" ? "divergence" : "categorical"}
              points={data.samples.map(s => ({
                id: s.sample_id,
                x: s.psi_UMAP1,
                y: s.psi_UMAP2,
                group: groupGetter(s),
                divergence: s.modality_divergence,
                meta: `${s.sample_id} | ${s.condition}${s.cell_type && s.cell_type !== s.condition ? ` · ${s.cell_type}` : ""}${s.donor_id ? ` | Donor ${s.donor_id}` : ""}${s.modality_divergence != null ? ` | Δ_kNN=${s.modality_divergence.toFixed(2)}` : ""}`,
              }))}
            />
          </div>

          <div className="mck-subcard" style={{marginTop: 14}}>
            <div className="mck-track-head">
              <strong>C. Ward-Ordered Depth-Stabilized PSI Z-Score Table (±{data.z_clip} SD · Top {Math.min(60, data.heatmap_rows.length)} of {heatmapUsed} Rows)</strong>
              <span>Cerulean (−2.5 SD skipping) ↔ White (mean) ↔ Burgundy (+2.5 SD inclusion) · Click GeneCard → for 4-Track molecular view</span>
            </div>
            <div className="report-table-wrap" style={{maxHeight: 340}}>
              <table className="mck-zmap-table">
                <thead>
                  <tr>
                    <th>Gene · Event</th>
                    <th>Subtype</th>
                    <th className="num">Bio-Variance</th>
                    <th>Ward-Clustered Sample / Cell Z-Score Strip ({data.ordered_samples.length} columns)</th>
                    {onInspectGeneCard && <th>Inspect</th>}
                  </tr>
                </thead>
                <tbody>
                  {data.heatmap_rows.slice(0, 60).map(row => (
                    <tr key={row.event_id}>
                      <td>
                        <strong>{row.gene_symbol}</strong>
                        <small style={{display: "block", color: "#64748B", fontSize: 10}}>{row.event_id.slice(0, 36)}</small>
                      </td>
                      <td>
                        <span className="mck-tag" style={{borderColor: EVENT_COLORS[row.event_type] || "#64748B", color: EVENT_COLORS[row.event_type] || "#475569"}}>
                          {row.event_type}
                        </span>
                      </td>
                      <td className="num">{row.excess_biovar.toFixed(4)}</td>
                      <td style={{minWidth: 280}}>
                        <div style={{display: "flex", height: 15, width: "100%", border: "1px solid #E2E8F0"}}>
                          {row.z_values.map((zv, ci) => (
                            <span
                              key={ci}
                              title={`${data.ordered_samples[ci]} (${data.ordered_conditions[ci]}): Z = ${zv >= 0 ? "+" : ""}${zv.toFixed(2)}, PSI = ${row.psi_values[ci]}`}
                              style={{flex: 1, background: zScoreColor(zv, data.z_clip)}}
                            />
                          ))}
                        </div>
                      </td>
                      {onInspectGeneCard && (
                        <td>
                          <button type="button" className="mck-link-btn" onClick={() => onInspectGeneCard(row.gene_symbol, row.event_id)}>
                            GeneCard →
                          </button>
                        </td>
                      )}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </>
      )}
    </section>
  );
}
