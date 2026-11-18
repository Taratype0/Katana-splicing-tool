import {useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {api} from "../api";

const PALETTE = ["#051C2C", "#2251FF", "#0077C8", "#B83A4B", "#0D9488", "#D04A02", "#6366F1", "#64748B"];

export function SingleCellManifoldExplorer({projectId}: {projectId: string}) {
  const [colorBy, setColorBy] = useState<string>("activation_arm");
  const [feature, setFeature] = useState<string>("");
  const [cellTypeFilter, setCellTypeFilter] = useState<string>("all");

  const query = useQuery({
    queryKey: ["single-cell-manifold", projectId, feature],
    queryFn: () => api.singleCellManifold(projectId, feature),
  });

  if (query.isLoading) return <p>Loading single-cell Smart-seq3 dual manifold…</p>;
  if (!query.data?.available) return null;

  const data = query.data;
  const cells = data.cells.filter(c => cellTypeFilter === "all" || c.cell_type === cellTypeFilter);
  const donorsCount = new Set(cells.map(c => c.donor_id)).size;

  const categories = Array.from(new Set(cells.map(c => String(c[colorBy] ?? "Other"))));
  const catColors = Object.fromEntries(categories.map((cat, idx) => [cat, PALETTE[idx % PALETTE.length]]));

  const isContinuous = colorBy.startsWith("module_") || (colorBy === "feature" && Boolean(feature));
  const contVals = cells.map(c => Number(colorBy === "feature" ? c.feature_value : c[colorBy])).filter(v => Number.isFinite(v));
  const vmin = contVals.length ? Math.min(...contVals) : 0;
  const vmax = contVals.length ? Math.max(...contVals) : 1;

  function pointColor(c: Record<string, unknown>) {
    if (!isContinuous) return catColors[String(c[colorBy] ?? "Other")] || "#051C2C";
    const raw = Number(colorBy === "feature" ? c.feature_value : c[colorBy]);
    if (!Number.isFinite(raw)) return "#CBD5E1";
    const t = Math.max(0, Math.min(1, (raw - vmin) / Math.max(1e-6, vmax - vmin)));
    // Interpolate #93C5FD -> #2251FF -> #051C2C
    const r = Math.round(147 * (1 - t) + 5 * t);
    const g = Math.round(197 * (1 - t) + 28 * t);
    const b = Math.round(253 * (1 - t) + 44 * t);
    return `rgb(${r},${g},${b})`;
  }

  function renderProjection(title: string, xKey: "expr_PC1" | "psi_PC1", yKey: "expr_PC2" | "psi_PC2") {
    const xs = cells.map(c => Number(c[xKey]) || 0);
    const ys = cells.map(c => Number(c[yKey]) || 0);
    const minX = Math.min(...xs, -1), maxX = Math.max(...xs, 1);
    const minY = Math.min(...ys, -1), maxY = Math.max(...ys, 1);
    const W = 440, H = 300, pad = 34;
    const mapX = (v: number) => pad + ((v - minX) / Math.max(1e-6, maxX - minX)) * (W - 2 * pad);
    const mapY = (v: number) => H - pad - ((v - minY) / Math.max(1e-6, maxY - minY)) * (H - 2 * pad);

    return (
      <div className="mck-subcard">
        <h4>{title}</h4>
        <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label={title} style={{width: "100%", height: "auto"}}>
          <rect x={pad} y={pad} width={W - 2 * pad} height={H - 2 * pad} fill="#F8FAFC" stroke="#E2E8F0" />
          {cells.map(c => (
            <circle
              key={c.cell_id}
              cx={mapX(Number(c[xKey]) || 0)}
              cy={mapY(Number(c[yKey]) || 0)}
              r={4.8}
              fill={pointColor(c)}
              stroke="#FFFFFF"
              strokeWidth={0.6}
            >
              <title>{`${c.cell_id} · Donor ${c.donor_id} · ${c.activation_arm} · ${c.cell_type}`}</title>
            </circle>
          ))}
          <text x={W / 2} y={H - 8} textAnchor="middle" fontSize={10} fill="#051C2C">{xKey.replace("_", " ")}</text>
          <text x={12} y={H / 2} textAnchor="middle" fontSize={10} fill="#051C2C" transform={`rotate(-90 12 ${H / 2})`}>{yKey.replace("_", " ")}</text>
        </svg>
      </div>
    );
  }

  return (
    <section className="result-card mck-panel" aria-label="Single-cell Smart-seq3 Dual Manifold Explorer">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">SINGLE-CELL SMART-SEQ3 · DONOR-AWARE DUAL MANIFOLD</span>
          <h3>Interactive Expression vs Splicing PSI Manifold</h3>
        </div>
        <span className="mck-stat-pill">
          {donorsCount} Biological Donors · {cells.length} Smart-seq3 Cells (cells are not independent organisms)
        </span>
      </div>

      <div className="mck-controls-bar">
        <label>
          Color Manifold By
          <select aria-label="Color manifold by" value={colorBy} onChange={e => setColorBy(e.target.value)}>
            <option value="activation_arm">Activation State (Activated vs Resting)</option>
            <option value="cell_type">CD4+ T-Cell Subtype</option>
            <option value="donor_id">Biological Donor</option>
            {data.module_columns.map(mc => (
              <option key={mc} value={mc}>SC-SpliceNet {mc.replace("module_", "")} Activity</option>
            ))}
            <option value="feature">Selected Gene Expression / Event PSI</option>
          </select>
        </label>

        {colorBy === "feature" && (
          <label>
            Gene / Event Feature
            <select aria-label="Manifold feature" value={feature} onChange={e => setFeature(e.target.value)}>
              <option value="">Choose RBP gene or splicing event…</option>
              <optgroup label="Host RBP Expression">
                {data.expr_genes.slice(0, 40).map(g => <option key={g} value={g}>{g}</option>)}
              </optgroup>
              <optgroup label="Single-Cell Splicing Event PSI">
                {data.psi_events.slice(0, 40).map(ev => <option key={ev} value={ev}>{ev}</option>)}
              </optgroup>
            </select>
          </label>
        )}

        <label>
          Filter Cell Subtype
          <select aria-label="Filter cell subtype" value={cellTypeFilter} onChange={e => setCellTypeFilter(e.target.value)}>
            <option value="all">All CD4+ subtypes ({data.cells.length})</option>
            {data.cell_types.map(ct => <option key={ct} value={ct}>{ct}</option>)}
          </select>
        </label>
      </div>

      {!isContinuous && (
        <div className="mck-legend-row">
          {categories.map(cat => (
            <span key={cat} className="mck-legend-item">
              <i style={{background: catColors[cat]}} /> {cat}
            </span>
          ))}
        </div>
      )}

      <div className="mck-split-grid">
        {renderProjection("A. Host RBP / Gene Expression Manifold", "expr_PC1", "expr_PC2")}
        {renderProjection("B. Single-Cell Splicing PSI Manifold", "psi_PC1", "psi_PC2")}
      </div>
    </section>
  );
}
