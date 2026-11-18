import {useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {reportGet, reportUrl, VegaFigure} from "./ReportFigures";

type Row = Record<string, unknown>;
type GraphSpec = {
  nodes: string[];
  node_kme?: Record<string, number>;
  node_signed_kme?: Record<string, number>;
  node_genes?: Record<string, string>;
  node_modules?: Record<string, string>;
  edges: number[][];
  total_nodes: number;
  total_edges: number;
  rule: string;
};
type NetworkData = {
  available?: boolean;
  reason?: string;
  v2_mode?: boolean;
  engine?: string;
  manifest: Row;
  modules: Row[];
  samples: Row[];
  scores: Row[];
  members: Row[];
  genes: Row[];
  associations: Row[];
  pathways: Row[];
  bridge_rbps?: Row[];
  event_usage: Record<string, (number | null)[]>;
  graphs: Record<string, GraphSpec>;
  module_files: Record<string, {genes: string; events: string}>;
  source_directory: string;
  downloads: string[];
};

const MODULE_PALETTE = ["#2251FF", "#0D9488", "#D97706", "#64748B", "#8B5CF6", "#C83E4D", "#0085A1", "#4F46E5"];
const MODULE_HALO = ["#EFF6FF", "#ECFDF5", "#FFFBEB", "#F8FAFC", "#F5F3FF", "#FFF1F2", "#ECFEFF", "#EEF2FF"];

function fmtP(val: unknown): string {
  if (val == null || val === "" || Number.isNaN(Number(val))) return "NA";
  const n = Number(val);
  if (n <= 0) return "< 1e-4";
  if (n < 0.001) return n.toExponential(2);
  return n.toFixed(4);
}

function DataTable({
  rows,
  columns,
  onRow,
  onGeneCard,
}: {
  rows: Row[];
  columns: string[];
  onRow?: (row: Row) => void;
  onGeneCard?: (gene: string, event?: string) => void;
}) {
  const [page, setPage] = useState(0);
  return (
    <>
      <div className="report-table-wrap">
        <table>
          <thead>
            <tr>
              {columns.map(c => <th key={c}>{c.replaceAll("_", " ")}</th>)}
              {onGeneCard && <th>4-Track Dossier</th>}
            </tr>
          </thead>
          <tbody>
            {rows.slice(page * 50, page * 50 + 50).map((r, i) => (
              <tr
                key={i}
                onClick={() => onRow?.(r)}
                tabIndex={onRow ? 0 : undefined}
                onKeyDown={e => e.key === "Enter" && onRow?.(r)}
              >
                {columns.map(c => (
                  <td key={c}>
                    {r[c] == null
                      ? "NA"
                      : typeof r[c] === "number"
                      ? c.includes("p_value") || c.includes("q_value") || c === "fdr" || c === "pvalue"
                        ? fmtP(r[c])
                        : Number(r[c]).toPrecision(4)
                      : String(r[c])}
                  </td>
                ))}
                {onGeneCard && (
                  <td>
                    <button
                      type="button"
                      className="mck-link-btn"
                      onClick={ev => {
                        ev.stopPropagation();
                        onGeneCard(String(r.gene_symbol || r.gene_id || ""), String(r.event_key || ""));
                      }}
                    >
                      GeneCard →
                    </button>
                  </td>
                )}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div className="target-actions">
        <button disabled={!page} onClick={() => setPage(page - 1)}>Previous</button>
        <span>{rows.length} rows</span>
        <button disabled={(page + 1) * 50 >= rows.length} onClick={() => setPage(page + 1)}>Next</button>
      </div>
    </>
  );
}

export function ReportNetwork({
  projectId,
  onGene,
  onInspectGeneCard,
}: {
  projectId: string;
  onGene?: (gene: string) => void;
  onInspectGeneCard?: (gene: string, eventId?: string) => void;
}) {
  const query = useQuery({
    queryKey: ["report-network", projectId],
    queryFn: () => reportGet<NetworkData>(projectId, "report-network"),
  });
  const [zoom, setZoom] = useState(1);
  const [minKme, setMinKme] = useState(0.0);
  const [viewScope, setViewScope] = useState<string>("ALL");
  const [event, setEvent] = useState("");
  const [search, setSearch] = useState("");
  const [section, setSection] = useState<"network" | "associations" | "genes" | "events" | "pathways">("network");

  if (query.isLoading) return <p>Loading co-splicing gene network &amp; module significance…</p>;
  if (query.error) return <p role="alert">{query.error.message}</p>;
  const d = query.data;
  if (!d) return null;

  if (d.available === false || d.modules.length === 0) {
    return (
      <section className="mck-card" aria-label="Co-splicing network availability">
        <div className="mck-kicker">{d.engine || "Katana-SpliceNet"}</div>
        <h2>Co-Splicing Network Module Not Executed for This Cohort</h2>
        <p className="mck-note" style={{ marginTop: 8 }}>
          {d.reason || "Unsupervised scale-free co-splicing module discovery requires multi-sample cohort variance and was not included in this report archive."}
        </p>
      </section>
    );
  }

  const hasAllGraph = Boolean(d.graphs["ALL"] && d.graphs["ALL"].nodes.length > 0);
  const activeScope = viewScope === "ALL" && !hasAllGraph ? String(d.modules[0]?.module || "M1") : viewScope;
  const isAllMode = activeScope === "ALL";

  const modColorMap: Record<string, string> = {};
  const modHaloMap: Record<string, string> = {};
  d.modules.forEach((m, idx) => {
    const mName = String(m.module);
    modColorMap[mName] = MODULE_PALETTE[idx % MODULE_PALETTE.length];
    modHaloMap[mName] = MODULE_HALO[idx % MODULE_HALO.length];
  });

  const focusedModule = d.modules.find(m => String(m.module) === activeScope) || d.modules[0];
  const focusedModName = String(focusedModule?.module || "M1");
  const graph = d.graphs[isAllMode ? "ALL" : focusedModName] || d.graphs[focusedModName];

  const activeNodeSet = new Set(
    (graph?.nodes || []).filter(id => !graph?.node_kme || (graph.node_kme[id] ?? 1) >= minKme)
  );
  const activeNodes = (graph?.nodes || []).filter(id => activeNodeSet.has(id));

  // Build spatial coordinates:
  // - In "ALL" mode: Multi-community sector layout (each module forms a labeled spatial cluster with radial outward gene names!)
  // - In single-module mode: Concentric 2-ring Hub-and-Spoke layout (Top 6 Core Hubs inner ring + outer Secondary Hubs, all labeled!)
  const W = isAllMode ? 1020 : 650;
  const H = isAllMode ? 670 : 470;
  const cxMain = W / 2;
  const cyMain = isAllMode ? 350 : H / 2 + 8;

  type NodePos = {x: number; y: number; lx: number; ly: number; anchor: "start" | "end" | "middle"; mod: string; isCore: boolean};
  const positions = new Map<string, NodePos>();
  const clusterHeaders: Array<{mod: string; cx: number; cy: number; rHalo: number; sigText: string; isSig: boolean}> = [];

  if (isAllMode) {
    const presentMods = Array.from(new Set(activeNodes.map(id => graph?.node_modules?.[id] || "M1")));
    const nMods = Math.max(1, presentMods.length);
    const R_cluster = nMods <= 2 ? 210 : nMods <= 4 ? 192 : 202;
    const rHalo = nMods <= 4 ? 102 : 84;

    presentMods.forEach((mName, mIdx) => {
      const thetaC = -Math.PI / 2 + (mIdx * 2 * Math.PI) / nMods;
      const ccx = cxMain + R_cluster * 1.20 * Math.cos(thetaC);
      const ccy = cyMain + R_cluster * 0.88 * Math.sin(thetaC);
      const mRec = d.modules.find(x => String(x.module) === mName);
      const bestFdr = Number(mRec?.omnibus_fdr ?? mRec?.best_fdr ?? 1);
      const bestP = Number(mRec?.best_p ?? mRec?.omnibus_p ?? 1);
      const isSig = bestFdr <= 0.05;
      const sigText = bestFdr <= 0.05
        ? `FDR=${fmtP(bestFdr)} ★`
        : bestP <= 0.055
        ? `P=${fmtP(bestP)} · FDR=${fmtP(bestFdr)}`
        : `FDR=${fmtP(bestFdr)} (NS)`;
      clusterHeaders.push({mod: mName, cx: ccx, cy: ccy, rHalo, sigText, isSig});

      const modNodes = activeNodes.filter(id => (graph?.node_modules?.[id] || "M1") === mName);
      modNodes.forEach((id, idx) => {
        const nHere = Math.max(1, modNodes.length);
        const ang = -Math.PI / 2 + (idx * 2 * Math.PI) / nHere;
        const rad = idx < 3 && nHere > 7 ? rHalo * 0.36 : rHalo * 0.78;
        const nx = ccx + rad * Math.cos(ang);
        const ny = ccy + rad * Math.sin(ang);
        const cosA = Math.cos(ang);
        const anchor: "start" | "end" | "middle" = Math.abs(cosA) < 0.18 ? "middle" : cosA > 0 ? "start" : "end";
        const lx = nx + (anchor === "start" ? 8 : anchor === "end" ? -8 : 0);
        const ly = ny + (anchor === "middle" ? (Math.sin(ang) > 0 ? 13 : -9) : 3.5);
        positions.set(id, {x: nx, y: ny, lx, ly, anchor, mod: mName, isCore: idx < 3});
      });
    });
  } else {
    const coreCount = Math.min(6, Math.max(3, Math.floor(activeNodes.length * 0.28)));
    const coreNodes = activeNodes.slice(0, coreCount);
    const outerNodes = activeNodes.slice(coreCount);
    coreNodes.forEach((id, i) => {
      const ang = -Math.PI / 2 + (i * 2 * Math.PI) / Math.max(1, coreNodes.length);
      const rCore = outerNodes.length > 0 ? 78 : 145;
      const nx = cxMain + rCore * Math.cos(ang);
      const ny = cyMain + rCore * Math.sin(ang);
      const cosA = Math.cos(ang);
      const anchor: "start" | "end" | "middle" = Math.abs(cosA) < 0.2 ? "middle" : cosA > 0 ? "start" : "end";
      positions.set(id, {
        x: nx,
        y: ny,
        lx: nx + (anchor === "start" ? 9 : anchor === "end" ? -9 : 0),
        ly: ny + (anchor === "middle" ? (Math.sin(ang) > 0 ? 14 : -10) : 4),
        anchor,
        mod: focusedModName,
        isCore: true,
      });
    });
    outerNodes.forEach((id, i) => {
      const ang = -Math.PI / 2 + ((i + 0.35) * 2 * Math.PI) / Math.max(1, outerNodes.length);
      const rOutX = 212;
      const rOutY = 172;
      const nx = cxMain + rOutX * Math.cos(ang);
      const ny = cyMain + rOutY * Math.sin(ang);
      const cosA = Math.cos(ang);
      const anchor: "start" | "end" | "middle" = Math.abs(cosA) < 0.16 ? "middle" : cosA > 0 ? "start" : "end";
      positions.set(id, {
        x: nx,
        y: ny,
        lx: nx + (anchor === "start" ? 8 : anchor === "end" ? -8 : 0),
        ly: ny + (anchor === "middle" ? (Math.sin(ang) > 0 ? 13 : -9) : 3.5),
        anchor,
        mod: focusedModName,
        isCore: false,
      });
    });
  }

  const membersInScope = isAllMode
    ? d.members.filter(r => Math.abs(Number(r.kME ?? 1)) >= minKme)
    : d.members.filter(r => r.module === focusedModName && Math.abs(Number(r.kME ?? 1)) >= minKme);
  const selected = membersInScope.find(r => r.event_key === event) || d.members.find(r => r.event_key === activeNodes[0]) || membersInScope[0];
  const selectedKey = String(selected?.event_key || "");
  const eigengeneMod = String(selected?.module || focusedModName);

  const scores = d.scores.map(r => ({
    sample: String(r.sample_id),
    condition: d.samples.find(s => s.sample_id === r.sample_id)?.condition,
    value: r[eigengeneMod],
  }));
  const usage = d.event_usage[selectedKey]
    ? d.samples.map((s, i) => ({sample: s.sample_id, condition: s.condition, value: d.event_usage[selectedKey]?.[i] ?? null}))
    : [];

  const chart = (values: Row[], title: string) => ({
    $schema: "https://vega.github.io/schema/vega-lite/v5.json",
    width: "container",
    height: 220,
    config: {font: "IBM Plex Sans, Arial, sans-serif", view: {stroke: null}},
    data: {values},
    mark: {type: "point", filled: true, size: 88},
    encoding: {
      x: {field: "sample", type: "nominal", axis: {title: "Biological sample / cell", labelAngle: -45, labels: values.length <= 36}},
      y: {field: "value", type: "quantitative", axis: {title}, scale: title.includes("PSI") ? {domain: [0, 1]} : {zero: false}},
      color: {field: "condition", type: "nominal", scale: {range: ["#051C2C", "#2251FF", "#C83E4D", "#0085A1", "#D97706"]}},
      tooltip: [{field: "sample"}, {field: "condition"}, {field: "value"}],
    },
  });

  const file = (n: string) => reportUrl(projectId, `${d.source_directory}/report_tables/${n}`);
  const modBridges = (d.bridge_rbps || []).filter(b =>
    isAllMode ? true : String(b.module) === focusedModName || `M${String(b.module)}` === focusedModName
  );
  const modAssocs = d.associations.filter(a => (isAllMode ? true : String(a.module) === eigengeneMod));

  return (
    <div className="native-network">
      {/* Header + Module Statistical Significance Scorecard Grid */}
      <section className="result-card mck-panel">
        <div className="mck-panel-header">
          <div>
            <span className="eyebrow">{d.engine || "KATANA-SPLICENET · CO-SPLICING GENE NETWORK"}</span>
            <h2>Co-Splicing Gene Network &amp; Module Statistical Significance</h2>
          </div>
          <span className="mck-stat-pill">
            {d.modules.length} Co-Splicing Modules · {d.members.length.toLocaleString()} Assigned Splicing Events
          </span>
        </div>
        <p className="mck-hint" style={{marginBottom: 12}}>
          <strong>What is Katana-SpliceNet?</strong> Unsupervised co-splicing network that groups splicing events whose exon inclusion (PSI) co-varies across all biological samples/cells into modules (<code>M1..M{d.modules.length}</code>), then tests each module&apos;s sample eigengene across conditions via ANOVA / Welch t-tests and BH FDR.
        </p>

        {/* Interactive Module Cards with Explicit Significance (FDR / P / Delta-ME) & Top Pathway! */}
        <div style={{display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(235px, 1fr))", gap: 10, marginBottom: 12}}>
          {hasAllGraph && (
            <button
              type="button"
              onClick={() => {setViewScope("ALL"); setSection("network");}}
              style={{
                textAlign: "left",
                padding: "11px 13px",
                borderRadius: 6,
                border: isAllMode ? "2px solid #2251FF" : "1px solid #CBD5E1",
                background: isAllMode ? "#EFF6FF" : "#FFFFFF",
                cursor: "pointer",
              }}
            >
              <div style={{display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 4}}>
                <strong style={{color: "#051C2C", fontSize: 13.5}}>Whole Network (All Modules)</strong>
                <span className="mck-tag mck-tag-blue">Default Overview</span>
              </div>
              <div style={{fontSize: 12, color: "#334155", marginBottom: 3}}>
                {d.modules.length} modules · {d.graphs["ALL"]?.nodes.length || 0} hub genes · {d.graphs["ALL"]?.total_edges || 0} PSI edges
              </div>
              <small style={{color: "#64748B", display: "block"}}>
                Simultaneous view of all co-splicing communities, hub gene symbols &amp; inter-module correlations
              </small>
            </button>
          )}

          {d.modules.map(m => {
            const mName = String(m.module);
            const isSel = activeScope === mName;
            const omFdr = m.omnibus_fdr != null ? Number(m.omnibus_fdr) : null;
            const bestFdr = m.best_fdr != null ? Number(m.best_fdr) : null;
            const bestP = m.best_p != null ? Number(m.best_p) : null;
            const fdrShow = omFdr ?? bestFdr;
            const sigClass = fdrShow != null && fdrShow <= 0.05 ? "sig" : bestP != null && bestP <= 0.05 ? "nom" : "ns";
            return (
              <button
                key={mName}
                type="button"
                aria-label={`Open module ${mName}`}
                onClick={() => {setViewScope(mName); setEvent(""); setSection("network");}}
                style={{
                  textAlign: "left",
                  padding: "11px 13px",
                  borderRadius: 6,
                  border: isSel ? `2px solid ${modColorMap[mName] || "#2251FF"}` : "1px solid #E2E8F0",
                  borderLeft: `5px solid ${modColorMap[mName] || "#051C2C"}`,
                  background: isSel ? "#F8FAFC" : "#FFFFFF",
                  cursor: "pointer",
                }}
              >
                <div style={{display: "flex", justifyContent: "space-between", alignItems: "center", gap: 6, marginBottom: 4}}>
                  <strong style={{color: "#051C2C", fontSize: 13.5}}>
                    {mName} · {String(m.unique_annotated_genes)} genes ({String(m.events)} events)
                  </strong>
                  <span
                    style={{
                      fontSize: 10.5,
                      fontWeight: 700,
                      padding: "2px 6px",
                      borderRadius: 4,
                      background: sigClass === "sig" ? "#DBEAFE" : sigClass === "nom" ? "#FEF3C7" : "#F1F5F9",
                      color: sigClass === "sig" ? "#1E40AF" : sigClass === "nom" ? "#92400E" : "#64748B",
                    }}
                  >
                    {sigClass === "sig" ? `FDR = ${fmtP(fdrShow)} ★` : sigClass === "nom" ? `P = ${fmtP(bestP)}` : `FDR = ${fmtP(fdrShow)} (NS)`}
                  </span>
                </div>
                {m.best_contrast ? (
                  <div style={{fontSize: 11.5, color: "#1E293B", marginBottom: 3}}>
                    <strong>Peak shift:</strong> {String(m.best_contrast)} (ΔME = {Number(m.best_effect ?? 0) >= 0 ? "+" : ""}{Number(m.best_effect ?? 0).toFixed(2)}, P = {fmtP(m.best_p)})
                  </div>
                ) : null}
                {m.top_pathway ? (
                  <div style={{fontSize: 11.2, color: "#0F766E", marginBottom: 3}}>
                    <strong>Pathway:</strong> {String(m.top_pathway).slice(0, 44)} (FDR {fmtP(m.top_pathway_fdr)})
                  </div>
                ) : null}
                <small style={{color: "#64748B", display: "block"}}>
                  <strong>Hubs:</strong> {String(m.top_hub_genes || m.top_hub_gene || "—")}
                </small>
              </button>
            );
          })}
        </div>

        {/* Network Controls + Significant Enriched RBP Regulators (only real FE >= 1.15, P <= 0.05!) */}
        <div className="mck-controls-bar">
          <label>
            Network View
            <select
              aria-label="CoSplice module"
              value={activeScope}
              onChange={e => {setViewScope(e.target.value); setEvent(""); setSearch("");}}
            >
              {hasAllGraph && <option value="ALL">All Modules · Whole Co-Splicing Network ({d.modules.length} modules)</option>}
              {d.modules.map(m => (
                <option key={String(m.module)} value={String(m.module)}>
                  {String(m.module)} · {String(m.unique_annotated_genes)} genes · {m.omnibus_fdr != null ? `ANOVA FDR=${fmtP(m.omnibus_fdr)}` : `${String(m.events)} events`}
                </option>
              ))}
            </select>
          </label>
          <label>
            Hub |kME| ≥ <strong>{minKme.toFixed(2)}</strong>
            <input type="range" aria-label="Minimum absolute kME" min="0" max="0.85" step="0.05" value={minKme} onChange={e => setMinKme(Number(e.target.value))} />
          </label>
          <label>
            Canvas Zoom
            <input type="range" aria-label="Network zoom" min="0.7" max="2.2" step="0.1" value={zoom} onChange={e => setZoom(Number(e.target.value))} />
          </label>
        </div>

        {!!modBridges.length && (
          <div className="mck-bridge-strip" style={{marginTop: 8}}>
            <strong>Enriched Upstream Cis-Motif RBP Regulators ({isAllMode ? "Study-Wide" : focusedModName}): </strong>
            {modBridges.slice(0, 8).map((b, idx) => (
              <span key={idx} className="mck-tag mck-tag-blue" style={{marginRight: 6}}>
                {isAllMode ? `${String(b.module)} ← ` : ""}{String(b.rbp || b.RBP)} ({b.fold_enrichment ? `${Number(b.fold_enrichment).toFixed(2)}×, ` : ""}P = {fmtP(b.p_value ?? b.fwer_p)})
              </span>
            ))}
          </div>
        )}

        {/* Clean Human-Readable Consulting Tabs */}
        <div className="artifact-tabs" role="tablist" aria-label="CoSplice views" style={{marginTop: 12}}>
          {([
            ["network", "Interactive Gene Network"],
            ["associations", "Contrast Significance (ΔME & FDR)"],
            ["genes", "Module Hub Genes (kME)"],
            ["events", "Assigned Splicing Events"],
            ["pathways", "GO & Pathway Enrichment"],
          ] as const).map(([id, label]) => (
            <button
              role="tab"
              aria-selected={section === id}
              className={section === id ? "active" : ""}
              key={id}
              onClick={() => setSection(id)}
            >
              {label}
            </button>
          ))}
        </div>
      </section>

      {section === "network" && (
        <>
          <section className="result-card mck-panel">
            <div className="mck-panel-header">
              <h3>
                {isAllMode
                  ? `Whole-Study Co-Splicing Gene Network (${activeNodes.length} Labeled Hub Genes Across ${d.modules.length} Modules)`
                  : `Module ${focusedModName} Hub-and-Spoke Gene Network (${activeNodes.length} Unique Hub Genes)`}
              </h3>
              <div className="mck-legend-inline">
                <span><i style={{background: "#2251FF", height: 3, width: 14, borderRadius: 0}} />Positive co-splicing (r &gt; 0)</span>
                <span><i style={{background: "#C83E4D", height: 3, width: 14, borderRadius: 0}} />Inverse / mutually-exclusive (r &lt; 0)</span>
              </div>
            </div>

            <div className={isAllMode ? undefined : "native-network-grid"}>
              <div style={{border: "1px solid #E2E8F0", borderRadius: 6, background: "#FFFFFF", padding: 6}}>
                <svg
                  viewBox={`${cxMain - W / (2 * zoom)} ${cyMain - H / (2 * zoom)} ${W / zoom} ${H / zoom}`}
                  role="img"
                  aria-label={`Interactive ${activeScope} co-splicing gene network`}
                  style={{width: "100%", height: "auto", maxHeight: isAllMode ? 590 : 470}}
                >
                  {/* Community Halos & Significance Headers in Whole-Network Mode */}
                  {isAllMode &&
                    clusterHeaders.map(ch => (
                      <g key={ch.mod} style={{cursor: "pointer"}} onClick={() => setViewScope(ch.mod)}>
                        <circle
                          cx={ch.cx}
                          cy={ch.cy}
                          r={ch.rHalo + 18}
                          fill={modHaloMap[ch.mod] || "#F8FAFC"}
                          stroke={modColorMap[ch.mod] || "#CBD5E1"}
                          strokeOpacity={0.35}
                          strokeDasharray="4,3"
                        />
                        <text
                          x={ch.cx}
                          y={ch.cy - ch.rHalo - 24}
                          textAnchor="middle"
                          fontSize={12}
                          fontWeight={700}
                          fill={modColorMap[ch.mod] || "#051C2C"}
                        >
                          {ch.mod} ({ch.sigText})
                        </text>
                      </g>
                    ))}

                  {/* Real Pairwise Sample-PSI Correlation Edges */}
                  {(graph?.edges || []).map((ed, i) => {
                    const [a, b, w, signedR] = ed;
                    const idA = graph.nodes[a], idB = graph.nodes[b];
                    if (!activeNodeSet.has(idA) || !activeNodeSet.has(idB)) return null;
                    const p = positions.get(idA), q = positions.get(idB);
                    if (!p || !q) return null;
                    const isNeg = (signedR ?? w) < 0;
                    return (
                      <line
                        key={i}
                        x1={p.x}
                        y1={p.y}
                        x2={q.x}
                        y2={q.y}
                        stroke={isNeg ? "#C83E4D" : "#3B82F6"}
                        strokeDasharray={isNeg ? "4,3" : undefined}
                        strokeOpacity={0.22 + Math.min(0.65, w * 0.58)}
                        strokeWidth={0.7 + w * 1.8}
                      />
                    );
                  })}

                  {/* Labeled Gene Hub Nodes (100% of Displayed Hub Genes Have Visible Gene Names!) */}
                  {activeNodes.map(id => {
                    const p = positions.get(id);
                    if (!p) return null;
                    const gLabel = String(graph?.node_genes?.[id] || id.slice(0, 14));
                    const kVal = graph?.node_kme?.[id] ?? 0.6;
                    const isSel = id === selectedKey;
                    const nodeCol = isSel ? "#C83E4D" : modColorMap[p.mod] || "#051C2C";
                    const rad = isSel ? 8.8 : p.isCore ? 6.8 : 5.2;
                    return (
                      <g
                        key={id}
                        role="button"
                        aria-label={`Select ${gLabel}`}
                        tabIndex={0}
                        style={{cursor: "pointer"}}
                        onClick={() => setEvent(id)}
                        onKeyDown={e => e.key === "Enter" && setEvent(id)}
                      >
                        <circle cx={p.x} cy={p.y} r={rad} fill={nodeCol} stroke="#FFFFFF" strokeWidth={1.3}>
                          <title>{`${gLabel} (${p.mod}) · |kME|=${kVal.toFixed(3)} · ${id}`}</title>
                        </circle>
                        <text
                          x={p.lx}
                          y={p.ly}
                          textAnchor={p.anchor}
                          fontSize={isSel ? 11.2 : p.isCore ? 10.2 : 9.3}
                          fontWeight={isSel || p.isCore ? 700 : 600}
                          fill={isSel ? "#C83E4D" : "#051C2C"}
                          style={{paintOrder: "stroke", stroke: "#FFFFFF", strokeWidth: 2.6, strokeLinejoin: "round"}}
                        >
                          {isAllMode ? gLabel : `${gLabel} (${kVal.toFixed(2)})`}
                        </text>
                      </g>
                    );
                  })}
                </svg>
                <p className="mck-hint" style={{margin: "6px 4px 2px"}}>{graph?.rule}</p>
              </div>

              {/* Side-by-side Eigengene Significance & Sample Trajectory */}
              <div style={{marginTop: isAllMode ? 14 : 0}}>
                <div className="mck-subcard" style={{marginBottom: 12}}>
                  <div className="mck-track-head">
                    <strong>Module {eigengeneMod} Statistical Significance Across Comparisons</strong>
                    <span>Omnibus ANOVA &amp; Pairwise Welch Eigengene Tests</span>
                  </div>
                  <div className="report-table-wrap" style={{maxHeight: 185}}>
                    <table>
                      <thead>
                        <tr>
                          <th>Module</th>
                          <th>Comparison</th>
                          <th className="num">ΔEigengene</th>
                          <th className="num">P-value</th>
                          <th className="num">BH FDR</th>
                        </tr>
                      </thead>
                      <tbody>
                        {modAssocs.map((a, idx) => {
                          const qv = Number(a.q_value ?? 1);
                          return (
                            <tr key={idx} style={{fontWeight: qv <= 0.05 ? 700 : 400}}>
                              <td>{String(a.module)}</td>
                              <td>{String(a.contrast_id)}</td>
                              <td className="num">{a.effect != null ? Number(a.effect).toFixed(3) : "—"}</td>
                              <td className="num">{fmtP(a.p_value)}</td>
                              <td className="num" style={{color: qv <= 0.05 ? "#2251FF" : "#475569"}}>
                                {fmtP(a.q_value)}{qv <= 0.05 ? " ★" : ""}
                              </td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </div>

                <div className="mck-subcard">
                  <div className="mck-track-head">
                    <strong>Module {eigengeneMod} Sample Eigengene Activity</strong>
                    <span>Per-sample co-splicing program score</span>
                  </div>
                  <VegaFigure spec={chart(scores, `${eigengeneMod} eigengene score`)} title="Module score by sample" />
                </div>
              </div>
            </div>
          </section>

          {/* Selected Hub Gene / Event Detail + Direct 4-Track GeneCard Link */}
          <section className="result-card mck-panel">
            <h3>
              Selected Hub Gene · {String(selected?.gene_symbol || selected?.gene_id || "")}{" "}
              <small style={{fontWeight: 500, color: "#64748B"}}>
                ({String(selected?.module || focusedModName)} · {String(selected?.event_type || "SE")} · kME = {Number(selected?.kME ?? 0).toFixed(3)})
              </small>
            </h3>
            <div className="mck-controls-bar">
              <select aria-label="CoSplice event" value={selectedKey} onChange={e => setEvent(e.target.value)}>
                {membersInScope.slice(0, 250).map(m => (
                  <option key={String(m.event_key)} value={String(m.event_key)}>
                    [{String(m.module)}] {String(m.gene_symbol || m.gene_id)} · {String(m.event_type)} · kME={Number(m.kME ?? 0).toFixed(2)} · {String(m.event_key).slice(0, 40)}
                  </option>
                ))}
              </select>
              <button type="button" className="secondary-action compact" onClick={() => onGene?.(String(selected?.gene_id || selected?.gene_symbol || ""))}>
                Follow in PS-DAS / DEG
              </button>
              {onInspectGeneCard && selected && (
                <button
                  type="button"
                  className="primary-action compact"
                  onClick={() => onInspectGeneCard(String(selected.gene_symbol || selected.gene_id || ""), String(selected.event_key || ""))}
                >
                  Open 4-Track GeneCard &amp; Sashimi →
                </button>
              )}
            </div>
            {!!usage.length && <VegaFigure spec={chart(usage, "Measured Sample PSI")} title="Selected event PSI by sample" />}
          </section>
        </>
      )}

      {section !== "network" && (
        <section className="result-card mck-panel">
          <h3>
            {section === "genes"
              ? `Module Hub Genes (${isAllMode ? "All Modules" : focusedModName})`
              : section === "events"
              ? `Assigned Splicing Events (${isAllMode ? "All Modules" : focusedModName})`
              : section === "associations"
              ? "Module Contrast Statistical Significance (ANOVA & Pairwise Welch Tests)"
              : `Module GO & Pathway Enrichment (${isAllMode ? "All Modules" : focusedModName})`}
          </h3>
          <input
            aria-label="Search module members"
            value={search}
            onChange={e => setSearch(e.target.value)}
            placeholder="Filter by gene symbol, module, comparison, or pathway…"
            style={{marginBottom: 10, padding: "6px 10px", width: "100%", maxWidth: 420}}
          />
          <DataTable
            key={`${activeScope}:${section}:${search}:${minKme}`}
            rows={(section === "genes"
              ? d.genes
              : section === "events"
              ? membersInScope
              : section === "associations"
              ? d.associations
              : d.pathways
            ).filter(
              r => (isAllMode || r.module === focusedModName) && JSON.stringify(r).toLowerCase().includes(search.toLowerCase())
            )}
            columns={
              section === "genes"
                ? ["module", "gene_symbol", "gene_id", "events", "max_abs_kME", "within_module_weight"]
                : section === "events"
                ? ["module", "gene_symbol", "event_type", "event_key", "kME", "within_module_weight"]
                : section === "associations"
                ? ["module", "kind", "contrast_id", "effect", "statistic", "p_value", "q_value"]
                : ["module", "clean_term", "pvalue", "fdr"]
            }
            onRow={
              section === "events"
                ? r => {
                    setEvent(String(r.event_key));
                    setSection("network");
                  }
                : undefined
            }
            onGeneCard={section === "genes" || section === "events" ? onInspectGeneCard : undefined}
          />
          {section === "pathways" && !d.pathways.length && <p className="mck-hint">No pathway rows were recorded for this module set.</p>}
        </section>
      )}

      {!!d.downloads.length && (
        <section className="result-card mck-panel">
          <h3>Complete Network Files &amp; Provenance</h3>
          <div className="target-actions">
            {d.module_files[focusedModName]?.genes && <a href={file(d.module_files[focusedModName].genes)}>All module genes</a>}
            {d.module_files[focusedModName]?.events && <a href={file(d.module_files[focusedModName].events)}>All module events</a>}
            <a href={file("all_module_lists.zip")}>All module lists ZIP</a>
            <a href={file("analysis_config.json")}>Analysis settings</a>
            <a href={file("network_edges.tsv")}>Complete network edges</a>
            <a href={file("module_members.tsv")}>Complete network nodes</a>
          </div>
        </section>
      )}
    </div>
  );
}
