import {useMemo, useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {api} from "../api";
import type {ManifestContrast} from "../types";

export function QuadrantExplorer({
  projectId,
  comparison,
  onInspectGeneCard,
}: {
  projectId: string;
  comparison: string;
  onInspectGeneCard: (geneSymbol: string, eventId?: string) => void;
}) {
  const [fcCut, setFcCut] = useState(1.0);
  const [dpsiCut, setDpsiCut] = useState(0.10);
  const [fdrCut, setFdrCut] = useState(0.05);
  const [regimeFilter, setRegimeFilter] = useState<"all" | "splicing_only" | "dual" | "nmd_domain">("splicing_only");

  const query = useQuery({
    queryKey: ["quadrant-explorer", projectId, comparison],
    queryFn: () => api.quadrantExplorer(projectId, comparison),
  });

  const rows = query.data?.rows || [];
  const classified = useMemo(() => {
    let spOnly = 0, dual = 0, degOnly = 0;
    const mapped = rows.map(r => {
      const sigDeg = r.deg_fdr_val <= fdrCut && Math.abs(r.log2fc_val) >= fcCut;
      const sigDas = r.rs && r.das_fdr_val <= fdrCut && Math.abs(r.dpsi_val) >= dpsiCut;
      const regime = sigDeg && sigDas ? "dual" : sigDas ? "splicing_only" : sigDeg ? "deg_only" : "bg";
      if (regime === "splicing_only") spOnly++;
      else if (regime === "dual") dual++;
      else if (regime === "deg_only") degOnly++;
      return {...r, regime};
    });
    const orthPct = spOnly + dual > 0 ? ((spOnly / (spOnly + dual)) * 100).toFixed(1) : "0.0";
    return {mapped, spOnly, dual, degOnly, orthPct};
  }, [rows, fcCut, dpsiCut, fdrCut]);

  if (query.isLoading) return <p>Loading Expression × Splicing quadrant…</p>;
  if (!rows.length) return null;

  const tableRows = classified.mapped.filter(r => {
    if (regimeFilter === "splicing_only") return r.regime === "splicing_only";
    if (regimeFilter === "dual") return r.regime === "dual";
    if (regimeFilter === "nmd_domain") {
      return (r.regime === "splicing_only" || r.regime === "dual") &&
        (Boolean(r.disrupted_protein_features) || (r.nmd_switch_class && r.nmd_switch_class.includes("NMD")));
    }
    return r.regime !== "bg";
  }).sort((a, b) => Math.abs(b.dpsi_val) - Math.abs(a.dpsi_val));

  const W = 560, H = 310, padL = 52, padR = 20, padT = 18, padB = 38;
  const plotW = W - padL - padR, plotH = H - padT - padB;
  const xMap = (fc: number) => padL + ((Math.max(-6, Math.min(6, fc)) + 6) / 12) * plotW;
  const yMap = (dp: number) => padT + ((1 - Math.max(-1, Math.min(1, dp))) / 2) * plotH;

  return (
    <section className="result-card mck-panel" aria-label="Expression vs Splicing Decoupling Quadrant">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">TRANSCRIPTIONAL VS POST-TRANSCRIPTIONAL DECOUPLING</span>
          <h3>Interactive DEG × PS-DAS Quadrant Explorer ({comparison})</h3>
        </div>
        <span className="mck-stat-pill">
          {classified.orthPct}% of PS-DAS genes are DEG-independent ({classified.spOnly} splicing-only · {classified.dual} dual)
        </span>
      </div>

      <div className="mck-controls-bar">
        <label>
          |log2FC| ≥ <strong>{fcCut.toFixed(2)}</strong>
          <input type="range" min={0.25} max={3.0} step={0.25} value={fcCut} onChange={e => setFcCut(Number(e.target.value))} />
        </label>
        <label>
          |ΔPSI| ≥ <strong>{dpsiCut.toFixed(2)}</strong>
          <input type="range" min={0.05} max={0.50} step={0.05} value={dpsiCut} onChange={e => setDpsiCut(Number(e.target.value))} />
        </label>
        <label>
          FDR ≤
          <select value={String(fdrCut)} onChange={e => setFdrCut(Number(e.target.value))}>
            <option value="0.01">0.01</option>
            <option value="0.05">0.05</option>
            <option value="0.10">0.10</option>
          </select>
        </label>
        <div className="mck-filter-pills">
          {([
            ["splicing_only", `Splicing-Only (${classified.spOnly})`],
            ["dual", `Dual DEG+DAS (${classified.dual})`],
            ["nmd_domain", "Domain / NMD Altered"],
            ["all", "All Significant"],
          ] as const).map(([id, lbl]) => (
            <button key={id} type="button" className={`mck-chip ${regimeFilter === id ? "active" : ""}`} onClick={() => setRegimeFilter(id)}>
              {lbl}
            </button>
          ))}
        </div>
      </div>

      <div className="mck-split-grid">
        <div className="mck-subcard">
          <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label="Interactive DEG vs DAS quadrant plot" style={{width: "100%", height: "auto"}}>
            <rect x={padL} y={padT} width={plotW} height={plotH} fill="#F8FAFC" stroke="#CBD5E1" />
            <line x1={xMap(0)} y1={padT} x2={xMap(0)} y2={padT + plotH} stroke="#CBD5E1" />
            <line x1={padL} y1={yMap(0)} x2={padL + plotW} y2={yMap(0)} stroke="#CBD5E1" />
            <line x1={xMap(-fcCut)} y1={padT} x2={xMap(-fcCut)} y2={padT + plotH} stroke="#94A3B8" strokeDasharray="4,3" />
            <line x1={xMap(fcCut)} y1={padT} x2={xMap(fcCut)} y2={padT + plotH} stroke="#94A3B8" strokeDasharray="4,3" />
            <line x1={padL} y1={yMap(dpsiCut)} x2={padL + plotW} y2={yMap(dpsiCut)} stroke="#94A3B8" strokeDasharray="4,3" />
            <line x1={padL} y1={yMap(-dpsiCut)} x2={padL + plotW} y2={yMap(-dpsiCut)} stroke="#94A3B8" strokeDasharray="4,3" />

            {classified.mapped.slice(0, 2200).map(r => {
              const fill = r.regime === "dual" ? "#B83A4B" : r.regime === "splicing_only" ? "#2251FF" : r.regime === "deg_only" ? "#94A3B8" : "#E2E8F0";
              const rad = r.regime === "dual" ? 4.2 : r.regime === "splicing_only" ? 3.6 : 2.1;
              return (
                <circle
                  key={r.gene_id}
                  cx={xMap(r.log2fc_val)}
                  cy={yMap(r.dpsi_val)}
                  r={rad}
                  fill={fill}
                  fillOpacity={r.regime === "bg" ? 0.45 : 0.88}
                  style={{cursor: "pointer"}}
                  onClick={() => onInspectGeneCard(r.gene_symbol, r.event_id)}
                >
                  <title>{`${r.gene_symbol} · log2FC=${r.log2fc_val.toFixed(2)} · ΔPSI=${r.dpsi_val.toFixed(3)}`}</title>
                </circle>
              );
            })}
            <text x={padL + plotW / 2} y={H - 6} textAnchor="middle" fontSize={10.5} fill="#051C2C">Host Gene Expression log2FC (clipped ±6)</text>
            <text x={14} y={padT + plotH / 2} textAnchor="middle" fontSize={10.5} fill="#051C2C" transform={`rotate(-90 14 ${padT + plotH / 2})`}>Strongest Event ΔPSI</text>
          </svg>
          <p className="mck-hint">Click any gene point to open its 4-Track GeneCard &amp; Sashimi dossier directly.</p>
        </div>

        <div className="mck-subcard">
          <h4>Matching Quadrant Genes ({tableRows.length})</h4>
          <div className="report-table-wrap" style={{maxHeight: 275}}>
            <table>
              <thead>
                <tr>
                  <th>Gene</th>
                  <th className="num">log2FC</th>
                  <th className="num">ΔPSI</th>
                  <th>Consequence</th>
                  <th>Inspect</th>
                </tr>
              </thead>
              <tbody>
                {tableRows.slice(0, 60).map(r => (
                  <tr key={r.gene_id}>
                    <td><strong>{r.gene_symbol}</strong></td>
                    <td className="num">{r.log2fc_val.toFixed(2)}</td>
                    <td className="num">{r.dpsi_val >= 0 ? `+${r.dpsi_val.toFixed(3)}` : r.dpsi_val.toFixed(3)}</td>
                    <td><small>{r.disrupted_protein_features || r.nmd_switch_class || r.regime}</small></td>
                    <td>
                      <button type="button" className="mck-link-btn" onClick={() => onInspectGeneCard(r.gene_symbol, r.event_id)}>
                        GeneCard →
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </section>
  );
}

export function ContrastRescueExplorer({
  projectId,
  contrasts,
  onInspectGeneCard,
}: {
  projectId: string;
  contrasts: ManifestContrast[];
  onInspectGeneCard: (geneSymbol: string, eventId?: string) => void;
}) {
  const [contrastA, setContrastA] = useState(contrasts[0]?.contrast_id || "");
  const [contrastB, setContrastB] = useState(contrasts[1]?.contrast_id || contrasts[0]?.contrast_id || "");
  const [relation, setRelation] = useState<"reversed" | "concordant" | "specific_to_a">("reversed");
  const [minDpsi, setMinDpsi] = useState(0.10);
  const [maxFdr, setMaxFdr] = useState(0.05);

  const query = useQuery({
    queryKey: ["contrast-rescue", projectId, contrastA, contrastB, relation, minDpsi, maxFdr],
    queryFn: () => api.contrastRescue(projectId, {
      contrast_a: contrastA,
      contrast_b: contrastB,
      min_abs_dpsi: minDpsi,
      max_fdr: maxFdr,
      relation,
      rs_only: true,
    }),
    enabled: contrasts.length >= 2 && Boolean(contrastA) && Boolean(contrastB),
  });

  if (contrasts.length < 2) return null;
  const rows = query.data?.rows || [];

  return (
    <section className="result-card mck-panel" aria-label="Multi-contrast rescue and concordance explorer">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">MULTI-COMPARISON DIRECTIONAL LOGIC</span>
          <h3>Interactive Rescue / Reversal &amp; Contrast-Specific Event Filter</h3>
        </div>
        <span className="mck-stat-pill">{query.data?.total ?? 0} matching PS-DAS events</span>
      </div>

      <div className="mck-controls-bar">
        <label>
          Comparison A
          <select aria-label="Comparison A" value={contrastA} onChange={e => setContrastA(e.target.value)}>
            {contrasts.map(c => <option key={c.contrast_id} value={c.contrast_id}>{c.label || c.contrast_id}</option>)}
          </select>
        </label>
        <label>
          Comparison B
          <select aria-label="Comparison B" value={contrastB} onChange={e => setContrastB(e.target.value)}>
            {contrasts.map(c => <option key={c.contrast_id} value={c.contrast_id}>{c.label || c.contrast_id}</option>)}
          </select>
        </label>
        <label>
          Directional Logic
          <select aria-label="Directional relation" value={relation} onChange={e => setRelation(e.target.value as typeof relation)}>
            <option value="reversed">Reversed / Rescued (ΔPSI opposite sign in B)</option>
            <option value="concordant">Concordant (ΔPSI same sign in A &amp; B)</option>
            <option value="specific_to_a">Specific to A (unchanged in B)</option>
          </select>
        </label>
        <label>
          |ΔPSI| ≥ <strong>{minDpsi.toFixed(2)}</strong>
          <input type="range" min={0.05} max={0.40} step={0.05} value={minDpsi} onChange={e => setMinDpsi(Number(e.target.value))} />
        </label>
        <label>
          FDR ≤
          <select value={String(maxFdr)} onChange={e => setMaxFdr(Number(e.target.value))}>
            <option value="0.01">0.01</option>
            <option value="0.05">0.05</option>
            <option value="0.10">0.10</option>
          </select>
        </label>
      </div>

      <div className="report-table-wrap" style={{maxHeight: 290}}>
        <table>
          <thead>
            <tr>
              <th>Gene</th>
              <th>Type</th>
              <th className="num">ΔPSI ({contrastA})</th>
              <th className="num">ΔPSI ({contrastB})</th>
              <th>NMD / Domain Surgery</th>
              <th>Action</th>
            </tr>
          </thead>
          <tbody>
            {rows.map(r => (
              <tr key={r.event_id}>
                <td><strong>{r.gene_symbol}</strong></td>
                <td>{r.event_type || "SE"}</td>
                <td className="num">{r.dpsi_a >= 0 ? `+${r.dpsi_a.toFixed(3)}` : r.dpsi_a.toFixed(3)}</td>
                <td className="num">{r.dpsi_b != null ? (r.dpsi_b >= 0 ? `+${r.dpsi_b.toFixed(3)}` : r.dpsi_b.toFixed(3)) : "NS / NA"}</td>
                <td><small>{r.disrupted_protein_features || r.nmd_switch_class || "—"}</small></td>
                <td>
                  <button type="button" className="mck-link-btn" onClick={() => onInspectGeneCard(r.gene_symbol, r.event_id)}>
                    Open GeneCard →
                  </button>
                </td>
              </tr>
            ))}
            {!rows.length && (
              <tr><td colSpan={6}>No PS-DAS events satisfy this directional comparison rule at the selected thresholds.</td></tr>
            )}
          </tbody>
        </table>
      </div>
    </section>
  );
}
