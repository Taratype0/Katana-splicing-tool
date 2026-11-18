import {useEffect, useRef, useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {api, type GeneCardPayload, type LiveTargetRankerResult} from "../api";
import type {ManifestContrast, TargetBasketEntry} from "../types";

function StripBoxPlot({
  title,
  subtitle,
  unit,
  points,
  domainZeroToOne = false,
}: {
  title: string;
  subtitle: string;
  unit: string;
  points: Array<{group: string; value: number; sample: string}>;
  domainZeroToOne?: boolean;
}) {
  if (!points.length) {
    return (
      <div className="mck-track-half">
        <div className="mck-track-head"><strong>{title}</strong><span>{subtitle}</span></div>
        <p className="mck-compact-note">Per-sample values were not recorded in this portable summary; group-level statistics are preserved above.</p>
      </div>
    );
  }
  const groups = Array.from(new Set(points.map(p => p.group)));
  const vals = points.map(p => p.value);
  const rawMin = domainZeroToOne ? 0 : Math.min(0, ...vals);
  const rawMax = domainZeroToOne ? 1 : Math.max(...vals) * 1.15 || 1;
  const W = 400, H = 175, padL = 52, padR = 16, padT = 16, padB = 34;
  const plotW = W - padL - padR, plotH = H - padT - padB;
  const yMap = (v: number) => padT + (1 - (v - rawMin) / Math.max(1e-6, rawMax - rawMin)) * plotH;
  const colors = ["#051C2C", "#2251FF", "#0077C8", "#B83A4B"];

  return (
    <div className="mck-track-half">
      <div className="mck-track-head"><strong>{title}</strong><span>{subtitle}</span></div>
      <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label={title} style={{width: "100%", height: "auto"}}>
        <rect x={padL} y={padT} width={plotW} height={plotH} fill="#F8FAFC" stroke="#E2E8F0" />
        {[0, 0.5, 1].map(frac => {
          const v = rawMin + frac * (rawMax - rawMin);
          const y = yMap(v);
          return (
            <g key={frac}>
              <line x1={padL} y1={y} x2={padL + plotW} y2={y} stroke="#E2E8F0" />
              <text x={padL - 6} y={y + 3} textAnchor="end" fontSize={8.5} fill="#64748B">
                {domainZeroToOne ? v.toFixed(2) : v >= 1000 ? `${(v / 1000).toFixed(1)}k` : v.toFixed(1)}
              </text>
            </g>
          );
        })}
        {groups.map((grp, gi) => {
          const cx = padL + ((gi + 0.5) / Math.max(1, groups.length)) * plotW;
          const gpts = points.filter(p => p.group === grp);
          const gvals = gpts.map(p => p.value).sort((a, b) => a - b);
          const mean = gvals.reduce((a, b) => a + b, 0) / gvals.length;
          const q1 = gvals[Math.floor(gvals.length * 0.25)] ?? mean;
          const q3 = gvals[Math.floor(gvals.length * 0.75)] ?? mean;
          const col = colors[gi % colors.length];
          return (
            <g key={grp}>
              {gvals.length >= 3 && (
                <rect
                  x={cx - 18}
                  y={Math.min(yMap(q1), yMap(q3))}
                  width={36}
                  height={Math.max(2, Math.abs(yMap(q3) - yMap(q1)))}
                  fill={col}
                  fillOpacity={0.14}
                  stroke={col}
                  strokeWidth={1}
                />
              )}
              <line x1={cx - 22} y1={yMap(mean)} x2={cx + 22} y2={yMap(mean)} stroke={col} strokeWidth={2.2} />
              {gpts.slice(0, 80).map((pt, pi) => {
                const jitter = ((pi % 7) - 3) * 3.8;
                return (
                  <circle key={`${pt.sample}:${pi}`} cx={cx + jitter} cy={yMap(pt.value)} r={3.8} fill={col} stroke="#FFFFFF" strokeWidth={0.7}>
                    <title>{`${pt.sample} (${grp}): ${pt.value}`}</title>
                  </circle>
                );
              })}
              <text x={cx} y={H - 10} textAnchor="middle" fontSize={9.2} fontWeight={600} fill="#051C2C">
                {grp.length > 16 ? `${grp.slice(0, 15)}…` : grp} (n={gpts.length})
              </text>
            </g>
          );
        })}
        <text x={12} y={padT + plotH / 2} textAnchor="middle" fontSize={8.5} fill="#051C2C" transform={`rotate(-90 12 ${padT + plotH / 2})`}>
          {unit}
        </text>
      </svg>
    </div>
  );
}

function QuantitativeSashimiTrack({data}: {data: NonNullable<GeneCardPayload["track2_sashimi"]>}) {
  const exons = data.exons || [];
  const junctions = data.junctions || [];
  if (!exons.length) {
    return <p className="mck-compact-note">Exact exon splice-junction geometry is recorded in external annotations; export the server Sashimi task below to render BAM coverage tracks.</p>;
  }
  const W = 860, H = 228;
  const nEx = exons.length;
  const startX = 95, endX = 720;
  const span = endX - startX;
  const exW = Math.min(115, (span * 0.46) / nEx);
  const gapW = nEx > 1 ? (span - nEx * exW) / (nEx - 1) : 0;
  const exonBoxes = exons.map((ex, idx) => ({
    ...ex,
    x1: startX + idx * (exW + gapW),
    x2: startX + idx * (exW + gapW) + exW,
  }));

  const maxReads = Math.max(
    1,
    ...data.arms.flatMap(a => [a.mean_ijc ?? 0, a.mean_sjc ?? 0])
  );

  return (
    <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label="Quantitative splice-junction Sashimi preview" style={{width: "100%", height: "auto"}}>
      {data.arms.map((arm, ai) => {
        const baseY = ai === 0 ? 74 : 158;
        const armColor = ai === 0 ? "#051C2C" : "#2251FF";
        const waveFill = ai === 0 ? "#CBD5E1" : "#BFDBFE";
        const psiVal = arm.mean_psi;

        return (
          <g key={arm.label}>
            <text x={12} y={baseY - 12} fontSize={10.5} fontWeight={700} fill={armColor}>
              {arm.label}
            </text>
            <text x={12} y={baseY + 4} fontSize={8.8} fill="#64748B">
              n = {arm.n_samples}
            </text>
            <line x1={startX - 8} y1={baseY} x2={endX + 8} y2={baseY} stroke="#64748B" strokeWidth={1.1} />
            {exonBoxes.slice(0, -1).map((eb, bi) => {
              const mid = (eb.x2 + exonBoxes[bi + 1].x1) / 2;
              return (
                <text key={bi} x={mid} y={baseY + 3.5} textAnchor="middle" fontSize={9} fill="#94A3B8">
                  {data.strand === "-" ? "◀  ◀" : "▶  ▶"}
                </text>
              );
            })}

            {exonBoxes.map((eb, ei) => {
              const isAlt = eb.role !== "constitutive";
              const relHeight = isAlt ? Math.max(5, (psiVal ?? 0.5) * 28) : 26;
              return (
                <g key={ei}>
                  <path
                    d={`M ${eb.x1} ${baseY} Q ${eb.x1 + 8} ${baseY - relHeight}, ${(eb.x1 + eb.x2) / 2} ${baseY - relHeight} Q ${eb.x2 - 8} ${baseY - relHeight}, ${eb.x2} ${baseY} Z`}
                    fill={isAlt ? (ai === 0 ? "#94A3B8" : "#60A5FA") : waveFill}
                    stroke={armColor}
                    strokeWidth={1.1}
                  />
                  <rect
                    x={eb.x1}
                    y={baseY}
                    width={eb.x2 - eb.x1}
                    height={9}
                    fill={isAlt ? "#2251FF" : "#051C2C"}
                    rx={1}
                  />
                  {ai === 1 && (
                    <text x={(eb.x1 + eb.x2) / 2} y={baseY + 23} textAnchor="middle" fontSize={8.5} fill="#334155">
                      {eb.label} ({eb.end - eb.start} nt)
                    </text>
                  )}
                </g>
              );
            })}

            {junctions.map((j, ji) => {
              const bFrom = exonBoxes[j.from_idx];
              const bTo = exonBoxes[j.to_idx];
              if (!bFrom || !bTo) return null;
              const x1 = bFrom.x2 - 2;
              const x2 = bTo.x1 + 2;
              const isInc = j.role === "inclusion";
              const count = isInc ? arm.mean_ijc : arm.mean_sjc;
              const frac = count != null ? Math.max(0.15, count / maxReads) : isInc ? (psiVal ?? 0.5) : 1 - (psiVal ?? 0.5);
              const strokeW = 1.1 + 3.4 * Math.pow(frac, 0.85);
              const arcH = isInc ? 34 : 52;
              const apexY = baseY - arcH;
              const arcCol = isInc ? "#2251FF" : "#B83A4B";
              const midX = (x1 + x2) / 2;
              const badgeText = count != null ? `${Math.round(count)}` : psiVal != null ? `${Math.round((isInc ? psiVal : 1 - psiVal) * 100)}%` : "";

              return (
                <g key={ji}>
                  <path
                    d={`M ${x1} ${baseY - 4} C ${x1 + 22} ${apexY}, ${x2 - 22} ${apexY}, ${x2} ${baseY - 4}`}
                    fill="none"
                    stroke={arcCol}
                    strokeWidth={strokeW}
                    strokeOpacity={0.9}
                  />
                  {badgeText && (
                    <g>
                      <rect x={midX - 18} y={apexY - 8} width={36} height={14} rx={3} fill="#FFFFFF" stroke={arcCol} strokeWidth={0.8} />
                      <text x={midX} y={apexY + 2} textAnchor="middle" fontSize={8.4} fontWeight={700} fill={arcCol}>
                        {badgeText}
                      </text>
                    </g>
                  )}
                </g>
              );
            })}

            <rect x={736} y={baseY - 34} width={114} height={44} rx={2} fill="#F8FAFC" stroke="#CBD5E1" />
            <text x={793} y={baseY - 16} textAnchor="middle" fontSize={10.5} fontWeight={700} fill={armColor}>
              PSI: {psiVal != null ? psiVal.toFixed(3) : "NA"}
            </text>
            <text x={793} y={baseY - 1} textAnchor="middle" fontSize={8.3} fill="#64748B">
              {arm.mean_ijc != null ? `IJC ${Math.round(arm.mean_ijc)} · SJC ${Math.round(arm.mean_sjc ?? 0)}` : "Single-cell PSI"}
            </text>
          </g>
        );
      })}
      <text x={startX} y={218} fontSize={8.8} fill="#64748B">
        Locus: {data.chrom} ({data.strand} strand) · Compressed-intron exon/junction geometry · Blue arcs = Inclusion · Burgundy arcs = Skipping
      </text>
    </svg>
  );
}

function DualTrackProteomeSvg({track4}: {track4: NonNullable<GeneCardPayload["track4_proteome"]>}) {
  if (track4.status === "noncoding") {
    return (
      <div className="mck-compact-strip">
        <span className="mck-tag">NON-CODING TRANSCRIPT</span>
        <span>
          <strong>{track4.transcript_id || "Annotated lncRNA / processed transcript"}</strong> — {track4.summary || "No translated CDS open reading frame; alternative splicing remodels mature RNA exon structure without protein translation."}
        </span>
      </div>
    );
  }
  if (track4.status !== "coding_architecture") {
    return track4.note ? <p className="mck-compact-note">{track4.note}</p> : null;
  }

  const lenInc = track4.aa_len_inclusion ?? 400;
  const lenExc = track4.aa_len_exclusion ?? 400;
  const maxAa = Math.max(100, lenInc, lenExc);
  const nmdCls = track4.nmd_switch_class || "coding_preserved";
  const incIsNmd = nmdCls === "inclusion_NMD" || nmdCls === "dual_NMD";
  const excIsNmd = nmdCls === "exclusion_NMD" || nmdCls === "dual_NMD";
  const domainsRaw = (track4.disrupted_protein_features || "")
    .split(/[;|]/)
    .map(s => s.trim())
    .filter(Boolean);

  const W = 900, H = 148, x0 = 148, maxBarW = 465;
  const wInc = Math.max(42, (lenInc / maxAa) * maxBarW);
  const wExc = Math.max(42, (lenExc / maxAa) * maxBarW);

  return (
    <div>
      <svg viewBox={`0 0 ${W} ${H}`} role="img" aria-label="Dual-track mature mRNA NMD and UniProt domain surgery" style={{width: "100%", height: "auto"}}>
        <text x={12} y={40} fontSize={10} fontWeight={700} fill="#051C2C">Inclusion Isoform</text>
        <text x={12} y={54} fontSize={8.5} fill="#64748B">{Math.round(lenInc)} aa {track4.uniprot_acc ? `· ${track4.uniprot_acc}` : ""}</text>
        <rect x={x0} y={28} width={wInc} height={22} rx={2} fill={incIsNmd ? "#FEE2E2" : "#DBEAFE"} stroke={incIsNmd ? "#B83A4B" : "#051C2C"} strokeWidth={1.2} />
        {domainsRaw.slice(0, 3).map((dom, di) => {
          const dw = Math.min(115, wInc * 0.28);
          const dx = x0 + Math.min(wInc - dw - 4, wInc * (0.28 + di * 0.24));
          return (
            <g key={di}>
              <rect x={Math.max(x0 + 2, dx)} y={30} width={Math.max(24, dw)} height={18} rx={1.5} fill="#2251FF" />
              <text x={Math.max(x0 + 2, dx) + Math.max(24, dw) / 2} y={42} textAnchor="middle" fontSize={7.6} fontWeight={700} fill="#FFFFFF">
                {dom.slice(0, 19)}
              </text>
            </g>
          );
        })}
        {incIsNmd ? (
          <g>
            <polygon points={`${x0 + wInc},22 ${x0 + wInc - 6},12 ${x0 + wInc + 6},12`} fill="#B83A4B" />
            <text x={x0 + wInc + 10} y={43} fontSize={9} fontWeight={700} fill="#B83A4B">
              PTC ≥ 50 nt upstream of last EJC → NMD Target
            </text>
          </g>
        ) : (
          <text x={x0 + wInc + 10} y={43} fontSize={9} fontWeight={600} fill="#051C2C">
            Coding intact ({Math.round(lenInc)} aa)
          </text>
        )}

        <text x={12} y={102} fontSize={10} fontWeight={700} fill="#051C2C">Skipping Isoform</text>
        <text x={12} y={116} fontSize={8.5} fill="#64748B">{Math.round(lenExc)} aa (Δ = {track4.delta_aa != null ? `${track4.delta_aa >= 0 ? "+" : ""}${Math.round(track4.delta_aa)} aa` : "0 aa"})</text>
        <rect x={x0} y={90} width={wExc} height={22} rx={2} fill={excIsNmd ? "#FEE2E2" : "#E2E8F0"} stroke={excIsNmd ? "#B83A4B" : "#334155"} strokeWidth={1.2} />
        {excIsNmd ? (
          <g>
            <polygon points={`${x0 + wExc},84 ${x0 + wExc - 6},74 ${x0 + wExc + 6},74`} fill="#B83A4B" />
            <text x={x0 + wExc + 10} y={105} fontSize={9} fontWeight={700} fill="#B83A4B">
              PTC ≥ 50 nt upstream of last EJC → NMD Target
            </text>
          </g>
        ) : (
          <text x={x0 + wExc + 10} y={105} fontSize={9} fontWeight={600} fill="#334155">
            {domainsRaw.length ? `Altered / host domain: ${domainsRaw[0].slice(0, 38)}` : `Coding intact (${Math.round(lenExc)} aa)`}
          </text>
        )}
      </svg>
      {track4.summary && <p className="mck-hint" style={{marginTop: 4}}>{track4.summary}</p>}
    </div>
  );
}

export function GeneCardDossier({
  projectId,
  contrasts,
  entries,
  initialGene = "",
  initialEvent = "",
  initialComparison = "",
}: {
  projectId: string;
  contrasts: ManifestContrast[];
  entries: TargetBasketEntry[];
  initialGene?: string;
  initialEvent?: string;
  initialComparison?: string;
}) {
  const basketCandidates = entries.flatMap(entry =>
    entry.labels.map((lbl, idx) => ({
      label: lbl,
      id: entry.row_ids[idx] || lbl,
      kind: entry.target_kind,
    }))
  );

  const [comparison, setComparison] = useState(initialComparison || contrasts[0]?.contrast_id || "");
  const [geneInput, setGeneInput] = useState(initialGene);
  const [activeGene, setActiveGene] = useState(initialGene);
  const [activeEvent, setActiveEvent] = useState(initialEvent);
  const [seqInput, setSeqInput] = useState("");
  const [committedSeq, setCommittedSeq] = useState("");
  const dossierRef = useRef<HTMLDivElement>(null);

  // Live 5-Pillar Multi-Evidence Target Ranker State
  const [showRanker, setShowRanker] = useState(true);
  const [weights, setWeights] = useState({
    w_dpsi: 1.0,
    w_isoform: 1.0,
    w_splicenet: 1.0,
    w_motif: 1.0,
    w_proteome: 1.5,
    deg_independent_only: false,
    require_domain_or_nmd: false,
  });

  useEffect(() => {
    if (initialGene || initialEvent) {
      setGeneInput(initialGene);
      setActiveGene(initialGene);
      setActiveEvent(initialEvent);
    }
  }, [initialGene, initialEvent]);

  const rankerQuery = useQuery<LiveTargetRankerResult | null>({
    queryKey: ["live-target-ranker", projectId, comparison, weights],
    queryFn: () =>
      typeof api.computeTargetRanker === "function"
        ? api.computeTargetRanker(projectId, {comparison, ...weights})
        : Promise.resolve(null),
    enabled: Boolean(typeof api.computeTargetRanker === "function" && comparison),
  });

  const query = useQuery<GeneCardPayload>({
    queryKey: ["genecard-preview", projectId, activeGene, activeEvent, comparison, committedSeq],
    queryFn: () =>
      typeof api.genecardPreview === "function"
        ? api.genecardPreview(projectId, {
            gene_query: activeGene,
            event_id: activeEvent,
            comparison,
            custom_sequence: committedSeq,
          })
        : Promise.resolve({found: false}),
    enabled: typeof api.genecardPreview === "function",
  });

  if (typeof api.genecardPreview !== "function") return null;

  function exportCardHtml() {
    if (!dossierRef.current) return;
    const htmlContent = `<!doctype html><html><head><meta charset="utf-8"><title>Katana GeneCard · ${query.data?.gene_symbol || "Target"}</title><style>body{font-family:'IBM Plex Sans',Arial,sans-serif;max-width:1050px;margin:28px auto;color:#051C2C}.mck-dossier-track{border:1px solid #E2E8F0;border-top:2px solid #051C2C;padding:16px;margin-bottom:16px}</style></head><body>${dossierRef.current.innerHTML}</body></html>`;
    const blob = new Blob([htmlContent], {type: "text/html"});
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `Katana_GeneCard_${query.data?.gene_symbol || "Target"}_${comparison}.html`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  const d = query.data;
  const ranked = rankerQuery.data;

  return (
    <section className="result-card mck-panel" aria-label="4-Track Adaptive GeneCard & Sashimi Preview">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">CH08 MULTI-EVIDENCE TARGET SCORER &amp; GUARANTEED 4-TRACK MOLECULAR DOSSIER</span>
          <h3>Live 5-Pillar Target Ranker, 4-Track GeneCard &amp; Quantitative Sashimi Studio</h3>
        </div>
        <div className="mck-filter-pills">
          {typeof api.computeTargetRanker === "function" && (
            <button type="button" className="secondary-action compact" onClick={() => setShowRanker(s => !s)}>
              {showRanker ? "Hide 5-Pillar Ranker" : "Show 5-Pillar Ranker"}
            </button>
          )}
          <button type="button" className="secondary-action compact" onClick={exportCardHtml} disabled={!d?.found}>
            Export Standalone GeneCard HTML
          </button>
        </div>
      </div>

      {/* LIVE ENGINE: Interactive 5-Pillar Multi-Evidence Weighting & Target Ranker */}
      {typeof api.computeTargetRanker === "function" && showRanker && ranked?.available && (
        <div className="mck-subcard" style={{marginBottom: 14}}>
          <div className="mck-toolbar-row">
            <h4>Live 5-Pillar Evidence Weighting Scorer ({ranked.total_qualified} qualified events in {comparison})</h4>
            <div style={{display: "flex", gap: 12, fontSize: 12}}>
              <label style={{display: "inline-flex", alignItems: "center", gap: 5}}>
                <input
                  type="checkbox"
                  checked={weights.deg_independent_only}
                  onChange={e => setWeights(w => ({...w, deg_independent_only: e.target.checked}))}
                />
                DEG-independent only
              </label>
              <label style={{display: "inline-flex", alignItems: "center", gap: 5}}>
                <input
                  type="checkbox"
                  checked={weights.require_domain_or_nmd}
                  onChange={e => setWeights(w => ({...w, require_domain_or_nmd: e.target.checked}))}
                />
                Require UniProt Domain / NMD
              </label>
            </div>
          </div>
          <div className="mck-slider-strip">
            {([
              ["w_dpsi", "|ΔPSI| Weight"],
              ["w_isoform", "DTU / DEXSeq"],
              ["w_splicenet", "SpliceNet |kME|"],
              ["w_motif", "Cis-RBP Motif"],
              ["w_proteome", "UniProt / NMD"],
            ] as const).map(([key, label]) => (
              <label key={key}>
                <span>{label}: <strong>{weights[key].toFixed(1)}×</strong></span>
                <input
                  type="range"
                  min={0}
                  max={3}
                  step={0.5}
                  value={weights[key]}
                  onChange={e => setWeights(w => ({...w, [key]: Number(e.target.value)}))}
                />
              </label>
            ))}
          </div>
          <div className="mck-filter-pills" style={{marginTop: 8}}>
            <small style={{fontWeight: 600, color: "#051C2C"}}>Top Ranked Priority Targets (1-Click Load):</small>
            {ranked.rows.slice(0, 10).map(r => (
              <button
                key={r.event_id}
                type="button"
                className={`mck-chip ${d?.event_id === r.event_id ? "active" : ""}`}
                onClick={() => {
                  setGeneInput(r.gene_symbol);
                  setActiveGene(r.gene_symbol);
                  setActiveEvent(r.event_id);
                }}
              >
                {r.gene_symbol} · {r.composite_score.toFixed(2)} pts
              </button>
            ))}
          </div>
        </div>
      )}

      {/* Quick Selector Bar: Basket Items + Instant Gene/Event Lookup + Comparison */}
      <div className="mck-controls-bar">
        <label>
          Comparison
          <select aria-label="GeneCard comparison" value={comparison} onChange={e => setComparison(e.target.value)}>
            {contrasts.map(c => <option key={c.contrast_id} value={c.contrast_id}>{c.label || c.contrast_id}</option>)}
          </select>
        </label>
        <label style={{minWidth: 240}}>
          Inspect Gene or Event ID
          <div style={{display: "flex", gap: 6}}>
            <input
              aria-label="Inspect gene symbol or event ID"
              placeholder="e.g. Macf1, PTPRC, SRSF7…"
              value={geneInput}
              onChange={e => setGeneInput(e.target.value)}
              onKeyDown={e => {
                if (e.key === "Enter") {
                  setActiveGene(geneInput.trim());
                  setActiveEvent("");
                }
              }}
            />
            <button type="button" className="primary-action compact" onClick={() => { setActiveGene(geneInput.trim()); setActiveEvent(""); }}>
              Load
            </button>
          </div>
        </label>
        {!!d?.host_events && d.host_events.length > 1 && (
          <label>
            Events on {d.gene_symbol} ({d.host_events.length})
            <select aria-label="Host gene splicing events" value={d.event_id} onChange={e => setActiveEvent(e.target.value)}>
              {d.host_events.map(ev => (
                <option key={ev.event_id} value={ev.event_id}>
                  {ev.event_type} · ΔPSI={ev.dpsi != null ? ev.dpsi.toFixed(3) : "NA"} · {ev.event_id.slice(0, 38)}
                </option>
              ))}
            </select>
          </label>
        )}
      </div>

      {!!basketCandidates.length && (
        <div className="mck-filter-pills" style={{marginBottom: 12}}>
          <small style={{fontWeight: 600, color: "#051C2C"}}>Selected Basket Targets:</small>
          {basketCandidates.slice(0, 18).map((c, idx) => (
            <button
              key={`${c.label}:${idx}`}
              type="button"
              className={`mck-chip ${d?.gene_symbol === c.label.split(" ")[0] ? "active" : ""}`}
              onClick={() => {
                setGeneInput(c.label.split(" ")[0]);
                setActiveGene(c.label.split(" ")[0]);
                setActiveEvent("");
              }}
            >
              {c.label}
            </button>
          ))}
        </div>
      )}

      {query.isLoading && <p>Assembling 4-Track GeneCard evidence…</p>}
      {d && !d.found && <p className="mck-compact-note">{d.reason}</p>}

      {d?.found && (
        <div ref={dossierRef} className="mck-dossier-stack">
          <div className="mck-dossier-banner">
            <div>
              <h4>{d.gene_symbol} <small>({d.gene_id})</small></h4>
              <span>Event: <code>{d.event_id}</code> ({d.event_type})</span>
            </div>
            <div className="mck-kpi-strip">
              <span><strong>ΔPSI:</strong> {d.dpsi != null ? `${d.dpsi >= 0 ? "+" : ""}${d.dpsi.toFixed(3)}` : "NA"}</span>
              <span><strong>DAS FDR:</strong> {d.fdr != null ? d.fdr.toExponential(2) : "NA"}</span>
              <span><strong>Host DEG log2FC:</strong> {d.track1_quantification?.deg.log2fc != null ? `${d.track1_quantification.deg.log2fc >= 0 ? "+" : ""}${d.track1_quantification.deg.log2fc.toFixed(2)}` : "NS / NA"}</span>
            </div>
          </div>

          {/* TRACK 1: Host Gene RNA Expression + Event Replicate PSI */}
          <article className="mck-dossier-track">
            <div className="mck-track-title">TRACK 1 · HOST GENE RNA EXPRESSION &amp; REPLICATE SPLICING PSI</div>
            <div className="mck-split-grid">
              <StripBoxPlot
                title="Host Gene RNA Expression"
                subtitle={d.track1_quantification?.deg.log2fc != null ? `log2FC = ${d.track1_quantification.deg.log2fc.toFixed(2)} · FDR = ${(d.track1_quantification.deg.fdr ?? 1).toExponential(2)}` : "Saved sample expression"}
                unit={d.track1_quantification?.expression_samples[0]?.unit || "Normalized expression"}
                points={(d.track1_quantification?.expression_samples || []).map(s => ({group: s.arm || s.condition, value: s.value, sample: s.sample_id}))}
              />
              <StripBoxPlot
                title={`Splicing Event Inclusion (PSI · ${d.event_type})`}
                subtitle={`ΔPSI = ${d.dpsi != null ? d.dpsi.toFixed(3) : "NA"} · FDR = ${d.fdr != null ? d.fdr.toExponential(2) : "NA"}`}
                unit="Inclusion level (PSI)"
                domainZeroToOne
                points={(d.track1_quantification?.psi_samples || []).map(s => ({group: s.condition, value: s.psi, sample: s.sample_id}))}
              />
            </div>
          </article>

          {/* TRACK 2: Quantitative Splice-Junction Sashimi Preview */}
          {d.track2_sashimi && (
            <article className="mck-dossier-track">
              <div className="mck-track-title">
                TRACK 2 · QUANTITATIVE SPLICE-JUNCTION SASHIMI PREVIEW
                <small style={{marginLeft: "auto", fontWeight: 400, color: "#64748B"}}>
                  Instant local preview from REPORT.zip · Use "Pipeline task / Run locally" below for full base-resolution BAM RPKM Sashimi
                </small>
              </div>
              <QuantitativeSashimiTrack data={d.track2_sashimi} />
            </article>
          )}

          {/* TRACK 3: Guaranteed Positional RBP Cis-Motif & Live 170-PWM Sequence Scanner */}
          {d.track3_motifs && (
            <article className="mck-dossier-track">
              <div className="mck-track-title">
                TRACK 3 · POSITIONAL RBP CIS-MOTIF &amp; LIVE 170-PWM SEQUENCE EVIDENCE
                {d.track3_motifs.mode && (
                  <small style={{marginLeft: "auto", fontWeight: 400, color: "#64748B"}}>
                    {d.track3_motifs.mode === "live_sequence_pwm_scan" ? "Live 170-PWM Sequence Scan" : d.track3_motifs.mode === "live_regulatory_context" ? "Bundled 170-PWM Regulatory Context" : "Precomputed Katana-SpliceTrace"}
                  </small>
                )}
              </div>
              {d.track3_motifs.hits && d.track3_motifs.hits.length > 0 ? (
                <div className="mck-motif-badges">
                  {d.track3_motifs.hits.map((h, idx) => (
                    <div key={idx} className="mck-motif-hit-card">
                      <strong>{h.rbp}</strong>
                      <span>{h.region}</span>
                      {h.delta_score != null && <small>Score / Shift: {h.delta_score >= 0 ? `+${h.delta_score.toFixed(2)}` : h.delta_score.toFixed(2)}</small>}
                    </div>
                  ))}
                </div>
              ) : (
                <p className="mck-compact-note">{d.track3_motifs.note || "No pre-saved high-threshold motif hit; paste exon/flanking sequence below to scan 170 RBP PWMs live."}</p>
              )}
              <div style={{display: "flex", gap: 6, marginTop: 10, alignItems: "center"}}>
                <input
                  aria-label="Live Track 3 RNA sequence scan"
                  placeholder="Optional live Track 3 scan: paste ±200 nt exon/intron RNA sequence (e.g. GGAUGGACUUGCAUGCUU…) to score 170 PWMs"
                  value={seqInput}
                  onChange={e => setSeqInput(e.target.value)}
                  style={{flex: 1, padding: "6px 10px", fontSize: 12, border: "1px solid #CBD5E1"}}
                />
                <button type="button" className="secondary-action compact" onClick={() => setCommittedSeq(seqInput.trim())}>
                  Scan 170 PWMs
                </button>
                {committedSeq && (
                  <button type="button" className="ghost-action" onClick={() => { setSeqInput(""); setCommittedSeq(""); }}>
                    Reset
                  </button>
                )}
              </div>
            </article>
          )}

          {/* TRACK 4: Guaranteed Mature mRNA 50-nt EJC NMD & Bundled UniProt / Pfam Domain Track */}
          {d.track4_proteome && (
            <article className="mck-dossier-track">
              <div className="mck-track-title">
                TRACK 4 · MATURE mRNA 50-NT EJC NMD &amp; UNIPROT / PFAM DOMAIN ARCHITECTURE
                {d.track4_proteome.computed_mode && (
                  <small style={{marginLeft: "auto", fontWeight: 400, color: "#64748B"}}>
                    Live Exon Coordinate &amp; Bundled UniProt/Pfam Inference
                  </small>
                )}
              </div>
              <DualTrackProteomeSvg track4={d.track4_proteome} />
            </article>
          )}
        </div>
      )}
    </section>
  );
}
