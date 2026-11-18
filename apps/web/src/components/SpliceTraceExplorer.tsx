import {useMemo, useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {api, type LiveMotifScanResult} from "../api";

const REGION_ORDER = [
  ["upstream_intron_85_200", "Up Int [-200..-85]"],
  ["upstream_intron_35_85", "Up Int [-85..-35]"],
  ["upstream_intron_0_35", "3′SS [-35..-3]"],
  ["upstream_exon_3p_50", "Up Exon 3′"],
  ["alt_exon_5p_50", "Alt Exon 5′"],
  ["alt_exon_3p_50", "Alt Exon 3′"],
  ["downstream_exon_5p_50", "Dn Exon 5′"],
  ["downstream_intron_0_35", "5′SS [+3..+35]"],
  ["downstream_intron_35_85", "Dn Int [+35..+85]"],
  ["downstream_intron_85_200", "Dn Int [+85..+200]"],
  ["mature_consequence_delta", "Mature Δ(I−E)"],
];

const BASE_COLORS: Record<string, string> = {
  A: "#0D9488",
  C: "#2251FF",
  G: "#D04A02",
  U: "#B83A4B",
};

export function PwmSequenceLogo({matrix, motifId}: {matrix: number[][]; motifId: string}) {
  const bases = ["A", "C", "G", "U"];
  const cols = matrix.map(row => {
    const sum = row.reduce((a, b) => a + Math.max(0, b), 0) || 1;
    const p = row.map(v => Math.max(1e-6, v / sum));
    const h = -p.reduce((acc, pi) => acc + pi * Math.log2(pi), 0);
    const ic = Math.max(0, 2 - h);
    const letters = bases
      .map((b, idx) => ({base: b, bits: p[idx] * ic}))
      .sort((a, b) => a.bits - b.bits);
    return {ic, letters};
  });
  const colW = 30;
  const width = Math.max(200, cols.length * colW + 44);
  const height = 112;
  const plotH = 76;

  return (
    <div className="mck-logo-card">
      <div className="mck-sublabel">Curated RNA PWM Sequence Logo · 5′ → 3′ ({motifId})</div>
      <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`PWM logo ${motifId}`} style={{maxWidth: 420, height: "auto"}}>
        <line x1={32} y1={12} x2={32} y2={12 + plotH} stroke="#64748B" strokeWidth={1} />
        <line x1={32} y1={12 + plotH} x2={width - 8} y2={12 + plotH} stroke="#64748B" strokeWidth={1} />
        <text x={26} y={18} textAnchor="end" fontSize={9} fill="#64748B">2.0</text>
        <text x={26} y={12 + plotH / 2} textAnchor="end" fontSize={9} fill="#64748B">1.0</text>
        <text x={26} y={12 + plotH} textAnchor="end" fontSize={9} fill="#64748B">0</text>
        <text x={10} y={12 + plotH / 2} textAnchor="middle" fontSize={9} fill="#051C2C" transform={`rotate(-90 10 ${12 + plotH / 2})`}>bits</text>
        {cols.map((col, ci) => {
          let yCursor = 12 + plotH;
          return (
            <g key={ci}>
              {col.letters.map(item => {
                const hPx = (item.bits / 2.0) * plotH;
                yCursor -= hPx;
                if (hPx < 2.5) return null;
                return (
                  <text
                    key={item.base}
                    x={36 + ci * colW + colW / 2}
                    y={yCursor + hPx * 0.86}
                    textAnchor="middle"
                    fontFamily="IBM Plex Sans, Arial, sans-serif"
                    fontWeight={800}
                    fontSize={Math.max(8, hPx * 1.08)}
                    fill={BASE_COLORS[item.base]}
                  >
                    {item.base}
                  </text>
                );
              })}
              <text x={36 + ci * colW + colW / 2} y={102} textAnchor="middle" fontSize={8.5} fill="#64748B">{ci + 1}</text>
            </g>
          );
        })}
      </svg>
    </div>
  );
}

export function SpliceTraceExplorer({
  projectId,
  comparison,
  onInspectGeneCard,
}: {
  projectId: string;
  comparison: string;
  onInspectGeneCard: (geneSymbol: string, eventId?: string) => void;
}) {
  const query = useQuery({
    queryKey: ["splicetrace-workbench", projectId, comparison],
    queryFn: () =>
      typeof api.splicetraceWorkbench === "function"
        ? api.splicetraceWorkbench(projectId, comparison)
        : Promise.resolve({available: false, comparison_id: comparison, rbps: [], regions: {}, targets: [], bridge: []}),
  });

  const [selectedRbp, setSelectedRbp] = useState("");
  const [consequenceFilter, setConsequenceFilter] = useState<"all" | "domain" | "nmd" | "motif">("all");
  const [search, setSearch] = useState("");

  // Live Custom IUPAC / 170-RBP PWM Scanner State
  const [motifInput, setMotifInput] = useState("UGGACU");
  const [seqInput, setSeqInput] = useState("");
  const [committedScan, setCommittedScan] = useState({motif: "UGGACU", seq: ""});

  const scanQuery = useQuery<LiveMotifScanResult | null>({
    queryKey: ["live-motif-scan", projectId, comparison, committedScan],
    queryFn: () =>
      typeof api.computeMotifScan === "function"
        ? api.computeMotifScan(projectId, {
            comparison,
            query_motif_or_rbp: committedScan.motif,
            custom_sequence: committedScan.seq,
          })
        : Promise.resolve(null),
  });

  const data = query.data;
  const activeRbp = data?.rbps?.find(r => r.rbp === selectedRbp) || data?.rbps?.[0];
  const rbpName = activeRbp?.rbp || "";
  const profile = (rbpName && data?.regions?.[rbpName]) || [];

  const filteredTargets = useMemo(() => {
    const targets = data?.targets || [];
    return targets.filter(t => {
      if (consequenceFilter === "domain" && !t.disrupted_protein_features) return false;
      if (consequenceFilter === "nmd" && (!t.nmd_switch_class.includes("NMD") || t.nmd_switch_class === "dual_NMD")) return false;
      if (consequenceFilter === "motif" && rbpName && !t.cis_rbps.toUpperCase().includes(rbpName.toUpperCase())) return false;
      if (search && !`${t.gene_symbol} ${t.event_id} ${t.disrupted_protein_features} ${t.cis_rbps}`.toLowerCase().includes(search.toLowerCase())) return false;
      return true;
    });
  }, [data?.targets, consequenceFilter, rbpName, search]);

  if (typeof api.splicetraceWorkbench !== "function" && typeof api.computeMotifScan !== "function") return null;
  if (query.isLoading && scanQuery.isLoading) return <p>Loading Katana-SpliceTrace RBP &amp; Motif Lab…</p>;

  const orderedBars = REGION_ORDER.map(([key, label]) => {
    const found = profile.find(p => p.region === key);
    return {
      key,
      label,
      z: found?.stat_z ?? 0,
      p: found?.nominal_p ?? 1,
    };
  });
  const zExtent = Math.max(3.2, ...orderedBars.map(b => Math.abs(b.z)));
  const scan = scanQuery.data;

  return (
    <section className="result-card mck-panel" aria-label="Interactive Katana-SpliceTrace workbench">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">KATANA-SPLICETRACE · SEQUENCE-FIRST CIS-MOTIF, 170-PWM SCANNER &amp; PROTEOME LAB</span>
          <h3>Interactive RBP Positional Map, Live IUPAC / PWM Scanner &amp; Proteome Consequences ({comparison})</h3>
        </div>
        {data?.available && (
          <span className="mck-stat-pill">
            {data.rbps.filter(r => (r.fwer_p ?? 1) <= 0.05).length} FWER ≤ 0.05 RBPs · {data.targets.length} |ΔPSI| ≥ 0.10 target events
          </span>
        )}
      </div>

      {/* LIVE ENGINE: Custom IUPAC Motif & 170-RBP PWM Scanner */}
      {typeof api.computeMotifScan === "function" && (
        <div className="mck-subcard" style={{marginBottom: 16}}>
          <div className="mck-toolbar-row">
            <h4>Live Custom IUPAC RNA Motif &amp; 170-RBP PWM Scanner</h4>
            <div className="mck-filter-pills">
              {([
                ["UGGACU", "UGGACU (m⁶A / YTHDC1)"],
                ["WGCAUG", "WGCAUG (RBFOX)"],
                ["YGCY", "YGCY (MBNL)"],
                ["GGAAGG", "GGAAGG (SRSF1)"],
                ["SRSF7", "SRSF7 PWM"],
                ["PTBP1", "PTBP1 PWM"],
              ] as const).map(([val, lbl]) => (
                <button
                  key={val}
                  type="button"
                  className={`mck-chip ${committedScan.motif === val ? "active" : ""}`}
                  onClick={() => {
                    setMotifInput(val);
                    setCommittedScan(s => ({...s, motif: val}));
                  }}
                >
                  {lbl}
                </button>
              ))}
            </div>
          </div>

          <div className="mck-controls-bar" style={{marginTop: 8}}>
            <label style={{minWidth: 210}}>
              IUPAC Motif or RBP Symbol (170 PWMs)
              <input
                aria-label="Custom IUPAC motif or RBP symbol"
                placeholder="e.g. UGGACU, WGCAUG, RBM10…"
                value={motifInput}
                onChange={e => setMotifInput(e.target.value)}
              />
            </label>
            <label style={{flex: 1, minWidth: 260}}>
              Optional Exon / Intron RNA Sequence (5′ → 3′ live log-odds scan)
              <input
                aria-label="Optional RNA sequence for live PWM scan"
                placeholder="Paste RNA/DNA sequence (e.g. GGAUGGACUUGCAUGCUU…) to score all 170 PWMs live"
                value={seqInput}
                onChange={e => setSeqInput(e.target.value)}
              />
            </label>
            <button
              type="button"
              className="primary-action compact"
              onClick={() => setCommittedScan({motif: motifInput.trim() || "UGGACU", seq: seqInput.trim()})}
              disabled={scanQuery.isFetching}
            >
              {scanQuery.isFetching ? "Scanning…" : "Scan Motif & PWMs"}
            </button>
          </div>

          {scan && (
            <div className="mck-split-grid" style={{marginTop: 10}}>
              <div>
                {scan.pwm_matrix && <PwmSequenceLogo matrix={scan.pwm_matrix} motifId={scan.motif_id} />}
                {!!scan.homologous_rbps.length && (
                  <div style={{marginTop: 8}}>
                    <small style={{fontWeight: 600, color: "#051C2C"}}>Matching Bundled RBP PWMs ({scan.catalog_size} catalog): </small>
                    {scan.homologous_rbps.slice(0, 6).map(h => (
                      <span key={h.rbp} className="mck-tag" style={{marginRight: 4}}>
                        {h.rbp} ({h.consensus})
                      </span>
                    ))}
                  </div>
                )}
              </div>
              <div>
                <div className="mck-kpi-grid" style={{gridTemplateColumns: "repeat(3, minmax(0, 1fr))", marginBottom: 8}}>
                  <div className="mck-kpi-box">
                    <small>Consensus</small>
                    <strong>{scan.consensus || scan.query}</strong>
                    <span>{scan.motif_id}</span>
                  </div>
                  <div className="mck-kpi-box">
                    <small>Inclusion / Skipping</small>
                    <strong>+{scan.inclusion_hits} / −{scan.skipping_hits}</strong>
                    <span>Binomial P = {scan.direction_binomial_p.toPrecision(2)}</span>
                  </div>
                  <div className="mck-kpi-box">
                    <small>Matched Targets</small>
                    <strong>{scan.matched_events.length}</strong>
                    <span>Click any gene for GeneCard</span>
                  </div>
                </div>
                {!!scan.sequence_scan_hits.length && (
                  <div className="mck-motif-badges" style={{marginBottom: 8}}>
                    {scan.sequence_scan_hits.slice(0, 6).map((sh, idx) => (
                      <div key={idx} className="mck-motif-hit-card">
                        <strong>{sh.rbp}</strong>
                        <span>{sh.region}</span>
                        <small>Log-odds: +{sh.delta_score.toFixed(2)} bits</small>
                      </div>
                    ))}
                  </div>
                )}
                {!!scan.matched_events.length && (
                  <div className="mck-filter-pills">
                    <small style={{fontWeight: 600, color: "#051C2C"}}>Top {scan.query} Target Events:</small>
                    {scan.matched_events.slice(0, 8).map(ev => (
                      <button
                        key={ev.event_id}
                        type="button"
                        className="mck-chip"
                        onClick={() => onInspectGeneCard(ev.gene_symbol, ev.event_id)}
                      >
                        {ev.gene_symbol} ({ev.dpsi != null ? `${ev.dpsi >= 0 ? "+" : ""}${ev.dpsi.toFixed(2)}` : ev.event_type}) →
                      </button>
                    ))}
                  </div>
                )}
              </div>
            </div>
          )}
        </div>
      )}

      {data?.available && (
        <>
          <div className="mck-split-grid">
            {/* Left: RBP Leaderboard */}
            <div className="mck-subcard">
              <h4>1. Gene-Block max-T FWER Leaderboard</h4>
              <p className="mck-hint">Click any RBP to inspect its 11-window positional RNA map, curated PWM logo, and target splicing events.</p>
              <div className="report-table-wrap" style={{maxHeight: 330}}>
                <table>
                  <thead>
                    <tr>
                      <th>RBP</th>
                      <th>Peak Window</th>
                      <th className="num">Z</th>
                      <th className="num">FWER P</th>
                      <th className="num">Host log2FC</th>
                    </tr>
                  </thead>
                  <tbody>
                    {data.rbps.map(r => {
                      const isSig = (r.fwer_p ?? 1) <= 0.05;
                      const active = r.rbp === rbpName;
                      return (
                        <tr
                          key={r.rbp}
                          className={active ? "mck-row-active" : ""}
                          onClick={() => setSelectedRbp(r.rbp)}
                          style={{cursor: "pointer"}}
                        >
                          <td>
                            <strong>{r.rbp}</strong>{" "}
                            {isSig && <span className="mck-tag mck-tag-blue">FWER</span>}
                          </td>
                          <td><small>{r.best_region}</small></td>
                          <td className="num">{r.stat_z != null ? r.stat_z.toFixed(2) : "NA"}</td>
                          <td className="num">{r.fwer_p != null ? r.fwer_p.toExponential(2) : r.nominal_p != null ? r.nominal_p.toExponential(2) : "NA"}</td>
                          <td className="num">{r.host_deg_log2fc != null ? `${r.host_deg_log2fc >= 0 ? "+" : ""}${r.host_deg_log2fc.toFixed(2)}` : "—"}</td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
            </div>

            {/* Right: 11-Channel Positional RNA Map + Live PWM Sequence Logo */}
            <div className="mck-subcard">
              <h4>2. Positional Splicing Map &amp; Motif Logo · {rbpName}</h4>
              {orderedBars.some(b => b.z !== 0) ? (
                <svg viewBox="0 0 620 205" role="img" aria-label={`11-channel positional map for ${rbpName}`} style={{width: "100%", height: "auto"}}>
                  <line x1={48} y1={95} x2={605} y2={95} stroke="#051C2C" strokeWidth={1.1} />
                  <line x1={48} y1={95 - (1.96 / zExtent) * 72} x2={605} y2={95 - (1.96 / zExtent) * 72} stroke="#94A3B8" strokeDasharray="4,3" strokeWidth={0.9} />
                  <line x1={48} y1={95 + (1.96 / zExtent) * 72} x2={605} y2={95 + (1.96 / zExtent) * 72} stroke="#94A3B8" strokeDasharray="4,3" strokeWidth={0.9} />
                  <text x={42} y={24} textAnchor="end" fontSize={9} fill="#2251FF">+{zExtent.toFixed(1)} Z</text>
                  <text x={42} y={98} textAnchor="end" fontSize={9} fill="#64748B">0</text>
                  <text x={42} y={168} textAnchor="end" fontSize={9} fill="#B83A4B">−{zExtent.toFixed(1)} Z</text>
                  {orderedBars.map((b, i) => {
                    const x = 58 + i * 49;
                    const barH = (Math.abs(b.z) / zExtent) * 72;
                    const y = b.z >= 0 ? 95 - barH : 95;
                    const fill = b.z >= 0 ? (Math.abs(b.z) >= 1.96 ? "#2251FF" : "#93C5FD") : (Math.abs(b.z) >= 1.96 ? "#B83A4B" : "#FDA4AF");
                    return (
                      <g key={b.key}>
                        <rect x={x} y={y} width={32} height={Math.max(1.5, barH)} fill={fill} rx={1.5} />
                        <text x={x + 16} y={182} textAnchor="end" fontSize={8} fill="#1E293B" transform={`rotate(-28 ${x + 16} 182)`}>
                          {b.label}
                        </text>
                      </g>
                    );
                  })}
                </svg>
              ) : (
                <p className="mck-hint">Summary RBP ranking recorded ({activeRbp?.best_region}); select target events below to inspect locus-level cis-motifs.</p>
              )}
              {activeRbp?.pwm?.matrix && <PwmSequenceLogo matrix={activeRbp.pwm.matrix} motifId={activeRbp.pwm.motif_id} />}
              {!!data.bridge.length && (
                <div className="mck-bridge-strip">
                  <strong>SpliceNet × SpliceTrace FWER Module Couplings: </strong>
                  {data.bridge.slice(0, 6).map((b, idx) => (
                    <span key={idx} className="mck-tag">
                      {String(b.module || b.module_id || "")} → {String(b.rbp || b.RBP || "")} (P={Number(b.fwer_p ?? b.p_value ?? 0).toPrecision(2)})
                    </span>
                  ))}
                </div>
              )}
            </div>
          </div>

          {/* Bottom: Target Event & Proteome Consequence Drill-down */}
          <div className="mck-subcard" style={{marginTop: 16}}>
            <div className="mck-toolbar-row">
              <h4>3. Downstream Target Splicing Events &amp; UniProt / 50-nt EJC NMD Surgery</h4>
              <div className="mck-filter-pills">
                {([
                  ["all", "All |ΔPSI| ≥ 0.10"],
                  ["motif", `Carrying ${rbpName} Motif`],
                  ["domain", "UniProt Domain Altered"],
                  ["nmd", "50-nt EJC NMD Switch"],
                ] as const).map(([id, lbl]) => (
                  <button
                    key={id}
                    type="button"
                    className={`mck-chip ${consequenceFilter === id ? "active" : ""}`}
                    onClick={() => setConsequenceFilter(id)}
                  >
                    {lbl}
                  </button>
                ))}
                <input
                  aria-label="Filter target genes or domains"
                  placeholder="Search gene or protein domain…"
                  value={search}
                  onChange={e => setSearch(e.target.value)}
                  style={{padding: "5px 10px", fontSize: 12}}
                />
              </div>
            </div>

            <div className="report-table-wrap" style={{maxHeight: 300}}>
              <table>
                <thead>
                  <tr>
                    <th>Gene</th>
                    <th>Type</th>
                    <th className="num">ΔPSI</th>
                    <th>Cis-RBP Motifs</th>
                    <th>Frame / 50-nt EJC NMD</th>
                    <th>Remodeled UniProt / Pfam Domains</th>
                    <th>Dossier</th>
                  </tr>
                </thead>
                <tbody>
                  {filteredTargets.slice(0, 80).map(t => (
                    <tr key={t.event_id}>
                      <td><strong>{t.gene_symbol}</strong></td>
                      <td>{t.event_type}</td>
                      <td className="num">{t.dpsi != null ? `${t.dpsi >= 0 ? "+" : ""}${t.dpsi.toFixed(3)}` : "NA"}</td>
                      <td><small>{t.cis_rbps || "—"}</small></td>
                      <td>
                        <span className={`mck-tag ${t.nmd_switch_class.includes("NMD") ? "mck-tag-coral" : ""}`}>
                          {t.nmd_switch_class || t.frame_status || "coding"}
                        </span>
                      </td>
                      <td><small>{t.disrupted_protein_features || "—"}</small></td>
                      <td>
                        <button
                          type="button"
                          className="mck-link-btn"
                          onClick={() => onInspectGeneCard(t.gene_symbol, t.event_id)}
                        >
                          Open 4-Track GeneCard →
                        </button>
                      </td>
                    </tr>
                  ))}
                  {!filteredTargets.length && (
                    <tr><td colSpan={7}>No target events match the current filter criteria.</td></tr>
                  )}
                </tbody>
              </table>
            </div>
          </div>
        </>
      )}
    </section>
  );
}
