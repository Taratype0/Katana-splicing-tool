import {useState} from "react";
import {useQuery} from "@tanstack/react-query";
import {api, type LiveOraResultPayload, type LiveOraRow} from "../api";

const FACET_META: Array<[string, string]> = [
  ["BP", "GO Biological Process"],
  ["CC", "GO Cellular Component"],
  ["MF", "GO Molecular Function"],
  ["KEGG", "KEGG Pathways"],
];

function FacetLollipopCard({
  code,
  title,
  rows,
}: {
  code: string;
  title: string;
  rows: LiveOraRow[];
}) {
  const sigCount = rows.filter(r => r.fdr <= 0.05).length;
  const maxVal = Math.max(2.2, ...rows.map(r => -Math.log10(Math.max(1e-30, r.fdr <= 0.05 ? r.fdr : r.pvalue)))) * 1.22;
  const fdrLine = -Math.log10(0.05);

  return (
    <div className="mck-subcard">
      <div className="mck-track-head">
        <strong>{title}</strong>
        <span className={`mck-tag ${sigCount > 0 ? "mck-tag-blue" : ""}`}>
          {sigCount > 0 ? `${sigCount} at BH FDR ≤ 0.05` : `0 at FDR ≤ 0.05 · Top nominal (${code})`}
        </span>
      </div>
      {!rows.length ? (
        <p className="mck-compact-note">0 overlapping genes in this branch for the current foreground selection.</p>
      ) : (
        <svg viewBox={`0 0 470 ${Math.max(115, rows.length * 28 + 34)}`} role="img" aria-label={title} style={{width: "100%", height: "auto"}}>
          {(() => {
            const padL = 225, padR = 54, padT = 12;
            const plotW = 470 - padL - padR;
            const xThr = padL + Math.min(1, fdrLine / maxVal) * plotW;
            return (
              <g>
                <line x1={padL} y1={padT - 4} x2={padL} y2={rows.length * 28 + 12} stroke="#CBD5E1" />
                {sigCount > 0 && (
                  <line x1={xThr} y1={padT - 4} x2={xThr} y2={rows.length * 28 + 12} stroke="#94A3B8" strokeDasharray="3,3" />
                )}
                {rows.map((r, idx) => {
                  const cy = padT + idx * 28 + 10;
                  const isSig = r.fdr <= 0.05;
                  const metricVal = -Math.log10(Math.max(1e-30, isSig ? r.fdr : r.pvalue));
                  const barEnd = padL + Math.max(0.03, Math.min(1, metricVal / maxVal)) * plotW;
                  const rad = Math.min(7.5, 3.6 + Math.sqrt(r.overlap_genes) * 0.9);
                  const lbl = r.term.length > 35 ? `${r.term.slice(0, 34)}…` : r.term;
                  return (
                    <g key={`${r.term}:${idx}`}>
                      <title>{`${r.term}\nOverlap: ${r.overlap_genes}/${r.term_size_in_bg} | P=${r.pvalue.toExponential(2)} | BH FDR=${r.fdr.toExponential(2)}\nGenes: ${r.hit_genes}`}</title>
                      <text x={padL - 8} y={cy + 3.5} textAnchor="end" fontSize={8.8} fill="#051C2C" fontWeight={isSig ? 600 : 400}>
                        {lbl}
                      </text>
                      <line
                        x1={padL}
                        y1={cy}
                        x2={barEnd}
                        y2={cy}
                        stroke={isSig ? "#2251FF" : "#94A3B8"}
                        strokeWidth={isSig ? 2.2 : 1.4}
                        strokeDasharray={isSig ? undefined : "3,2"}
                      />
                      <circle
                        cx={barEnd}
                        cy={cy}
                        r={rad}
                        fill={isSig ? (code === "KEGG" ? "#051C2C" : "#2251FF") : "#FFFFFF"}
                        stroke={code === "KEGG" ? "#051C2C" : "#0077C8"}
                        strokeWidth={1.5}
                      />
                      <text x={barEnd + rad + 5} y={cy + 3} fontSize={8.0} fill={isSig ? "#051C2C" : "#64748B"}>
                        n={r.overlap_genes} {isSig ? "" : `(P=${r.pvalue.toExponential(1)})`}
                      </text>
                    </g>
                  );
                })}
                <text x={padL + plotW / 2} y={rows.length * 28 + 26} textAnchor="middle" fontSize={8.5} fill="#64748B">
                  {sigCount > 0 ? "−log10(BH FDR) · Solid = FDR ≤ 0.05" : "−log10(Nominal P) · Exploratory hollow markers (FDR > 0.05)"}
                </text>
              </g>
            );
          })()}
        </svg>
      )}
    </div>
  );
}

export function LiveOraWorkbench({
  projectId,
  comparison,
  onInspectGeneCard,
}: {
  projectId: string;
  comparison: string;
  onInspectGeneCard?: (gene: string) => void;
}) {
  const [preset, setPreset] = useState("ps_das");
  const [fdrCutoff, setFdrCutoff] = useState(0.05);
  const [effectCutoff, setEffectCutoff] = useState(0.10);
  const [customText, setCustomText] = useState("");
  const [tableFilter, setTableFilter] = useState("ALL");

  const [committed, setCommitted] = useState({
    preset: "ps_das",
    fdr_cutoff: 0.05,
    effect_cutoff: 0.10,
    custom_genes: [] as string[],
  });

  const query = useQuery<LiveOraResultPayload>({
    queryKey: ["live-go-kegg-ora", projectId, comparison, committed],
    queryFn: () =>
      typeof api.computeLiveOra === "function"
        ? api.computeLiveOra(projectId, {comparison, ...committed})
        : Promise.resolve({available: false}),
  });

  if (typeof api.computeLiveOra !== "function") return null;

  function triggerRun(nextPreset?: string) {
    const usePreset = nextPreset ?? preset;
    if (nextPreset) setPreset(nextPreset);
    setCommitted({
      preset: usePreset,
      fdr_cutoff: fdrCutoff,
      effect_cutoff: effectCutoff,
      custom_genes: usePreset === "custom" ? [customText] : [],
    });
  }

  const data = query.data;
  const filteredRows = (data?.top_rows || []).filter(r => tableFilter === "ALL" || r.category_code === tableFilter);

  return (
    <section className="result-card mck-panel" aria-label="Live GO and KEGG Hypergeometric ORA Lab">
      <div className="mck-panel-header">
        <div>
          <span className="eyebrow">ZERO-CONFIG BUNDLED MOUSE &amp; HUMAN GO (BP / CC / MF) + KEGG ENGINE</span>
          <h3>Live 4-Facet GO &amp; KEGG Pathway Over-Representation Lab ({comparison || "Study"})</h3>
        </div>
        {data?.files && (
          <div className="mck-filter-pills">
            <a className="secondary-action compact" href={data.files.figure_png} target="_blank" rel="noreferrer">4-Facet PNG ↗</a>
            <a className="secondary-action compact" href={data.files.figure_svg} target="_blank" rel="noreferrer">Vector SVG ↗</a>
            <a className="secondary-action compact" href={data.files.results_tsv} download>Complete ORA TSV ↓</a>
          </div>
        )}
      </div>

      {/* 1-Click Preset Selector Pills */}
      <div className="mck-filter-pills" style={{marginBottom: 10}}>
        {([
          ["ps_das", "PS-DAS Splicing Genes"],
          ["deg_sig", "Significant DEG Genes"],
          ["splicing_only", "Splicing-Only Decoupled (Non-DEG)"],
          ["dual_deg_das", "Dual Regulated (DEG + PS-DAS)"],
          ["module:M1", "SpliceNet M1 Hubs"],
          ["custom", "Custom Gene Symbols"],
        ] as const).map(([id, lbl]) => (
          <button
            key={id}
            type="button"
            className={`mck-chip ${preset === id ? "active" : ""}`}
            onClick={() => {
              if (id === "custom") setPreset("custom");
              else triggerRun(id);
            }}
          >
            {lbl}
          </button>
        ))}
      </div>

      <div className="mck-controls-bar">
        <label>
          FDR Threshold
          <input
            aria-label="Live ORA FDR threshold"
            type="number"
            step={0.01}
            min={0.001}
            max={1.0}
            value={fdrCutoff}
            onChange={e => setFdrCutoff(Number(e.target.value) || 0.05)}
            style={{width: 90}}
          />
        </label>
        <label>
          Min |ΔPSI| / |log2FC|
          <input
            aria-label="Live ORA effect cutoff"
            type="number"
            step={0.05}
            min={0}
            max={4}
            value={effectCutoff}
            onChange={e => setEffectCutoff(Number(e.target.value) || 0.10)}
            style={{width: 95}}
          />
        </label>
        {preset === "custom" && (
          <label style={{flex: 1, minWidth: 280}}>
            Custom Gene Symbols (comma / space / newline separated)
            <input
              aria-label="Custom gene symbols for ORA"
              placeholder="e.g. Macf1, Ank2, Add3, Sptbn1, Ryr2, Camk2b…"
              value={customText}
              onChange={e => setCustomText(e.target.value)}
            />
          </label>
        )}
        <button type="button" className="primary-action compact" onClick={() => triggerRun()} disabled={query.isFetching}>
          {query.isFetching ? "Testing GO + KEGG…" : "Run Live GO + KEGG ORA"}
        </button>
      </div>

      {query.isLoading && <p className="mck-hint">Testing hypergeometric overlap across bundled GO (BP, CC, MF) &amp; KEGG catalogs against measured gene universe…</p>}
      {data && !data.available && <p className="mck-compact-note">{data.reason}</p>}

      {data?.available && data.facets && (
        <>
          <div className="mck-kpi-grid">
            <div className="mck-kpi-box">
              <small>Catalog Species</small>
              <strong style={{textTransform: "capitalize"}}>{data.species} GO + KEGG</strong>
              <span>Exact BH FDR per branch</span>
            </div>
            <div className="mck-kpi-box">
              <small>Foreground Genes Mapped</small>
              <strong>{data.foreground_count?.toLocaleString()}</strong>
              <span>{data.preset_label}</span>
            </div>
            <div className="mck-kpi-box">
              <small>Measured Background Universe</small>
              <strong>{data.background_count?.toLocaleString()}</strong>
              <span>Detected/testable workspace genes</span>
            </div>
            <div className="mck-kpi-box">
              <small>BH FDR ≤ 0.05 Terms</small>
              <strong>{data.significant_fdr05_count?.toLocaleString()}</strong>
              <span>Across BP · CC · MF · KEGG</span>
            </div>
          </div>

          {/* 2x2 Guaranteed 4-Facet Lollipop Grid (GO BP, GO CC, GO MF, KEGG) */}
          <div className="mck-split-grid">
            {FACET_META.map(([code, title]) => (
              <FacetLollipopCard key={code} code={code} title={title} rows={data.facets?.[code] || []} />
            ))}
          </div>

          {/* Interactive Traceable Hit-Gene Table */}
          <div className="mck-subcard" style={{marginTop: 14}}>
            <div className="mck-toolbar-row">
              <h4>Traceable GO &amp; KEGG Pathway Evidence Table</h4>
              <div className="mck-filter-pills">
                {["ALL", "BP", "CC", "MF", "KEGG"].map(cat => (
                  <button
                    key={cat}
                    type="button"
                    className={`mck-chip ${tableFilter === cat ? "active" : ""}`}
                    onClick={() => setTableFilter(cat)}
                  >
                    {cat === "ALL" ? "All Branches" : cat}
                  </button>
                ))}
              </div>
            </div>
            <div className="report-table-wrap" style={{maxHeight: 260}}>
              <table>
                <thead>
                  <tr>
                    <th>Ontology</th>
                    <th>Pathway / Term</th>
                    <th className="num">Hits</th>
                    <th className="num">Fold</th>
                    <th className="num">Nominal P</th>
                    <th className="num">BH FDR</th>
                    <th>Contributing Genes (Click for 4-Track GeneCard)</th>
                  </tr>
                </thead>
                <tbody>
                  {filteredRows.slice(0, 60).map((r, idx) => (
                    <tr key={`${r.category_code}:${r.term}:${idx}`}>
                      <td><span className={`mck-tag ${r.category_code === "KEGG" ? "mck-tag-blue" : ""}`}>{r.category_code}</span></td>
                      <td><strong>{r.term}</strong></td>
                      <td className="num">{r.overlap_genes} / {r.term_size_in_bg}</td>
                      <td className="num">{r.fold_enrichment.toFixed(2)}×</td>
                      <td className="num">{r.pvalue.toExponential(2)}</td>
                      <td className="num">{r.fdr.toExponential(2)}</td>
                      <td>
                        <div style={{display: "flex", flexWrap: "wrap", gap: 4}}>
                          {r.hit_genes.split(",").filter(Boolean).slice(0, 8).map(g => (
                            <button
                              key={g}
                              type="button"
                              className="mck-link-btn"
                              onClick={() => onInspectGeneCard?.(g)}
                            >
                              {g}
                            </button>
                          ))}
                        </div>
                      </td>
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
