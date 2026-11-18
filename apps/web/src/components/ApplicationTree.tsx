import {useState} from "react";

export type ApplicationRoute = {
  id: string;
  label: string;
  workspace: string;
  match?: string;
  chapters?: number[];
  figure?: string;
  eventType?: string;
  parent?: string;
  workflow?: string;
  assayScope?: "all" | "bulk" | "single-cell";
};

const baseRoutes: ApplicationRoute[] = [
  // 01 · Study Design, QC & 5'-3' Coverage
  {id: "project", label: "01 · Study Design, QC & 5′–3′ Coverage", workspace: "project", chapters: [1], match: "project.audit|quality", figure: "coverage|sample|quality|sequencing read|donor", assayScope: "all"},
  {id: "samples", label: "01.1 Samples & Replicate QC", workspace: "quality", parent: "project", match: "project.audit|quality", chapters: [1], figure: "sample|quality|sequencing read|coverage|donor", assayScope: "all"},
  {id: "comparisons", label: "01.2 Comparisons & Pairing", workspace: "project", parent: "project", match: "project.audit", assayScope: "all"},
  {id: "inputs", label: "01.3 Inputs & Resource Setup", workspace: "project", parent: "project", assayScope: "all"},

  // 02 · PSI Landscape, Z-Heatmap, PCA & UMAP Lab (100% contiguous numbering in both Bulk & Single-Cell!)
  {id: "patterns", label: "02 · PSI Landscape, Z-Heatmap, PCA & UMAP Lab", workspace: "landscape", match: "landscape|quality", chapters: [2], figure: "heatmap|PCA|UMAP|similarity|PSI|manifold|clustermap", assayScope: "all"},
  {id: "overview", label: "02.1 Subtype Composition & PS-DAS Qualification", workspace: "landscape", parent: "patterns", chapters: [2], match: "landscape|quality", figure: "composition|qualification|attrition|raw-versus", assayScope: "all"},
  {id: "sc-manifold", label: "02.2 Single-Cell Dual Manifold & Subtypes", workspace: "landscape", parent: "patterns", chapters: [2], match: "single_cell", figure: "manifold|embedding|UMAP", assayScope: "single-cell"},

  // 03 · Differential Suite (PS-DAS, DEG & Isoforms — strictly sequential 03 -> 03.1 -> 03.2)
  {id: "das", label: "03 · Differential AS · PS-DAS Volcano", workspace: "differential", match: "^differential.as$", figure: "DAS|Differential AS volcano", assayScope: "all"},
  ...["SE", "RI", "A3SS", "A5SS", "MXE"].map(eventType => ({
    id: `event-${eventType}`,
    label: `03.${eventType} · ${eventType} Subtype Analysis`,
    workspace: "differential",
    match: "differential.as",
    parent: "das",
    eventType,
    assayScope: "all" as const,
  })),
  {id: "deg", label: "03.1 DEG · Interactive Volcano", workspace: "differential", match: "^differential.deg$", figure: "DEG volcano", assayScope: "all"},
  {id: "isoform", label: "03.2 Isoform · DTU / DEXSeq", workspace: "differential", match: "dtu|dexseq|isoform|expression_response", chapters: [3, 4], figure: "DTU|DEXSeq|transcript|exon|isoform", assayScope: "bulk"},

  // 04 · Cross-Layer (DEG × PS-DAS) & Rescue Logic
  {id: "integrated", label: "04 · Cross-Layer (DEG × PS-DAS) & Rescue", workspace: "prioritize", workflow: "interpretation.joint_response", match: "evidence|joint|cross_comparison", chapters: [4], figure: "integration|candidate|evidence|overlap|decoupling|quadrant", assayScope: "all"},
  {id: "cross", label: "04.1 Multi-Contrast Rescue & Concordance", workspace: "differential", parent: "integrated", match: "cross_comparison|interaction", chapters: [4], figure: "overlap|shared|specific|similarity|interaction|cross|decoupling", assayScope: "all"},

  // 05 · GO & KEGG Live Enrichment Lab
  {id: "enrichment", label: "05 · GO & KEGG Live Enrichment Lab", workspace: "prioritize", match: "enrichment|pathways", chapters: [5], figure: "GO|KEGG|pathway|enrichment", assayScope: "all"},

  // 06 · Katana-SpliceNet · Co-Splicing Gene Network
  {id: "network", label: "06 · Katana-SpliceNet · Gene Network", workspace: "programs", match: "programs|splicenet", chapters: [6], assayScope: "all"},

  // 07 · Katana-SpliceTrace · RBP, Motif & Proteome Lab
  {id: "motif", label: "07 · Katana-SpliceTrace · RBP & Motif Lab", workspace: "prioritize", match: "motif|mapp|rbp|binding|splicetrace|proteome", chapters: [7], figure: "SpliceTrace|11-channel|logo|NMD|domain|MAPP|CLIP|RBP|bridge", assayScope: "all"},

  // 08 · Target Scorer & 4-Track GeneCard / Sashimi
  {id: "targets", label: "08 · Target Scorer & 4-Track GeneCard", workspace: "targets", assayScope: "all"},
  {id: "candidates", label: "08.1 Multi-Evidence Scorecard", workspace: "prioritize", parent: "targets", match: "evidence|cross_comparison|integrative", chapters: [8], figure: "scorecard|candidate|evidence", assayScope: "all"},
  {id: "sashimi", label: "08.2 Server / Local BAM Sashimi", workspace: "validate", parent: "targets", match: "sashimi", figure: "sashimi", assayScope: "all"},

  // 09 · Standalone Report & Complete Exports
  {id: "reports", label: "09 · Standalone Report & Complete Exports", workspace: "reports", assayScope: "all"},
];

const children: Record<string, [string, string, string][]> = {
  motif: [
    ["sequence", "07.1 Positional 11-Channel & PWM Logos", "motif|sequence|11-channel|logo"],
    ["clip", "07.2 50-nt EJC NMD & UniProt Domains", "CLIP|RBP|NMD|domain|surgery"],
    ["cross", "07.3 Module × RBP Regulatory Bridge", "Cross RBP|candidate membership|bridge"],
    ["mapp", "07.4 Native MAPP Activities", "MAPP"],
  ],
  isoform: [
    ["dtu", "03.2a Transcript Usage (DTU)", "DTU|transcript"],
    ["dexseq", "03.2b Exon Usage (DEXSeq)", "DEXSeq|exon"],
  ],
};

export const applicationRoutes: ApplicationRoute[] = baseRoutes.flatMap(r => [
  r,
  ...(children[r.id] || []).map(([id, label, figure]) => ({
    id: `${r.id}-${id}`,
    label,
    workspace: r.workspace,
    match: r.match,
    chapters: r.chapters,
    figure,
    parent: r.id,
    assayScope: r.assayScope ?? ("all" as const),
  })),
]);

export function ApplicationTree({
  active,
  onSelect,
  assayMode = "bulk",
}: {
  active: string;
  onSelect: (id: string) => void;
  assayMode?: string;
}) {
  const [open, setOpen] = useState<Record<string, boolean>>({});
  const current = applicationRoutes.find(r => r.id === active);
  const visibleRoutes = applicationRoutes.filter(
    r => !r.assayScope || r.assayScope === "all" || r.assayScope === assayMode
  );
  const subs = (id: string) => visibleRoutes.filter(c => c.parent === id);

  return (
    <aside className="workspace-nav" aria-label="Analyses">
      <div className="brand-lockup">
        <img className="brand-mark" src="/katana_icon.png" alt="Katana logo" width={40} height={40} />
        <div>
          <strong>KATANA</strong>
          <small>{assayMode === "single-cell" ? "SINGLE-CELL SPLICING SUITE" : "BULK SPLICING & ISOFORM SUITE"}</small>
        </div>
      </div>
      <nav>
        {visibleRoutes
          .filter(r => !r.parent)
          .map(r => {
            const items = subs(r.id);
            const expanded = open[r.id] ?? current?.parent === r.id;
            return (
              <div key={r.id} className="nav-item-group">
                <button
                  type="button"
                  className={`nav-item chapter-root ${active === r.id ? "active" : ""}`}
                  aria-current={active === r.id ? "page" : undefined}
                  onClick={() => onSelect(r.id)}
                >
                  <span>{r.label}</span>
                </button>
                {!!items.length && (
                  <>
                    <button
                      type="button"
                      className="tree-expand"
                      aria-expanded={expanded}
                      onClick={() => setOpen(s => ({...s, [r.id]: !expanded}))}
                    >
                      {expanded ? "▾" : "▸"} {r.id === "das" ? "Event-type subanalyses" : "Subanalyses"}
                    </button>
                    {expanded &&
                      items.map(c => (
                        <button
                          type="button"
                          key={c.id}
                          className={`nav-item nested ${active === c.id ? "active" : ""}`}
                          aria-current={active === c.id ? "page" : undefined}
                          onClick={() => onSelect(c.id)}
                        >
                          <span>{c.label}</span>
                        </button>
                      ))}
                  </>
                )}
              </div>
            );
          })}
      </nav>
    </aside>
  );
}
