import {FormEvent, useState} from "react";
import {ArrowRight, FolderInput, FolderOpen} from "lucide-react";
import {BulkInputForm} from "./BulkInputForm";
import {PathPicker} from "./PathPicker";

type MigrationRequest = {legacy_config_path: string; output_root: string};
type Props = {
  defaultPath: string;
  loading: boolean;
  error?: string;
  onImportReport?: (request: {report_path: string; output_root: string}) => void;
  onOpen: (path: string) => void;
  onMigrate?: (request: MigrationRequest) => void;
  onConfigureResources?: (request: {manifest_path: string; resource_bundle_path: string; output_root: string}) => void;
};

const WORKFLOW_MODES: Array<{
  title: string;
  badge: string;
  desc: string;
  targetMode: "new" | "report" | "open" | "resources";
}> = [
  {
    title: "Portable Report · Folder or REPORT.zip",
    badge: "INTERACTIVE WORKSPACE · BULK & SC",
    desc: "Explore completed reports with live GO+KEGG ORA, Beyond-Jutils PSI PCA/UMAP, and 4-Track GeneCards",
    targetMode: "report",
  },
  {
    title: "Standalone Bulk Inputs · rMATS + DEG",
    badge: "PYTHON PACKAGE · PS-DAS & SPLICENET",
    desc: "Register existing rMATS, DESeq2, DEXSeq, and MultiQC tables directly without rerunning alignment",
    targetMode: "new",
  },
  {
    title: "Single-Cell & Existing Projects",
    badge: "SMART-SEQ3 · EVENT-PSI UMAP",
    desc: "Open a katana.project.json workspace for bulk cohorts or Smart-seq3 donor-aware splicing manifolds",
    targetMode: "open",
  },
  {
    title: "Reference & Annotation Libraries",
    badge: "FASTA · GTF · 170-PWM & UNIPROT",
    desc: "Bind local genome FASTA, GTF, or custom pathway GMTs for SpliceTrace and Sashimi target follow-up",
    targetMode: "resources",
  },
];

export function ProjectGate({
  defaultPath,
  loading,
  error,
  onOpen,
  onImportReport,
  onMigrate,
  onConfigureResources,
}: Props) {
  const [mode, setMode] = useState<"open" | "migrate" | "resources" | "report" | "new">(onImportReport ? "report" : "open");
  const [path, setPath] = useState(defaultPath);
  const [reportPath, setReportPath] = useState("");
  const [legacyPath, setLegacyPath] = useState("");
  const [outputRoot, setOutputRoot] = useState("");
  const [bundlePath, setBundlePath] = useState("");

  function submit(event: FormEvent) {
    event.preventDefault();
    if (mode === "report" && reportPath.trim()) onImportReport?.({report_path: reportPath.trim(), output_root: ""});
    if (mode === "open" && path.trim()) onOpen(path.trim());
    if (mode === "migrate" && legacyPath.trim() && outputRoot.trim()) {
      onMigrate?.({legacy_config_path: legacyPath.trim(), output_root: outputRoot.trim()});
    }
    if (mode === "resources" && path.trim() && bundlePath.trim() && outputRoot.trim()) {
      onConfigureResources?.({manifest_path: path.trim(), resource_bundle_path: bundlePath.trim(), output_root: outputRoot.trim()});
    }
  }

  return (
    <main className="project-gate">
      <section className="gate-copy">
        <span className="eyebrow">KATANA TRANSCRIPTOMIC INTELLIGENCE · EXECUTIVE EDITION</span>
        <h1>Katana</h1>
        <p>
          Publication-grade Bulk &amp; Single-Cell RNA splicing discovery, Precision-Kernel Event-PSI UMAP,
          live GO/KEGG enrichment, 170-RBP PWM scanning, and 4-Track molecular GeneCards in one local workspace.
        </p>

        <div className="mck-quick-grid" aria-label="Supported analysis entrypoints">
          {WORKFLOW_MODES.map(item => (
            <button
              key={item.title}
              type="button"
              className="mck-quick-card"
              disabled={loading}
              onClick={() => setMode(item.targetMode)}
            >
              <span className="mck-quick-badge">{item.badge}</span>
              <strong>{item.title}</strong>
              <small>{item.desc}</small>
            </button>
          ))}
        </div>
      </section>

      <section className="open-card">
        <div className="gate-mode-tabs">
          <button type="button" className={mode === "new" ? "active" : ""} onClick={() => setMode("new")}>New bulk inputs</button>
          <button type="button" className={mode === "report" ? "active" : ""} onClick={() => setMode("report")}>Open report</button>
          <button type="button" className={mode === "open" ? "active" : ""} onClick={() => setMode("open")}><FolderOpen size={15} />Existing project</button>
          <button type="button" className={mode === "migrate" ? "active" : ""} onClick={() => setMode("migrate")}><FolderInput size={15} />Import older project</button>
          <button type="button" className={mode === "resources" ? "active" : ""} onClick={() => setMode("resources")}>GO / MotifTrace resources</button>
        </div>
        {mode === "new" ? (
          <BulkInputForm onCreated={onOpen} />
        ) : (
          <>
            <h2>
              {mode === "report"
                ? "Open a report folder or ZIP"
                : mode === "resources"
                ? "Attach annotation resources"
                : mode === "open"
                ? "Open your analysis project"
                : "Import an older project"}
            </h2>
            <form onSubmit={submit}>
              {mode === "report" && (
                <PathPicker
                  label="Report folder, ZIP or manifest"
                  value={reportPath}
                  onChange={setReportPath}
                  extensions={["report-manifest.json", ".zip"]}
                  placeholder="Select extracted report folder or REPORT.zip"
                />
              )}
              {(mode === "open" || mode === "resources") && (
                <PathPicker
                  label="Project manifest path"
                  kind="file"
                  extensions={[".json"]}
                  value={path}
                  onChange={setPath}
                  placeholder="Select a project"
                />
              )}
              {mode === "migrate" && (
                <>
                  <PathPicker
                    label="Legacy .katana_project.yaml"
                    kind="file"
                    extensions={[".yaml", ".yml"]}
                    value={legacyPath}
                    onChange={setLegacyPath}
                    placeholder="Select an older project"
                  />
                  <PathPicker
                    label="New external project directory"
                    kind="directory"
                    value={outputRoot}
                    onChange={setOutputRoot}
                    placeholder="Select an empty workspace folder"
                  />
                </>
              )}
              {mode === "resources" && (
                <>
                  <PathPicker
                    label="Prepared resources.json path"
                    kind="file"
                    extensions={[".json"]}
                    value={bundlePath}
                    onChange={setBundlePath}
                    placeholder="Select resources.json"
                  />
                  <PathPicker
                    label="New empty project directory"
                    kind="directory"
                    value={outputRoot}
                    onChange={setOutputRoot}
                    placeholder="Select an empty workspace folder"
                  />
                </>
              )}
              {error && <div className="form-error" role="alert">{error}</div>}
              <button
                className="primary-action"
                type="submit"
                disabled={
                  loading ||
                  (mode === "resources"
                    ? !path.trim() || !bundlePath.trim() || !outputRoot.trim()
                    : mode === "report"
                    ? !reportPath.trim()
                    : mode === "open"
                    ? !path.trim()
                    : !legacyPath.trim() || !outputRoot.trim())
                }
              >
                {loading
                  ? "Validating project..."
                  : mode === "resources"
                  ? "Attach resources and open"
                  : mode === "report"
                  ? "Import report and explore"
                  : mode === "open"
                  ? "Open project"
                  : "Migrate and open"}
                <ArrowRight size={17} />
              </button>
            </form>
          </>
        )}
      </section>
    </main>
  );
}
