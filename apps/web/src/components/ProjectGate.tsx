import { FormEvent, useState } from "react";
import { ArrowRight, FolderInput, FolderOpen, ShieldCheck } from "lucide-react";

type MigrationRequest = { legacy_config_path: string; output_root: string };
type Props = {
  defaultPath: string;
  loading: boolean;
  error?: string;
  onOpen: (path: string) => void;
  onMigrate?: (request: MigrationRequest) => void;
  onConfigureResources?: (request: { manifest_path: string; resource_bundle_path: string; output_root: string }) => void;
};

export function ProjectGate({ defaultPath, loading, error, onOpen, onMigrate, onConfigureResources }: Props) {
  const [mode, setMode] = useState<"open" | "migrate" | "resources">("open");
  const [path, setPath] = useState(defaultPath);
  const [legacyPath, setLegacyPath] = useState("");
  const [outputRoot, setOutputRoot] = useState("");
  const [bundlePath, setBundlePath] = useState("");

  function submit(event: FormEvent) {
    event.preventDefault();
    if (mode === "open" && path.trim()) onOpen(path.trim());
    if (mode === "migrate" && legacyPath.trim() && outputRoot.trim()) {
      onMigrate?.({ legacy_config_path: legacyPath.trim(), output_root: outputRoot.trim() });
    }
    if (mode === "resources" && path.trim() && bundlePath.trim() && outputRoot.trim()) {
      onConfigureResources?.({ manifest_path: path.trim(), resource_bundle_path: bundlePath.trim(), output_root: outputRoot.trim() });
    }
  }

  return (
    <main className="project-gate">
      <section className="gate-copy">
        <span className="eyebrow">Katana</span>
        <h1>Analyze RNA splicing.<br />Explain the results.</h1>
        <p>Explore differential splicing, gene expression, exon and transcript usage, GO/KEGG, RBP binding and co-splicing networks. Use bulk or completed single-cell results with sample metadata.</p>
        <div className="principle-row">
          <div><span>01</span><strong>Select an analysis</strong><small>Choose comparisons and set thresholds.</small></div>
          <div><span>02</span><strong>Inspect the results</strong><small>View plots, gene lists and source tables.</small></div>
          <div><span>03</span><strong>Export a report</strong><small>Share figures and Chinese/English explanations.</small></div>
        </div>
      </section>

      <section className="open-card">
        <div className="gate-mode-tabs">
          <button type="button" className={mode === "open" ? "active" : ""} onClick={() => setMode("open")}><FolderOpen size={15} />Existing project</button>
          <button type="button" className={mode === "migrate" ? "active" : ""} onClick={() => setMode("migrate")}><FolderInput size={15} />Import older project</button>
          <button type="button" className={mode === "resources" ? "active" : ""} onClick={() => setMode("resources")}>GO / MotifTrace resources</button>
        </div>
        <h2>{mode === "resources" ? "Attach annotation resources" : mode === "open" ? "Open your analysis project" : "Import an older project"}</h2>
        <p>{mode === "resources" ? "Attach GO/KEGG and RBP databases for the correct species to a new project copy. Your original project stays unchanged." : mode === "open" ? "Select the katana.project.json file containing your inputs, sample groups and comparisons." : "Import an older Katana configuration into a new project directory. Existing inputs and results stay untouched."}</p>
        <form onSubmit={submit}>
          {mode !== "migrate" ? (
            <>
              <label htmlFor="manifest-path">Project manifest path</label>
              <input id="manifest-path" value={path} onChange={(event) => setPath(event.target.value)} placeholder="Full path to katana.project.json" autoComplete="off" />
            </>
          ) : (
            <>
              <label htmlFor="legacy-path">Legacy .katana_project.yaml</label>
              <input id="legacy-path" value={legacyPath} onChange={(event) => setLegacyPath(event.target.value)} placeholder="Full path to .katana_project.yaml" autoComplete="off" />
              <label htmlFor="output-root">New external project directory</label>
              <input id="output-root" value={outputRoot} onChange={(event) => setOutputRoot(event.target.value)} placeholder="New project output directory" autoComplete="off" />
            </>
          )}
          {mode === "resources" && <>
            <label htmlFor="resource-bundle">Prepared resources.json path</label>
            <input id="resource-bundle" value={bundlePath} onChange={(event) => setBundlePath(event.target.value)} placeholder="Full path to resources.json" autoComplete="off" />
            <label htmlFor="resource-output">New empty project directory</label>
            <input id="resource-output" value={outputRoot} onChange={(event) => setOutputRoot(event.target.value)} placeholder="New project output directory" autoComplete="off" />
          </>}
          {error && <div className="form-error" role="alert">{error}</div>}
          <button className="primary-action" type="submit" disabled={loading || (mode === "resources" ? !path.trim() || !bundlePath.trim() || !outputRoot.trim() : mode === "open" ? !path.trim() : !legacyPath.trim() || !outputRoot.trim())}>
            {loading ? "Validating project..." : mode === "resources" ? "Attach resources and open" : mode === "open" ? "Open project" : "Migrate and open"}<ArrowRight size={17} />
          </button>
        </form>
        <div className="privacy-note"><ShieldCheck size={16} />Inputs stay local and analysis outputs never default to the source tree.</div>
      </section>
    </main>
  );
}
