import { FormEvent, useState } from "react";
import { ArrowRight, FolderInput, FolderOpen } from "lucide-react";

type MigrationRequest = { legacy_config_path: string; output_root: string };
type Props = {
  defaultPath: string;
  loading: boolean;
  error?: string;
  onImportReport?: (request: {report_path:string;output_root:string}) => void;
  onOpen: (path: string) => void;
  onMigrate?: (request: MigrationRequest) => void;
  onConfigureResources?: (request: { manifest_path: string; resource_bundle_path: string; output_root: string }) => void;
};

export function ProjectGate({ defaultPath, loading, error, onOpen, onImportReport, onMigrate, onConfigureResources }: Props) {
  const [mode, setMode] = useState<"open" | "migrate" | "resources" | "report">("open");
  const [path, setPath] = useState(defaultPath);
  const [legacyPath, setLegacyPath] = useState("");
  const [outputRoot, setOutputRoot] = useState("");
  const [bundlePath, setBundlePath] = useState("");

  function submit(event: FormEvent) {
    event.preventDefault();
    if (mode === "report" && path.trim() && outputRoot.trim()) onImportReport?.({report_path:path.trim(),output_root:outputRoot.trim()});
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
      <section className="gate-copy"><h1>Katana</h1></section>

      <section className="open-card">
        <div className="gate-mode-tabs"><button type="button" className={mode === "report" ? "active" : ""} onClick={()=>setMode("report")}>Open extracted report</button>
          <button type="button" className={mode === "open" ? "active" : ""} onClick={() => setMode("open")}><FolderOpen size={15} />Existing project</button>
          <button type="button" className={mode === "migrate" ? "active" : ""} onClick={() => setMode("migrate")}><FolderInput size={15} />Import older project</button>
          <button type="button" className={mode === "resources" ? "active" : ""} onClick={() => setMode("resources")}>GO / MotifTrace resources</button>
        </div>
        <h2>{mode === "report" ? "Open an extracted report" : mode === "resources" ? "Attach annotation resources" : mode === "open" ? "Open your analysis project" : "Import an older project"}</h2>
        <form onSubmit={submit}>
          {mode === "report" && <><label>New report workspace directory<input value={outputRoot} onChange={e=>setOutputRoot(e.target.value)}/></label></>}
          {mode !== "migrate" ? (
            <>
              <label htmlFor="manifest-path">{mode === "report" ? "Report directory or report-manifest.json" : "Project manifest path"}</label>
              <input id="manifest-path" value={path} onChange={(event) => setPath(event.target.value)} placeholder={mode === "report" ? "Path to extracted report" : "Path to katana.project.json"} autoComplete="off" />
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
          <button className="primary-action" type="submit" disabled={loading || (mode === "resources" ? !path.trim() || !bundlePath.trim() || !outputRoot.trim() : (mode === "open" || mode === "report") ? !path.trim() || (mode === "report" && !outputRoot.trim()) : !legacyPath.trim() || !outputRoot.trim())}>
            {loading ? "Validating project..." : mode === "resources" ? "Attach resources and open" : mode === "report" ? "Import report and explore" : mode === "open" ? "Open project" : "Migrate and open"}<ArrowRight size={17} />
          </button>
        </form>

      </section>
    </main>
  );
}
