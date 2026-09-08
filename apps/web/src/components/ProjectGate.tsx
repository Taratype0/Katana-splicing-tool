import { FormEvent, useState } from "react";
import { ArrowRight, FolderInput, FolderOpen, ShieldCheck } from "lucide-react";

type MigrationRequest = { legacy_config_path: string; output_root: string };
type Props = {
  defaultPath: string;
  loading: boolean;
  error?: string;
  onOpen: (path: string) => void;
  onMigrate?: (request: MigrationRequest) => void;
};

export function ProjectGate({ defaultPath, loading, error, onOpen, onMigrate }: Props) {
  const [mode, setMode] = useState<"open" | "migrate">("open");
  const [path, setPath] = useState(defaultPath);
  const [legacyPath, setLegacyPath] = useState("");
  const [outputRoot, setOutputRoot] = useState("");

  function submit(event: FormEvent) {
    event.preventDefault();
    if (mode === "open" && path.trim()) onOpen(path.trim());
    if (mode === "migrate" && legacyPath.trim() && outputRoot.trim()) {
      onMigrate?.({ legacy_config_path: legacyPath.trim(), output_root: outputRoot.trim() });
    }
  }

  return (
    <main className="project-gate">
      <section className="gate-copy">
        <span className="eyebrow">Katana 2</span>
        <h1>Splicing evidence,<br />kept in context.</h1>
        <p>Move from input audit to stable events, differential evidence, read validation, programs, single-cell analysis, and a reproducible report without hidden execution.</p>
        <div className="principle-row">
          <div><span>01</span><strong>Explicit runs</strong><small>Opening a project reads compatible cache only.</small></div>
          <div><span>02</span><strong>Traceable evidence</strong><small>Every result retains inputs, parameters, and parents.</small></div>
          <div><span>03</span><strong>Valid replicates</strong><small>Samples and donors remain inferential units.</small></div>
        </div>
      </section>

      <section className="open-card">
        <div className="gate-mode-tabs">
          <button type="button" className={mode === "open" ? "active" : ""} onClick={() => setMode("open")}><FolderOpen size={15} />Open v2 project</button>
          <button type="button" className={mode === "migrate" ? "active" : ""} onClick={() => setMode("migrate")}><FolderInput size={15} />Migrate legacy</button>
        </div>
        <h2>{mode === "open" ? "Open a manifest" : "Create a v2 project"}</h2>
        <p>{mode === "open" ? "The manifest fixes paths, references, samples, contrasts, and the external output root." : "Legacy inputs are discovered once and written to a new external run directory. Existing results stay untouched."}</p>
        <form onSubmit={submit}>
          {mode === "open" ? (
            <>
              <label htmlFor="manifest-path">Project manifest path</label>
              <input id="manifest-path" value={path} onChange={(event) => setPath(event.target.value)} placeholder="E:\\Analysis\\project\\katana.project.json" autoComplete="off" />
            </>
          ) : (
            <>
              <label htmlFor="legacy-path">Legacy .katana_project.yaml</label>
              <input id="legacy-path" value={legacyPath} onChange={(event) => setLegacyPath(event.target.value)} placeholder="E:\\Project\\results\\.katana_project.yaml" autoComplete="off" />
              <label htmlFor="output-root">New external project directory</label>
              <input id="output-root" value={outputRoot} onChange={(event) => setOutputRoot(event.target.value)} placeholder="E:\\KatanaRuns\\project-v2" autoComplete="off" />
            </>
          )}
          {error && <div className="form-error" role="alert">{error}</div>}
          <button className="primary-action" type="submit" disabled={loading || (mode === "open" ? !path.trim() : !legacyPath.trim() || !outputRoot.trim())}>
            {loading ? "Validating project..." : mode === "open" ? "Open project" : "Migrate and open"}<ArrowRight size={17} />
          </button>
        </form>
        <div className="privacy-note"><ShieldCheck size={16} />Inputs stay local and analysis outputs never default to the source tree.</div>
      </section>
    </main>
  );
}
