import { FormEvent, useState } from "react";
import { ArrowRight, FolderOpen, ShieldCheck } from "lucide-react";

type Props = {
  defaultPath: string;
  loading: boolean;
  error?: string;
  onOpen: (path: string) => void;
};

export function ProjectGate({ defaultPath, loading, error, onOpen }: Props) {
  const [path, setPath] = useState(defaultPath);

  function submit(event: FormEvent) {
    event.preventDefault();
    if (path.trim()) onOpen(path.trim());
  }

  return (
    <main className="project-gate">
      <section className="gate-copy">
        <span className="eyebrow">Katana Next / Platform alpha</span>
        <h1>Splicing evidence,<br />kept in context.</h1>
        <p>
          One local workspace connects data quality, statistical evidence, read validation,
          biological programs, and reproducible reporting.
        </p>
        <div className="principle-row">
          <div><span>01</span><strong>No hidden runs</strong><small>Opening a project only reads compatible cache.</small></div>
          <div><span>02</span><strong>One evidence trail</strong><small>Every plot resolves to source rows and parameters.</small></div>
          <div><span>03</span><strong>Valid units</strong><small>Samples and donors remain the statistical replicates.</small></div>
        </div>
      </section>

      <section className="open-card">
        <div className="open-card-icon"><FolderOpen size={25} /></div>
        <span className="eyebrow">Open a local project</span>
        <h2>Start from a manifest</h2>
        <p>The manifest fixes paths, genome build, samples, contrasts, and the external output root.</p>
        <form onSubmit={submit}>
          <label htmlFor="manifest-path">Project manifest path</label>
          <input
            id="manifest-path"
            value={path}
            onChange={(event) => setPath(event.target.value)}
            placeholder="E:\\Analysis\\project\\katana.project.json"
            autoComplete="off"
          />
          {error && <div className="form-error" role="alert">{error}</div>}
          <button className="primary-action" type="submit" disabled={loading || !path.trim()}>
            {loading ? "Auditing manifest..." : "Open project"}
            {!loading && <ArrowRight size={17} />}
          </button>
        </form>
        <div className="privacy-note"><ShieldCheck size={16} />No input is uploaded.</div>
      </section>
    </main>
  );
}
