import {render, screen, fireEvent} from "@testing-library/react";
import {QueryClient, QueryClientProvider} from "@tanstack/react-query";
import {vi, test, expect} from "vitest";
import {ReportWorkspace} from "./ReportWorkspace";
import {ApplicationTree, applicationRoutes} from "./ApplicationTree";
import type {ProjectPayload, Artifact} from "../types";

vi.mock("./ReportFigures", () => ({
  reportGet: vi.fn(async () => ({
    panels: [
      {id: "p", title: "Sample PSI PCA", chapter: 2, contrast_id: "all"},
      {id: "m", title: "Native MAPP activity", chapter: 7, contrast_id: "C1"},
    ],
  })),
  reportUrl: (pid: string, rel: string) => `/api/projects/${pid}/report-source/${rel}`,
  ReportFigures: ({panels}: {panels: {title: string}[]}) => (
    <div>{panels.map(p => <span key={p.title}>{p.title}</span>)}</div>
  ),
}));
vi.mock("./ReportExplorer", () => ({
  ReportExplorer: ({artifacts, eventType}: {artifacts: Artifact[]; eventType?: string}) => (
    <div>
      {artifacts.map(a => <span key={a.artifact_id}>{a.name}</span>)}
      {eventType}
    </div>
  ),
}));
vi.mock("./ReportNetwork", () => ({ReportNetwork: () => <div>Inline interactive network</div>}));
vi.mock("./PsiLandscapeLab", () => ({PsiLandscapeLab: () => <div>Live PSI Landscape Lab</div>}));
vi.mock("./LiveOraWorkbench", () => ({LiveOraWorkbench: () => <div>Live GO KEGG ORA Lab</div>}));

const project = {
  snapshot: {project_id: "P", genome_build: "GRCm39"},
  manifest: {samples: [], contrasts: [{contrast_id: "C1"}, {contrast_id: "C2"}]},
} as unknown as ProjectPayload;

function show(route: string, artifacts: Artifact[] = []) {
  return render(
    <QueryClientProvider client={new QueryClient({defaultOptions: {queries: {retry: false}}})}>
      <ReportWorkspace
        project={project}
        route={applicationRoutes.find(r => r.id === route)!}
        artifacts={artifacts}
        onAdd={vi.fn()}
        onNavigate={vi.fn()}
        onSetup={vi.fn()}
      />
    </QueryClientProvider>
  );
}

test("ordered 01-09 Consulting navigation includes QC, PSI Lab, DAS subtypes, SpliceNet, SpliceTrace and GeneCard", () => {
  const select = vi.fn();
  render(<ApplicationTree active="das" onSelect={select} />);
  expect(screen.getByRole("button", {name: /03 · Differential AS/})).toHaveAttribute("aria-current", "page");
  fireEvent.click(screen.getByRole("button", {name: /Event-type subanalyses/}));
  fireEvent.click(screen.getByRole("button", {name: /03\.SE · SE Subtype Analysis/}));
  expect(select).toHaveBeenCalledWith("event-SE");
  expect(screen.getByRole("button", {name: /02 · PSI Landscape/})).toBeInTheDocument();
  expect(screen.getByRole("button", {name: /05 · GO & KEGG Live Enrichment Lab/})).toBeInTheDocument();
  expect(screen.getByRole("button", {name: /06 · Katana-SpliceNet/})).toBeInTheDocument();
  expect(screen.getByRole("button", {name: /07 · Katana-SpliceTrace/})).toBeInTheDocument();
  expect(screen.getByRole("button", {name: /08 · Target Scorer & 4-Track GeneCard/})).toBeInTheDocument();
});

test("navigation adapts between bulk and single-cell assay modes", () => {
  const {rerender} = render(<ApplicationTree active="das" assayMode="bulk" onSelect={vi.fn()} />);
  expect(screen.queryByRole("button", {name: /02\.2 Single-Cell Dual Manifold/})).not.toBeInTheDocument();
  expect(screen.getByRole("button", {name: /03\.2 Isoform/})).toBeInTheDocument();
  rerender(<ApplicationTree active="sc-manifold" assayMode="single-cell" onSelect={vi.fn()} />);
  expect(screen.getByRole("button", {name: /02\.2 Single-Cell Dual Manifold/})).toBeInTheDocument();
  expect(screen.queryByRole("button", {name: /03\.2 Isoform/})).not.toBeInTheDocument();
});

test("shared figures and Live PSI Lab mount on 02 patterns without comparison leakage", async () => {
  show("patterns");
  expect(screen.getByText("Live PSI Landscape Lab")).toBeInTheDocument();
  expect(await screen.findByText("Sample PSI PCA")).toBeInTheDocument();
  expect(screen.queryByText("Native MAPP activity")).not.toBeInTheDocument();
});

test("motif figures are visible directly alongside SpliceTrace Lab", async () => {
  show("motif");
  expect(await screen.findByText("Native MAPP activity")).toBeInTheDocument();
});

test("DAS route selects only its method and comparison", () => {
  show("das", [
    {artifact_id: "1", name: "C1 DAS", metadata: {source_analysis_id: "differential.as"}, contrast_ids: ["C1"]},
    {artifact_id: "2", name: "C2 DAS", metadata: {source_analysis_id: "differential.as"}, contrast_ids: ["C2"]},
    {artifact_id: "3", name: "C1 DEG", metadata: {source_analysis_id: "differential.deg"}, contrast_ids: ["C1"]},
  ] as unknown as Artifact[]);
  expect(screen.getByText("C1 DAS")).toBeInTheDocument();
  expect(screen.queryByText("C2 DAS")).not.toBeInTheDocument();
  expect(screen.queryByText("C1 DEG")).not.toBeInTheDocument();
  fireEvent.change(screen.getByLabelText("Analysis comparison"), {target: {value: "C2"}});
  expect(screen.getByText("C2 DAS")).toBeInTheDocument();
});

test("SpliceNet mounts directly without an iframe", () => {
  const {container} = show("network");
  expect(screen.getByText("Inline interactive network")).toBeInTheDocument();
  expect(container.querySelector("iframe")).toBeNull();
});
