import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, expect, it, vi } from "vitest";
import { api } from "../api";
import { TargetWorkbench } from "./TargetWorkbench";
import type { TargetBasketEntry } from "../types";

vi.mock("../api", () => ({api: {exportTargets: vi.fn(), planTargets: vi.fn(), runTargets: vi.fn(), targetRun: vi.fn(),
  importTargetResults: vi.fn(), targetFileUrl: (project: string, run: string, path: string) => `/api/${project}/${run}/${path}`}}));
const entries: TargetBasketEntry[] = [{basket_id: "basket1", source_artifact_id: "source1", source_name: "DAS results",
  row_ids: ["B_vs_A::exact-event"], target_kind: "event", labels: ["Gene1 SE exact-event"],
  filters: {source_parameters: {fdr: .05}, browse_settings: {fdr: .2, event_set: "stable_significant"}}}];
beforeEach(() => {vi.clearAllMocks(); localStorage.clear();});
function mount(items = entries) {
  const client = new QueryClient({defaultOptions: {queries: {retry: false}, mutations: {retry: false}}});
  render(<QueryClientProvider client={client}><TargetWorkbench projectId="study" projectFingerprint="snapshot"
    entries={items} contrasts={[{contrast_id: "B_vs_A", label: "B versus A", numerator: {label: "B", sample_ids: ["B1"]}, denominator: {label: "A", sample_ids: ["A1"]}},
      {contrast_id: "C_vs_A", label: "C versus A", numerator: {label: "C", sample_ids: ["C1"]}, denominator: {label: "A", sample_ids: ["A1"]}}]}
    onRemove={vi.fn()} onClear={vi.fn()}/></QueryClientProvider>);
}

it("exports exact target context, explicitly binds resources, runs and browses verified output", async () => {
  vi.mocked(api.exportTargets).mockResolvedValue({request_id: "request1", request_path: "/owned/request.json", archive_path: "/owned/request.zip",
    bundle_directory: "/owned/bundle", download_url: "/api/download/request1", target_count: 1, targets: [], resources: [{resource_id: "bam:sample1", kind: "bam", sample_id: "sample1"}]});
  vi.mocked(api.planTargets).mockResolvedValue({ready: true, request_id: "request1", errors: [],
    warnings: [{resource_id: "bam:sample2", message: "Missing verified BAM binding", affected_module: "sashimi"}], tasks: []});
  vi.mocked(api.runTargets).mockResolvedValue({run_id: "run1", request_id: "request1", status: "queued"});
  vi.mocked(api.targetRun).mockResolvedValue({run_id: "run1", request_id: "request1", status: "partial", report_url: "/verified/INDEX.html", archive_url: "/verified/REPORT.zip",
    result: {status: "partial", targets: [{target_id: "event1", status: "partial", report: "targets/event1/REPORT.html", error: "Sashimi unavailable: missing verified BAM"}]}});
  mount();
  expect(api.exportTargets).not.toHaveBeenCalled();
  expect(api.runTargets).not.toHaveBeenCalled();
  fireEvent.click(screen.getByRole("button", {name: "Export portable target request"}));
  expect(await screen.findByRole("link", {name: "Download request ZIP"})).toHaveAttribute("href", "/api/download/request1");
  expect(api.exportTargets).toHaveBeenCalledWith("study", {selections: [{source_artifact_id: "source1", row_ids: ["B_vs_A::exact-event"], target_kind: "event", filters: entries[0].filters}], modules: ["gene_card", "event_card", "sashimi"]});
  fireEvent.change(screen.getByLabelText("Resource path bam:sample1"), {target: {value: "/server/sample.bam"}});
  fireEvent.click(screen.getByRole("button", {name: "Check target readiness"}));
  await screen.findByText("Ready");
  expect(screen.getByText("bam:sample2: Missing verified BAM binding")).toBeInTheDocument();
  expect(api.planTargets).toHaveBeenCalledWith("study", {request_id: "request1", resources: {"bam:sample1": "/server/sample.bam"}});
  expect(api.runTargets).not.toHaveBeenCalled();
  fireEvent.click(screen.getByRole("button", {name: "Run requested targets"}));
  expect(await screen.findByRole("link", {name: "Open target index"})).toHaveAttribute("href", "/verified/INDEX.html");
  expect(screen.getByRole("link", {name: "Open report and data"})).toHaveAttribute("href", "/api/study/run1/targets/event1/REPORT.html");
  expect(screen.getByText(/Sashimi unavailable/)).toBeInTheDocument();
  expect(screen.getByText("Target run: partial")).toBeInTheDocument();
  expect(api.runTargets).toHaveBeenCalledTimes(1);
});

it("imports returned results only after explicit action and shows verified reports", async () => {
  vi.mocked(api.importTargetResults).mockResolvedValue({run_id: "imported1", request_id: "request1", status: "complete", imported: true,
    report_url: "/imported/INDEX.html", result: {status: "complete", targets: []}});
  vi.mocked(api.targetRun).mockResolvedValue({run_id: "imported1", request_id: "request1", status: "complete", report_url: "/imported/INDEX.html"});
  mount([]);
  expect(api.importTargetResults).not.toHaveBeenCalled();
  fireEvent.change(screen.getByLabelText("Result manifest or directory"), {target: {value: "/returned/result-manifest.json"}});
  fireEvent.click(screen.getByRole("button", {name: "Verify and import target results"}));
  await waitFor(() => expect(api.importTargetResults).toHaveBeenCalledWith("study", "/returned/result-manifest.json"));
  expect(await screen.findByRole("link", {name: "Open target index"})).toHaveAttribute("href", "/imported/INDEX.html");
  expect(api.runTargets).not.toHaveBeenCalled();
});

it("does not label a tampered imported result as complete", async () => {
  vi.mocked(api.importTargetResults).mockRejectedValue(new Error("Member hash mismatch"));
  mount([]);
  fireEvent.change(screen.getByLabelText("Result manifest or directory"), {target: {value: "/changed/result-manifest.json"}});
  fireEvent.click(screen.getByRole("button", {name: "Verify and import target results"}));
  expect(await screen.findByRole("alert")).toHaveTextContent("Member hash mismatch");
  expect(screen.queryByRole("link", {name: "Open target index"})).not.toBeInTheDocument();
});

it("exports explicit all-project comparison scope and displays the verified target scope", async () => {
  vi.mocked(api.exportTargets).mockResolvedValue({request_id: "all1", request_path: "/owned/request.json", archive_path: "/owned/request.zip",
    bundle_directory: "/owned/bundle", download_url: "/download/all1", target_count: 1,
    targets: [{target_id: "event1", target_kind: "event", event_id: "exact-event", comparison_ids: ["B_vs_A", "C_vs_A"]}], resources: []});
  mount();
  expect(screen.getByLabelText("Target comparison scope")).toHaveValue("source");
  fireEvent.change(screen.getByLabelText("Target comparison scope"), {target: {value: "all"}});
  fireEvent.click(screen.getByRole("button", {name: "Export portable target request"}));
  await screen.findByRole("link", {name: "Download request ZIP"});
  expect(api.exportTargets).toHaveBeenCalledWith("study", expect.objectContaining({contrast_ids: ["B_vs_A", "C_vs_A"]}));
  expect(screen.getByText("exact-event: B_vs_A, C_vs_A")).toBeInTheDocument();
});

it("requires a nonempty explicit custom comparison selection and never treats empty as all", async () => {
  vi.mocked(api.exportTargets).mockResolvedValue({request_id: "custom1", request_path: "/owned/request.json", archive_path: "/owned/request.zip",
    bundle_directory: "/owned/bundle", download_url: "/download/custom1", target_count: 1,
    targets: [{target_id: "event1", target_kind: "event", event_id: "exact-event", comparison_ids: ["B_vs_A"]}], resources: []});
  mount();
  fireEvent.change(screen.getByLabelText("Target comparison scope"), {target: {value: "custom"}});
  const exportButton = screen.getByRole("button", {name: "Export portable target request"});
  expect(exportButton).toBeDisabled();
  fireEvent.click(exportButton);
  expect(api.exportTargets).not.toHaveBeenCalled();
  fireEvent.click(screen.getByRole("checkbox", {name: "B versus A"}));
  expect(exportButton).toBeEnabled();
  fireEvent.click(exportButton);
  await screen.findByRole("link", {name: "Download request ZIP"});
  expect(api.exportTargets).toHaveBeenCalledWith("study", expect.objectContaining({contrast_ids: ["B_vs_A"]}));
  expect(screen.getByText("exact-event: B_vs_A")).toBeInTheDocument();
});
