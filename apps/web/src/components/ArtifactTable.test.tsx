import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { expect, it, vi } from "vitest";
import { api } from "../api";
import { ArtifactTable } from "./ArtifactTable";
import type { Artifact } from "../types";

vi.mock("../api", () => ({api: {artifactRows: vi.fn()}}));

it("passes exact selected source rows and the actual table browsing state to the target boundary", async () => {
  const row = {row_id: "B_vs_A::SE:chr1:100:200", event_id: "SE:chr1:100:200", event_type: "SE", gene_id: "G1", FDR: .17};
  vi.mocked(api.artifactRows).mockResolvedValue({columns: Object.keys(row), rows: [row], offset: 0, limit: 50, total: 55});
  const toggle = vi.fn(), browse = vi.fn();
  const query = new QueryClient({defaultOptions: {queries: {retry: false}}});
  const artifact = {artifact_id: "source", kind: "table"} as Artifact;
  render(<QueryClientProvider client={query}><ArtifactTable projectId="study" artifact={artifact} selectedIds={new Set()}
    onToggle={toggle} onFocus={vi.fn()} onBrowseSettingsChange={browse}/></QueryClientProvider>);
  fireEvent.click(await screen.findByRole("checkbox", {name: `Select ${row.row_id}`}));
  expect(toggle).toHaveBeenCalledWith(row.row_id, row);
  fireEvent.change(screen.getByRole("textbox", {name: "Search table"}), {target: {value: "G1"}});
  fireEvent.click(screen.getByRole("button", {name: "FDR"}));
  await waitFor(() => expect(browse).toHaveBeenLastCalledWith({search: "G1", sort_by: "FDR", descending: false, page_offset: 0, page_size: 50}));
  fireEvent.click(screen.getByRole("button", {name: "Next"}));
  await waitFor(() => expect(api.artifactRows).toHaveBeenLastCalledWith("study", "source", {
    offset: 50, limit: 50, search: "G1", sortBy: "FDR", descending: false}));
  expect(browse).toHaveBeenLastCalledWith({search: "G1", sort_by: "FDR", descending: false, page_offset: 50, page_size: 50});
});
