import { act, renderHook } from "@testing-library/react";
import { beforeEach, expect, it } from "vitest";
import { browsingSnapshot, makeTargetEntry, mergeTargetEntries, useTargetBasket } from "./targetBasket";
import type { Artifact } from "./types";

const artifact: Artifact = {artifact_id: "original-source", name: "DAS complete", fingerprint: "content-sha", run_fingerprint: "native-fit",
  contrast_ids: ["B_vs_A"], parameters: {fdr: .05, count_mode: "JC"}, project_id: "study", analysis_id: "differential.as",
  kind: "table", entity_type: "splice_event", status: "complete", path: "saved-events.parquet", parent_artifact_ids: [],
  metadata: {}, warnings: [], created_at: "2026-10-23T00:00:00Z"};
const rows = [{row_id: "B_vs_A:event1", event_id: "SE:event1", event_type: "SE", gene_id: "G1", gene_symbol: "Gene1"},
  {row_id: "B_vs_A:event2", event_id: "SE:event2", event_type: "SE", gene_id: "G1", gene_symbol: "Gene1"}];
beforeEach(() => localStorage.clear());

it("keeps actual source parameters separate from user browsing thresholds and exact rows", () => {
  const settings = {fdr: .2, min_abs_dpsi: .03, event_set: "stable_significant", contrast_ids: ["B_vs_A"]};
  const snapshot = browsingSnapshot(artifact, settings, {search: "Gene1", sort_by: "FDR"});
  const entry = makeTargetEntry(artifact, rows, "event", snapshot);
  settings.fdr = .9;
  expect(entry.row_ids).toEqual(["B_vs_A:event1", "B_vs_A:event2"]);
  expect(entry.target_kind).toBe("event");
  expect(entry.filters.source_parameters).toEqual({fdr: .05, count_mode: "JC"});
  expect(entry.filters.browse_settings).toEqual({fdr: .2, min_abs_dpsi: .03, event_set: "stable_significant", contrast_ids: ["B_vs_A"], table: {search: "Gene1", sort_by: "FDR"}});
  expect(() => makeTargetEntry(artifact, [{row_id: "numeric-only", gene_symbol: "Ambiguous"}], "event", {})).toThrow(/exact event identities/);
});

it("merges duplicate rows only in the same kind, source and filtering context", () => {
  const first = makeTargetEntry(artifact, [rows[0]], "gene", {browse_settings: {fdr: .05}});
  const combined = mergeTargetEntries([first], makeTargetEntry(artifact, rows, "gene", first.filters));
  expect(combined).toHaveLength(1);
  expect(combined[0].row_ids).toEqual(rows.map(row => row.row_id));
  expect(combined[0].labels).toEqual(["Gene1", "Gene1"]);
  expect(mergeTargetEntries(combined, makeTargetEntry(artifact, rows, "event", first.filters))).toHaveLength(2);
  expect(mergeTargetEntries(combined, makeTargetEntry(artifact, rows, "gene", {browse_settings: {fdr: .2}}))).toHaveLength(2);
});

it("preserves targets across remounts and isolates projects and source-manifest versions", () => {
  const entry = makeTargetEntry(artifact, rows, "event", {});
  const first = renderHook(() => useTargetBasket("study1", "snapshot1"));
  act(() => first.result.current.add(entry));
  first.unmount();
  const reopened = renderHook(({project, fingerprint}) => useTargetBasket(project, fingerprint),
    {initialProps: {project: "study1", fingerprint: "snapshot1"}});
  expect(reopened.result.current.entries[0].row_ids).toEqual(entry.row_ids);
  reopened.rerender({project: "study2", fingerprint: "snapshot1"});
  expect(reopened.result.current.entries).toEqual([]);
  reopened.rerender({project: "study1", fingerprint: "snapshot2"});
  expect(reopened.result.current.entries).toEqual([]);
  reopened.rerender({project: "study1", fingerprint: "snapshot1"});
  expect(reopened.result.current.entries).toHaveLength(1);
  act(() => reopened.result.current.clear());
  expect(reopened.result.current.entries).toEqual([]);
});
