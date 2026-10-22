import { expect, it } from "vitest";
import { DEFAULT_DAS_EVENT_SET, parametersWithDasSelection } from "./dasSelection";
import type { WorkflowDefinition } from "./types";

function workflow(mode: NonNullable<WorkflowDefinition["das_selection"]>["mode"]): WorkflowDefinition {
  return {das_selection: {mode, label: mode, description: mode, default_event_set: null}} as WorkflowDefinition;
}

it("shares the chosen foreground only across explicitly switchable workflows", () => {
  expect(DEFAULT_DAS_EVENT_SET).toBe("stable_significant");
  const saved = {quality_artifact_id: "qc-original", event_set: "stable_significant"};
  expect(parametersWithDasSelection(workflow("switchable"), saved, "significant"))
    .toEqual({quality_artifact_id: "qc-original", event_set: "significant"});
  expect(saved.event_set).toBe("stable_significant");
  for (const mode of ["rs_required", "measurement_population", "not_applicable"] as const) {
    const own = {sequence_scope: "event_consequences"};
    expect(parametersWithDasSelection(workflow(mode), own, "significant")).toBe(own);
    expect(own).not.toHaveProperty("event_set");
  }
});
