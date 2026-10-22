import { fireEvent, render, screen } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import type { WorkflowDefinition } from "../types";
import { DasSelectionNote } from "./DasSelectionNote";
import { WorkflowControls } from "./WorkflowControls";

function workflow(mode: NonNullable<WorkflowDefinition["das_selection"]>["mode"]): WorkflowDefinition {
  return {workflow_id: "differential.as", workspace_id: "differential", label: "DAS", question: "DAS",
    method: "Completed native caller", output_summary: "DAS", analysis_unit: "event", expected_artifacts: 1,
    required_asset_kinds: [], optional_asset_kinds: [], supports_contrast_selection: false,
    enabled: true, limitations: "Unchanged caller statistics", das_selection: {
      mode, label: mode === "rs_required" ? "RS-DAS required" : "DAS event selection", description: "Recorded selection contract", default_event_set: "stable_significant",
    }, parameters: mode === "switchable" ? [{key: "event_set", label: "DAS event set", kind: "select", default: "stable_significant",
      help: "Choose the foreground", options: ["stable_significant", "significant"]}] : []};
}

it("shows readable DAS choices and sends the exact canonical value", () => {
  const change = vi.fn();
  render(<WorkflowControls workflow={workflow("switchable")} parameters={{}} contrasts={[]} selectedContrasts={new Set()}
    selections={[]} conditions={[]} onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()}
    onReset={vi.fn()} onParameterChange={change} />);
  expect(screen.getByRole("note", {name: "DAS event selection policy"})).toHaveTextContent("RS-DAS: significant + read support (default)");
  const select = screen.getByRole("combobox", {name: /DAS event set/});
  expect(select).toHaveValue("stable_significant");
  fireEvent.change(select, {target: {value: "significant"}});
  expect(change).toHaveBeenCalledWith("event_set", "significant");
  expect(screen.getByRole("option", {name: "Ordinary DAS: significant events"})).toBeInTheDocument();
});

it("plainly labels RS-required modules without an ordinary-DAS control", () => {
  render(<DasSelectionNote workflow={workflow("rs_required")} parameters={{event_set: "stable_significant"}} />);
  const note = screen.getByRole("note", {name: "DAS event selection policy"});
  expect(note).toHaveTextContent("RS-DAS required.");
  expect(note).toHaveTextContent("does not change this module's eligible input");
  expect(screen.queryByRole("combobox")).not.toBeInTheDocument();
});

it("keeps network measurement populations separate from significant foregrounds", () => {
  render(<DasSelectionNote workflow={workflow("measurement_population")} parameters={{}} />);
  expect(screen.getByRole("note")).toHaveTextContent("The significant DAS / RS-DAS foreground does not replace this module's declared measured-input population.");
});

it("labels fixed completed report imports truthfully", () => {
  render(<DasSelectionNote workflow={{...workflow("not_applicable"), workflow_id: "reports.interpreted_bundle"}} parameters={{}} />);
  expect(screen.getByRole("note")).toHaveTextContent("Completed result scope.");
  expect(screen.getByRole("note")).toHaveTextContent("not recalculated");
});
