import {fireEvent, render, screen, within} from "@testing-library/react";
import {describe, expect, it, vi} from "vitest";
import {RbpEvidenceGuide, RbpLayerNote} from "./RbpEvidenceGuide";
import {RBP_EVIDENCE, rbpEvidenceLayer} from "../rbpEvidence";
import type {Artifact, WorkflowDefinition} from "../types";

const workflows = RBP_EVIDENCE.layers.map(layer => ({workflow_id:layer.workflow_id,enabled:true})) as WorkflowDefinition[];
const saved: Artifact = {
  artifact_id:"fixture",project_id:"fixture",analysis_id:"prioritize.rbp_mechanism",
  name:"Saved fixture",kind:"file",entity_type:"evidence",status:"complete",fingerprint:"same",
  run_fingerprint:"same",path:"not-read.tsv",parent_artifact_ids:[],contrast_ids:["PBS"],
  parameters:{},metadata:{},warnings:[],created_at:"2026-10-17T00:00:00Z",
};
function guide(artifacts:Artifact[] = []) {
  const select=vi.fn();
  render(<RbpEvidenceGuide workflows={workflows} artifacts={artifacts}
    activeWorkflowId="prioritize.rbp_mechanism" parameters={{sequence_scope:"boundary_association"}} onSelect={select}/>);
  return select;
}

describe("Three complementary RBP evidence layers", () => {
  it("keeps identity and upstream clues visible alongside MotifTrace", () => {
    guide();
    for(const layer of RBP_EVIDENCE.layers) {
      expect(screen.getByRole("heading", {name:`${layer.number}. ${layer.title_en}`})).toBeInTheDocument();
      expect(screen.getByText(layer.title_zh)).toBeInTheDocument();
      expect(screen.getByText(layer.boundary_en)).toBeInTheDocument();
    }
    expect(screen.getAllByText(/No saved output loaded/)).toHaveLength(3);
  });
  it("selects only an existing route and never dispatches a run", () => {
    const select=guide();
    expect(select).not.toHaveBeenCalled();
    fireEvent.click(screen.getByRole("button",{name:"View 1. RBP identity annotations"}));
    expect(select).toHaveBeenLastCalledWith("prioritize.rbp_support",undefined);
    fireEvent.click(screen.getByRole("button",{name:"View 2. Upstream motif / CLIP regulatory clues"}));
    expect(select).toHaveBeenLastCalledWith("prioritize.rbp_mechanism",{sequence_scope:"boundary_association"});
    fireEvent.click(screen.getByRole("button",{name:"View 3. MotifTrace RNA sequence consequences"}));
    expect(select).toHaveBeenLastCalledWith("prioritize.motiftrace_consequences",undefined);
    expect(screen.queryByRole("button",{name:/Run analysis/i})).not.toBeInTheDocument();
  });
  it("classifies saved mode not merely the main workflow name", () => {
    const artifacts = [
      {...saved,parameters:{sequence_scope:"boundary_association"}},
      {...saved,artifact_id:"second",parameters:{sequence_scope:"completed_parallel_models"}},
    ];
    guide(artifacts);
    expect(screen.getAllByText("1 saved outputs; not gene/event counts.")).toHaveLength(2);
    expect(artifacts[1].parameters.sequence_scope).toBe("completed_parallel_models");
  });
  it("shows locus-specific boundaries even for a reused MAPP bundle", () => {
    render(<RbpLayerNote workflowId="prioritize.rbp_mechanism" parameters={{sequence_scope:"completed_parallel_models"}}/>);
    const note=screen.getByRole("complementary",{name:"Current RBP evidence layer"});
    expect(within(note).getByText(/cohort motif-window P is not an event\/site P/)).toBeInTheDocument();
  });
  it("does not invent an RBP layer for unrelated workflows or unknown modes", () => {
    expect(rbpEvidenceLayer("programs.cosplice")).toBeUndefined();
    expect(rbpEvidenceLayer("prioritize.rbp_mechanism",{sequence_scope:"future"})).toBeUndefined();
  });
});
