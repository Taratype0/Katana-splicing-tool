import contract from "../../../katana_next/science/rbp_evidence_layers.json";

export const RBP_EVIDENCE = contract;
export type RbpLayer = typeof contract.layers[number];

export function rbpEvidenceLayer(workflowId: string, parameters: Record<string, unknown> = {}): RbpLayer | undefined {
  const key = workflowId === "prioritize.rbp_mechanism"
    ? (contract.sequence_scope_layers as Record<string,string>)[String(parameters.sequence_scope ?? "boundary_association")]
    : (contract.workflow_layers as Record<string,string>)[workflowId];
  return contract.layers.find(layer => layer.id === key);
}
