import type { WorkflowDefinition } from "./types";

export type DasEventSet = "stable_significant" | "significant";
export const DEFAULT_DAS_EVENT_SET: DasEventSet = "stable_significant";

export function isDasEventSet(value: unknown): value is DasEventSet {
  return value === "stable_significant" || value === "significant";
}

export function isDasSwitchable(workflow?: WorkflowDefinition): boolean {
  return workflow?.das_selection?.mode === "switchable";
}

/** Share display/run selection only across modules declaring that capability.
 * Historical artifact parameters are displayed separately and never rewritten.
 */
export function parametersWithDasSelection(
  workflow: WorkflowDefinition | undefined,
  parameters: Record<string, unknown>,
  eventSet: DasEventSet,
): Record<string, unknown> {
  return isDasSwitchable(workflow) ? {...parameters, event_set: eventSet} : parameters;
}
