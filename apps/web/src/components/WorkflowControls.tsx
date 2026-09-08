import { Check, ChevronDown, RotateCcw, SlidersHorizontal } from "lucide-react";
import type { ManifestContrast, Selection, WorkflowDefinition } from "../types";

type Props = {
  workflow: WorkflowDefinition;
  contrasts: ManifestContrast[];
  selectedContrasts: Set<string>;
  parameters: Record<string, unknown>;
  selections: Selection[];
  conditions: string[];
  onToggleContrast: (contrastId: string) => void;
  onSelectAllContrasts: () => void;
  onClearContrasts: () => void;
  onParameterChange: (key: string, value: unknown) => void;
  onReset: () => void;
};

export function WorkflowControls({
  workflow,
  contrasts,
  selectedContrasts,
  parameters,
  selections,
  conditions,
  onToggleContrast,
  onSelectAllContrasts,
  onClearContrasts,
  onParameterChange,
  onReset,
}: Props) {
  const showComparisons = contrasts.length > 0 && workflow.supports_contrast_selection;

  return (
    <section className="workflow-controls" aria-label="Workflow parameters">
      <div className="control-heading">
        <div>
          <SlidersHorizontal size={16} />
          <strong>Run configuration</strong>
        </div>
        <button type="button" className="text-action" onClick={onReset}>
          <RotateCcw size={14} />Reset defaults
        </button>
      </div>

      {showComparisons && (
        <div className="comparison-control">
          <div className="control-label-row">
            <div>
              <span>Comparisons</span>
              <small>{selectedContrasts.size || contrasts.length} of {contrasts.length} included</small>
            </div>
            <div>
              <button type="button" className="text-action" onClick={onSelectAllContrasts}>Select all</button>
              <button type="button" className="text-action" onClick={onClearContrasts}>Clear</button>
            </div>
          </div>
          <div className="comparison-chips">
            {contrasts.map((contrast) => {
              const checked = selectedContrasts.has(contrast.contrast_id);
              return (
                <button
                  type="button"
                  key={contrast.contrast_id}
                  className={checked ? "comparison-chip selected" : "comparison-chip"}
                  onClick={() => onToggleContrast(contrast.contrast_id)}
                  title={`${contrast.numerator.label} minus ${contrast.denominator.label}`}
                >
                  <span className="chip-check">{checked && <Check size={12} />}</span>
                  <span>{contrast.label}</span>
                </button>
              );
            })}
          </div>
        </div>
      )}

      {workflow.parameters.length > 0 && (
        <div className="parameter-grid">
          {workflow.parameters.map((parameter) => {
            const value = parameters[parameter.key] ?? parameter.default;
            return (
              <label className={parameter.kind === "boolean" ? "parameter-field toggle-field" : "parameter-field"} key={parameter.key}>
                <span>{parameter.label}</span>
                {parameter.kind === "boolean" ? (
                  <button
                    type="button"
                    role="switch"
                    aria-checked={Boolean(value)}
                    className={value ? "toggle on" : "toggle"}
                    onClick={() => onParameterChange(parameter.key, !Boolean(value))}
                  >
                    <span />
                  </button>
                ) : parameter.kind === "select" ? (
                  <span className="select-shell">
                    <select value={String(value)} onChange={(event) => onParameterChange(parameter.key, event.target.value)}>
                      {parameter.options.map((option) => <option key={option} value={option}>{option}</option>)}
                    </select>
                    <ChevronDown size={14} />
                  </span>
                ) : parameter.kind === "selection" ? (
                  <span className="select-shell">
                    <select value={String(value)} onChange={(event) => onParameterChange(parameter.key, event.target.value)}>
                      <option value="">Choose a saved selection</option>
                      {selections.map((selection) => (
                        <option key={selection.selection_id} value={selection.selection_id}>
                          {selection.name} ({selection.entity_ids.length})
                        </option>
                      ))}
                    </select>
                    <ChevronDown size={14} />
                  </span>
                ) : parameter.kind === "condition" ? (
                  <span className="select-shell">
                    <select value={String(value)} onChange={(event) => onParameterChange(parameter.key, event.target.value)}>
                      <option value="">Choose condition</option>
                      {conditions.map((condition) => <option key={condition} value={condition}>{condition}</option>)}
                    </select>
                    <ChevronDown size={14} />
                  </span>
                ) : (
                  <input
                    type={parameter.kind === "number" ? "number" : "text"}
                    value={String(value)}
                    min={parameter.minimum}
                    max={parameter.maximum}
                    step={parameter.step}
                    onChange={(event) => onParameterChange(
                      parameter.key,
                      parameter.kind === "number" ? Number(event.target.value) : event.target.value,
                    )}
                  />
                )}
                <small>{parameter.help}</small>
              </label>
            );
          })}
        </div>
      )}
    </section>
  );
}
