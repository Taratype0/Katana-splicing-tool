import { Check, ChevronDown, RotateCcw, SlidersHorizontal } from "lucide-react";
import type { Artifact, ManifestContrast, Selection, WorkflowDefinition } from "../types";
import { displayParameterOption } from "../branding";
import { RbpLayerNote } from "./RbpEvidenceGuide";

type Props = {
  workflow: WorkflowDefinition;
  contrasts: ManifestContrast[];
  selectedContrasts: Set<string>;
  parameters: Record<string, unknown>;
  selections: Selection[];
  conditions: string[];
  artifacts?: Artifact[];
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
  artifacts = [],
  onToggleContrast,
  onSelectAllContrasts,
  onClearContrasts,
  onParameterChange,
  onReset,
}: Props) {
  const showComparisons = contrasts.length > 0 && workflow.supports_contrast_selection;
  const isReadSupport = workflow.workflow_id === "quality.stable_as";
  const visibleParameters = isReadSupport
    ? workflow.parameters.filter((parameter) => ["fdr", "min_abs_dpsi"].includes(parameter.key))
    : workflow.parameters;

  return (
    <section className="workflow-controls" aria-label="Workflow parameters">
      <RbpLayerNote workflowId={workflow.workflow_id} parameters={parameters} />
      <details style={{ margin: "0 0 16px", lineHeight: 1.6 }}>
        <summary>Statistical method and interpretation</summary>
        <p><strong>Method:</strong> {workflow.method}</p>
        <p><strong>Output:</strong> {workflow.output_summary}</p>
        <p><strong>Unit:</strong> {workflow.analysis_unit}</p>
        <p><strong>Limits:</strong> {workflow.limitations}</p>
        <p>Run executes only this module. Missing prerequisite evidence must be generated explicitly.</p>
      </details>
      <div className="control-heading">
        <div>
          <SlidersHorizontal size={16} />
          <strong>Analysis settings</strong>
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
              <small>{selectedContrasts.size} of {contrasts.length} included</small>
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

      {visibleParameters.length > 0 && (
        <div className="parameter-grid">
          {visibleParameters.map((parameter) => {
            const value = parameters[parameter.key] ?? parameter.default;
            if (workflow.workflow_id === "prioritize.rbp_mechanism") {
              const completedMaps = parameters.sequence_scope === "completed_sequence_maps";
              const consequences = parameters.sequence_scope === "event_consequences";
              if ((completedMaps || consequences) && ["window_bp", "motif_evidence", "site_pvalue", "binding_mode"].includes(parameter.key)) return null;
              if (!completedMaps && ["bundle_asset_id", "window_bundle_asset_id", "view_scope", "custom_motif_mode"].includes(parameter.key)) return null;
              if (!completedMaps && !consequences && parameter.key === "custom_motif_path") return null;
              if (!consequences && ["model_ids", "cluster_permutations"].includes(parameter.key)) return null;
            }
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
                      {parameter.options.map((option) => <option key={option} value={option}>{displayParameterOption(parameter.key, option)}</option>)}
                    </select>
                    <ChevronDown size={14} />
                  </span>
                ) : parameter.kind === "artifact" ? (
                  <span className="select-shell">
                    <select value={String(value)} onChange={(event) => onParameterChange(parameter.key, event.target.value)}>
                      <option value="">Latest compatible completed result</option>
                      {artifacts.filter((artifact) => parameter.options[0]?.split("|").includes(artifact.analysis_id) && artifact.name === parameter.options[1]).map((artifact) => (
                        <option key={artifact.artifact_id} value={artifact.artifact_id}>
                          {artifact.name} · {new Date(artifact.created_at).toLocaleString()} · {artifact.run_fingerprint.slice(0, 8)}
                        </option>
                      ))}
                    </select><ChevronDown size={14} />
                  </span>
                ) : parameter.kind === "selection" ? (
                  <span className="select-shell">
                    <select value={String(value)} onChange={(event) => onParameterChange(parameter.key, event.target.value)}>
                      <option value="">Choose a saved selection</option>
                      {selections.filter((selection) => selection.compatible_workflows ? selection.compatible_workflows.includes(workflow.workflow_id) : workflow.workflow_id === "validate.native_sashimi" ? selection.entity_type === "splice_event"
                        : workflow.workflow_id === "single_cell.sashimi" ? ["cell_splice_event", "donor_cell_type_event", "cell_type_splice_event"].includes(selection.entity_type)
                        : ["gene", "splice_event", "module_member"].includes(selection.entity_type)).map((selection) => (
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
                    value={typeof value === "object" ? JSON.stringify(value) : String(value)}
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
      {isReadSupport && <p className="preset-description" role="note" aria-label="AS read-support policy">Read support: IJC + SJC must be at least 10 in every selected sample across both groups. Missing counts fail this check. FDR and |dPSI| thresholds are applied separately.</p>}
    </section>
  );
}
