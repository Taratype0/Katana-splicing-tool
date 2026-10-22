"""One public mode catalog over the frozen MotifTrace workflows.

This file contains routing and product descriptions, never scientific defaults.
The existing workflow definitions remain the source of parameter values.
"""
from __future__ import annotations

from copy import deepcopy

PACKAGE_VERSION = "1.0.0.post1"
SCIENTIFIC_RELEASE = "1.0.0-frozen.20261017"

_MODES = (
    dict(mode_id="boundary", label="Boundary motifs and reference binding",
         workflow_id="prioritize.rbp_mechanism",
         fixed_parameters={"sequence_scope": "boundary_association"},
         description="Selected DAS boundary motifs and configured reference CLIP evidence; their statistics stay separate.",
         operation="new_analysis", requires_selection=True),
    dict(mode_id="se_impact", label="SE fixed-site sequence analysis",
         workflow_id="prioritize.se_rbp_impact", fixed_parameters={},
         description="Existing SE sequence-impact method with matched controls, native FIMO and split genomic coordinates.",
         operation="new_analysis", requires_selection=True),
    dict(mode_id="completed_maps", label="Completed sequence maps",
         workflow_id="prioritize.rbp_sequence_results", fixed_parameters={},
         description="Verify and browse completed SE evidence, short RNA and original model/window statistics; no native refit.",
         operation="completed_result_import", requires_selection=True),
    dict(mode_id="consequences", label="Five-event RNA consequences",
         workflow_id="prioritize.motiftrace_consequences", fixed_parameters={},
         description="SE, A3SS, A5SS, MXE and RI local RNA forms, motif changes and the existing paired-window/region tests.",
         operation="new_analysis", requires_selection=True),
    dict(mode_id="usage", label="Transcript and exon usage to motifs",
         workflow_id="prioritize.motiftrace_usage", fixed_parameters={},
         description="Reuse existing DTU, DEXSeq, DEG, abundance and sequence evidence; no implicit quantification or new single-case motif P.",
         operation="existing_result_integration", requires_selection=False),
    dict(mode_id="parallel", label="Completed MAPP and MotifTrace",
         workflow_id="prioritize.motiftrace_parallel", fixed_parameters={},
         description="Read a verified completed parallel-model bundle; retain original methods, test families and native P values separately.",
         operation="completed_result_import", requires_selection=False),
)


def mode_catalog() -> dict:
    """Return a detached, JSON-serializable catalog without importing a GUI."""
    return {"package_version": PACKAGE_VERSION, "scientific_release": SCIENTIFIC_RELEASE,
            "modes": deepcopy(list(_MODES))}


def resolve_mode(mode: str) -> dict:
    for item in _MODES:
        if item["mode_id"] == mode:
            return deepcopy(item)
    raise ValueError(f"Unknown MotifTrace mode {mode!r}. Choose: " + ", ".join(m["mode_id"] for m in _MODES))


def resolve_parameters(mode: str, parameters: dict | None = None) -> tuple[dict, dict]:
    spec = resolve_mode(mode)
    if parameters is not None and not isinstance(parameters, dict):
        raise ValueError("MotifTrace parameters must be a JSON object / Python dict")
    values = dict(parameters or {})
    for key, value in spec["fixed_parameters"].items():
        if key in values and values[key] != value:
            raise ValueError(f"Mode {mode!r} fixes {key}={value!r}; choose the corresponding mode instead")
        values[key] = value
    return spec, values
