"""Validate explicitly selected contrasts and source-specific design boundaries.

This is sample-metadata preflight, not fitting or a replacement for native model
validation. No formula is evaluated and no comparison is generated from labels.
"""
from __future__ import annotations

from collections import Counter
from importlib.resources import files
import json
from pathlib import Path
import re

from .sample_identity import _META_REQUIRED, _read_csv, _rows_shape, validate_sample_rows

CONTRAST_PLAN_VERSION = "1.0.0"
_SAFE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_CONDITION = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_PAIR_KEYS = {"donor_id", "biological_unit_id"}
_REQUIRED = {"contrast", "treatment", "control"}
_OPTIONAL = {"design", "formula", "paired_by", "paired", "cell_type"}


def method_design_support():
    """Return a new copy of the source-audited, versioned capability table."""
    return json.loads(files("katana_next.contracts").joinpath("data/method-design-support-v1.json").read_text(encoding="utf-8"))


def _report(route, design, intent, assay_mode):
    return {"contrast_plan_version": CONTRAST_PLAN_VERSION, "valid": False,
            "scope": "selected_contrast_sample_metadata_preflight", "route": route,
            "requested_design": design, "intent": intent, "assay_mode": assay_mode,
            "execution_authorized": False, "scientific_tools_executed": 0,
            "statistical_models_fitted": 0, "all_pairs_generated": False,
            "contrast_count": 0, "selected_contrasts": [], "issues": [],
            "native_model_and_measurement_eligibility_checked": False,
            "runtime_ready": False}


def _issue(report, code, path, message, severity="error"):
    report["issues"].append({"code": code, "path": path, "message": message, "severity": severity})


def _finish(report):
    report["valid"] = not any(issue["severity"] == "error" for issue in report["issues"])
    report["contrast_count"] = len(report["selected_contrasts"])
    return report


def _terms(formula, report, path):
    if not isinstance(formula, str) or not formula.strip().startswith("~"):
        _issue(report, "DESIGN_FORMULA", path, "Declare a one-sided additive formula explicitly; no filename inference or formula evaluation.")
        return None
    rhs = formula.strip()[1:].strip()
    if any(x in rhs for x in [":", "*", "(", ")", "/", "|", "^", "-", "~"]):
        _issue(report, "METHOD_SPECIFIC_DESIGN_REQUIRED", path, "Interaction, numeric transformations, nested/random or non-intercept models need their specific formal adapter; they are not silently reduced to ~ condition.")
        return None
    terms = [term.strip() for term in rhs.split("+")]
    if not terms or any(not term for term in terms) or len(terms) != len(set(terms)) or "condition" not in terms:
        _issue(report, "DESIGN_TERMS", path, "Declare condition exactly once and distinct additive covariates; no automatic term removal.")
        return None
    if set(terms) - {"condition", "batch", "subject"}:
        _issue(report, "UNSUPPORTED_COVARIATE", path, "Current source adapter qualifies categorical condition, batch and subject only; arbitrary/numeric covariates need an explicit method adapter.")
        return None
    return terms


def _design_matrix(selected, terms, treatment, pair_by, report, path):
    # Lazy: T09/T10 metadata checks retain their small import/dependency boundary.
    import numpy as np
    columns, names = [np.ones(len(selected))], ["Intercept"]
    condition_column = None
    usable = True
    for term in terms:
        if term == "condition":
            condition_column = len(columns)
            columns.append(np.asarray([s["group"] == treatment for s in selected], dtype=float))
            names.append("condition:treatment_minus_control")
            continue
        key = pair_by if term == "subject" else "batch"
        values = [s.get(key) if key else None for s in selected]
        if any(not v for v in values):
            _issue(report, "MISSING_DESIGN_METADATA", path + "/" + term, "Every selected sample needs explicitly declared " + term + "; missing values are not imputed.")
            usable = False
            continue
        levels = sorted(set(values))
        if len(levels) < 2:
            _issue(report, "SINGLE_LEVEL_DESIGN_TERM", path + "/" + term, "Selected covariate has only one level; choose a valid reduced model explicitly.")
            usable = False
        for level in levels[1:]:
            columns.append(np.asarray([v == level for v in values], dtype=float))
            names.append(term + ":" + level)
    if not usable or condition_column is None:
        return {"matrix_columns": names, "rank": None, "row_count": len(selected),
                "column_count": len(columns), "residual_df": None, "contrast_estimable": None}
    matrix = np.column_stack(columns)
    rank = int(np.linalg.matrix_rank(matrix))
    coefficient = np.zeros(matrix.shape[1])
    coefficient[condition_column] = 1.0
    estimable = bool(np.allclose(coefficient, coefficient @ np.linalg.pinv(matrix) @ matrix, atol=1e-9, rtol=1e-9))
    residual = len(selected) - rank
    if rank != matrix.shape[1]:
        _issue(report, "DESIGN_RANK_DEFICIENT", path, "Selected design is confounded/rank deficient; no covariate is silently removed.")
    if not estimable:
        _issue(report, "CONTRAST_NOT_ESTIMABLE", path, "Requested effect is outside the sample design row space.")
    if residual <= 0:
        _issue(report, "NO_RESIDUAL_DF", path, "Selected design has no residual degrees of freedom; no statistical fit is attempted.")
    return {"matrix_columns": names, "rank": rank, "row_count": len(selected),
            "column_count": matrix.shape[1], "residual_df": residual, "contrast_estimable": estimable}


def _replication(numerator, denominator, pair_by, report, path):
    pairs = {}
    for arm, samples in [("treatment", numerator), ("control", denominator)]:
        units = [s["biological_unit_id"] for s in samples]
        if any(not u for u in units):
            _issue(report, "BIOLOGICAL_UNIT_REQUIRED", path + "/" + arm, "New computation requires declared biological units; lanes and sample names cannot establish replication.")
        if len(set(u for u in units if u)) < 2:
            _issue(report, "INSUFFICIENT_REPLICATION", path + "/" + arm, "At least two declared independent biological units per selected arm are required by this profile.")
        if any(count > 1 for count in Counter(u for u in units if u).values()) or any(s["technical_replicate_of"] for s in samples):
            _issue(report, "TECHNICAL_OR_REPEATED_OBSERVATIONS", path + "/" + arm, "Technical/repeated observations in one arm require explicit upstream aggregation or a qualified model; they are not independent sample rows.")
        donors = [s["donor_id"] for s in samples if s["donor_id"]]
        if len(set(donors)) < len(donors):
            _issue(report, "NESTED_DONOR_MODEL_REQUIRED", path + "/" + arm, "Multiple observations from one donor in an arm require a separately qualified nested model or aggregation.")
    for key in ["biological_unit_id", "donor_id"]:
        shared = {s[key] for s in numerator if s[key]} & {s[key] for s in denominator if s[key]}
        if shared and pair_by is None:
            _issue(report, "REPEATED_UNITS_REQUIRE_BLOCKING", path, "Shared donors/biological units across conditions require an explicit paired/blocking model; no auto pairing is inferred.")
    if pair_by is not None:
        if pair_by not in _PAIR_KEYS:
            _issue(report, "PAIRING_FIELD", path + "/paired_by", "Choose donor_id or biological_unit_id explicitly; arbitrary sample ordering is not a pair map.")
            return pairs
        left, right = {}, {}
        for samples, target in [(numerator, left), (denominator, right)]:
            for sample in samples:
                value = sample.get(pair_by)
                if not value or value in target:
                    _issue(report, "PAIR_MAPPING_REQUIRED", path, "Paired design requires one nonempty pair ID per arm and no duplicate observations within a pair.")
                else:
                    target[value] = sample["sample_id"]
        if set(left) != set(right) or len(left) < 2:
            _issue(report, "INCOMPLETE_PAIRS", path, "Selected paired design needs at least two complete one-to-one pairs; unmatched samples are never dropped implicitly.")
        for pair in sorted(left.keys() & right.keys()):
            pairs[pair] = {"treatment_sample_id": left[pair], "control_sample_id": right[pair]}
    return pairs


def validate_contrast_rows(metadata_rows, contrast_rows, *, route="bulk_complete", design="~ condition", paired_by=None,
                           intent="compute", assay_mode="bulk_rna", cell_type=None):
    """Validate only the supplied contrast rows; native science is never run."""
    report = _report(route, design, intent, assay_mode)
    support = method_design_support()
    profiles = support["profiles"]
    if not isinstance(route, str) or route not in profiles:
        _issue(report, "UNKNOWN_ROUTE", "/route", "Select an explicitly documented route; method names alone do not qualify an upstream wrapper.")
    if not isinstance(intent, str) or intent not in {"compute", "import"}:
        _issue(report, "INVALID_INTENT", "/intent", "Select compute or import explicitly.")
    if route == "import_only" and intent != "import":
        _issue(report, "IMPORT_ONLY_ROUTE", "/intent", "The import-only route cannot request new computation.")
    if not isinstance(assay_mode, str) or assay_mode not in {"bulk_rna", "single_cell_rna", "single_nucleus_rna"}:
        _issue(report, "UNKNOWN_ASSAY", "/assay_mode", "Preserve the explicit bulk/scRNA/snRNA source assay.")
    elif intent == "compute" and assay_mode != "bulk_rna":
        _issue(report, "BULK_UPSTREAM_SCOPE", "/assay_mode", "These upstream compute profiles are bulk-only; single-cell/nucleus splicing needs its own qualified downstream adapter.")
    if not isinstance(design, str) or not design.strip():
        _issue(report, "DESIGN_FORMULA", "/design", "Record an explicit nonempty design string.")
    if paired_by is not None and (not isinstance(paired_by, str) or paired_by not in _PAIR_KEYS):
        _issue(report, "PAIRING_FIELD", "/paired_by", "Choose donor_id or biological_unit_id explicitly.")
    if cell_type is not None and (not isinstance(cell_type, str) or not cell_type.strip()):
        _issue(report, "CELL_TYPE_SELECTION", "/cell_type", "Cell-type selection must be an explicit nonempty label.")
    meta_ok = _rows_shape(metadata_rows, "metadata", report)
    contrast_ok = _rows_shape(contrast_rows, "contrasts", report)
    if not meta_ok or not contrast_ok or report["issues"]:
        return _finish(report)
    samples_report = validate_sample_rows([{"sample": row.get("sample", "")} for row in metadata_rows], metadata_rows,
                                          input_format="sample_ids", assay_mode=assay_mode)
    for issue in samples_report["issues"]:
        report["issues"].append(dict(issue))
    if not samples_report["valid"]:
        return _finish(report)
    samples = samples_report["samples"]
    seen_ids, seen_pairs = set(), set()
    profile = profiles[route]
    report["profile_status"] = profile["status"]
    for i, row in enumerate(contrast_rows):
        path = f"/contrasts/{i}"
        if not _REQUIRED.issubset(row):
            _issue(report, "MISSING_CONTRAST_COLUMN", path, "Use the native contrast,treatment,control columns explicitly.")
            continue
        cid, treatment, control = (row[k] for k in ["contrast", "treatment", "control"])
        if not _SAFE.fullmatch(cid) or cid in seen_ids:
            _issue(report, "DUPLICATE_OR_INVALID_CONTRAST_ID", path + "/contrast", "Contrast IDs must be unique safe components; no automatic renaming is performed.")
        seen_ids.add(cid)
        if not treatment.strip() or not control.strip() or treatment == control or treatment != treatment.strip() or control != control.strip():
            _issue(report, "CONTRAST_ARMS", path, "Treatment and control must be different, nonempty explicitly declared condition labels.")
        pair = (treatment, control, row.get("cell_type") or cell_type)
        if pair in seen_pairs:
            _issue(report, "DUPLICATE_DIRECTED_COMPARISON", path, "The same directed comparison was selected twice; it is not a second statistical test.")
        if (control, treatment, pair[2]) in seen_pairs:
            _issue(report, "RECIPROCAL_COMPARISON", path, "Both directions were explicitly selected; native effect/statistics remain authoritative and are not negated or recomputed here.", "warning")
        seen_pairs.add(pair)
        unknown = set(row) - _REQUIRED - _OPTIONAL
        if any(row[k] for k in unknown):
            _issue(report, "UNSUPPORTED_CONTRAST_FIELD", path, "Additional design/weight/blocking fields need an explicit adapter; they are not silently ignored.")
        formula = row.get("design") or row.get("formula") or design
        if not formula.strip():
            _issue(report, "DESIGN_FORMULA", path + "/design", "A declared design cannot be blank; imported native models remain explicit unverified provenance.")
        if row.get("design") and row.get("formula") and row["design"] != row["formula"]:
            _issue(report, "CONFLICTING_DESIGN_FIELDS", path, "design and formula disagree; neither is silently preferred.")
        chosen_pair = row.get("paired_by") or paired_by
        if chosen_pair is not None and chosen_pair not in _PAIR_KEYS:
            _issue(report, "PAIRING_FIELD", path + "/paired_by", "Row pairing field must explicitly name donor_id or biological_unit_id.")
        paired_text = row.get("paired", "").lower()
        bool_values = {"true": True, "1": True, "yes": True, "false": False, "0": False, "no": False}
        declared_paired = bool_values.get(paired_text)
        if paired_text and paired_text not in bool_values:
            _issue(report, "PAIRED_FLAG_SYNTAX", path + "/paired", "paired must use true/false, yes/no or 1/0 explicitly, including for imported provenance.")
        chosen_cell = row.get("cell_type") or cell_type
        scope = [s for s in samples if s["group"] in {treatment, control} and (chosen_cell is None or s["cell_type"] == chosen_cell)]
        numerator = [s for s in scope if s["group"] == treatment]
        denominator = [s for s in scope if s["group"] == control]
        if not numerator or not denominator:
            _issue(report, "UNKNOWN_OR_EMPTY_CONDITION", path, "Both selected conditions need samples within the explicitly selected cell-type scope.")
        if len({s["cell_type"] for s in scope}) > 1:
            _issue(report, "MIXED_CELL_TYPES", path, "Several cell types occur in the selected arms; select a cell type explicitly rather than pooling them.")
        record = {"contrast_id": cid, "treatment": treatment, "control": control,
                  "effect_direction": "treatment_minus_control", "cell_type": chosen_cell,
                  "effect_direction_scope": "declared_comparison_only; native columns remain authoritative",
                  "source_contrast_fields": dict(row),
                  "numerator_sample_ids": [s["sample_id"] for s in numerator],
                  "denominator_sample_ids": [s["sample_id"] for s in denominator],
                  "selected_sample_ids": [s["sample_id"] for s in scope],
                  "fit_sample_ids": None if intent == "import" else [s["sample_id"] for s in scope],
                  "fit_scope": "native_fit_unverified" if intent == "import" else "planned_per_selected_contrast_source_wrapper",
                  "design": {"requested_formula": formula, "terms": None, "paired_by": chosen_pair,
                             "declared_paired": declared_paired, "pair_map": {}, "matrix_columns": [], "rank": None, "row_count": len(scope),
                             "column_count": None, "residual_df": None, "contrast_estimable": None},
                  "methods": [{"method_id": name, "native_effect": support["methods"][name]["native_effect"],
                               "profile_status": "imported_native_design_unverified" if intent == "import" else profile["status"],
                               "direction_transform_applied": False,
                               "native_semantics_scope": "pinned_wrapper_convention; imported_source_unverified" if intent == "import" else "pinned_wrapper_convention; no_execution_performed",
                               "native_input_assignment_scope": "expected_pinned_wrapper; imported_assignment_unverified" if intent == "import" else "expected_pinned_wrapper; no_execution_performed",
                               "native_input_assignment": {"deseq2": ["condition", treatment, control],
                                    "rmats": {"group1": treatment, "group2": control},
                                    "suppa": {"input1": control, "input2": treatment}}.get(name, "native_headers_authoritative"),
                               "runtime_ready": False} for name in profile["methods"]]}
        report["selected_contrasts"].append(record)
        if intent == "import":
            _issue(report, "NATIVE_IMPORT_UNVERIFIED", path, "Imported native design is recorded without validating the original fit, changing P/FDR or authorising new computation.", "warning")
            continue
        for label in [treatment, control]:
            if not _CONDITION.fullmatch(label):
                _issue(report, "NATIVE_CONDITION_ID", path, "This pinned upstream requires condition IDs starting with a letter and containing letters, digits or underscores; display labels need an explicit mapping.")
        terms = _terms(formula, report, path + "/design")
        if terms is None:
            continue
        record["design"]["terms"] = terms
        if set(terms) - set(profile["allowed_terms"]):
            _issue(report, "BLOCKED_CURRENT_ROUTE", path + "/design", "The selected route does not implement this design; underlying method capabilities do not extend its wrapper automatically.")
        if chosen_pair and not profile["paired"]:
            _issue(report, "PAIRED_ROUTE_UNSUPPORTED", path, "This route does not implement paired testing; select the explicitly qualified adapter instead.")
        if ("subject" in terms) != bool(chosen_pair):
            _issue(report, "PAIR_DESIGN_MISMATCH", path, "A subject term requires an explicit pair mapping, and a selected paired mapping requires the subject term.")
        if declared_paired is not None and declared_paired != bool(chosen_pair):
            _issue(report, "PAIRED_FLAG_CONFLICT", path, "paired must be true/false and agree with the explicit pairing/model selection.")
        record["design"]["pair_map"] = _replication(numerator, denominator, chosen_pair, report, path)
        if scope:
            record["design"].update(_design_matrix(scope, terms, treatment, chosen_pair, report, path + "/design"))
    return _finish(report)


def validate_contrastsheets(metadata_path, contrast_path, **kwargs):
    """Read CSV metadata and only the selected contrasts; never modify either."""
    report = _report(kwargs.get("route", "bulk_complete"), kwargs.get("design", "~ condition"),
                     kwargs.get("intent", "compute"), kwargs.get("assay_mode", "bulk_rna"))
    metadata = _read_csv(metadata_path, _META_REQUIRED, report, "metadata")
    contrasts = _read_csv(contrast_path, _REQUIRED, report, "contrasts")
    if report["issues"]:
        return _finish(report)
    return validate_contrast_rows(metadata, contrasts, **kwargs)
