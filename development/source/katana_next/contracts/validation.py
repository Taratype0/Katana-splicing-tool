"""Read-only structural, identity, provenance and optional local-file checks.

The exchange contract version is independent of the historical platform schema.
Passing validation is not design estimability, method eligibility, or run approval.
"""
from __future__ import annotations

from collections import Counter
from importlib.resources import files
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any
from urllib.parse import urlsplit
import hashlib
import json
import math
import os

from jsonschema import Draft202012Validator

CONTRACT_VERSION = "1.0.0"
KINDS = {"project", "asset", "run", "result"}


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON object key: {key}")
        result[key] = value
    return result


def load_document(path: str | Path) -> dict:
    """Load strict JSON without creating a project/catalog/output directory."""
    def reject_constant(value):
        raise ValueError(f"Non-JSON numeric constant: {value}")
    return json.loads(Path(path).read_text(encoding="utf-8-sig"),
                      object_pairs_hook=_pairs, parse_constant=reject_constant)


def _schema(kind):
    resource = files("katana_next.contracts").joinpath("schemas/contract-v1.schema.json")
    schema = json.loads(resource.read_text(encoding="utf-8"))
    schema["$ref"] = f"#/$defs/{kind}"
    Draft202012Validator.check_schema(schema)
    return schema


class _Check:
    def __init__(self, base_dir, check_files, verify_hashes):
        self.base = Path(base_dir).resolve() if base_dir is not None else None
        self.check_files = check_files or verify_hashes
        self.verify_hashes = verify_hashes
        self.issues = []
        self.file_checks = {"local_files_checked": 0, "hashes_checked": 0, "remote_assets_not_fetched": 0}

    def issue(self, code, path, message, severity="error"):
        self.issues.append({"code": code, "path": path, "message": message, "severity": severity})

    def refs(self, values, known, path, kind):
        for i, value in enumerate(values):
            if value not in known:
                self.issue("UNKNOWN_" + kind.upper(), f"{path}/{i}", f"Unknown {kind} ID: {value}")

    def index(self, values, key, path):
        counts = Counter(x[key] for x in values)
        for value, count in counts.items():
            if count > 1:
                self.issue("DUPLICATE_ID", path, f"{key} is repeated: {value}")
        return {x[key]: x for x in values}

    def finite(self, value, path=""):
        if isinstance(value, float) and not math.isfinite(value):
            self.issue("NONFINITE_NUMBER", path or "/", "JSON numbers must be finite; use null for unavailable values.")
        elif isinstance(value, dict):
            for key, item in value.items():
                self.finite(item, path + "/" + str(key))
        elif isinstance(value, list):
            for i, item in enumerate(value):
                self.finite(item, path + "/" + str(i))

    def asset(self, asset, path):
        loc = asset["location"]
        text = loc["path"]
        kind = loc["kind"]
        resolved = None
        if kind == "project_relative":
            value = PurePosixPath(text)
            if "\\" in text or value.is_absolute() or PureWindowsPath(text).drive or ".." in value.parts or text.startswith("~") or ":" in text or text in ["", "."]:
                self.issue("INVALID_RELATIVE_PATH", path + "/location/path", "Use a portable relative file path with forward slashes and no parent traversal.")
                return
            if self.base is not None:
                resolved = (self.base / text).resolve()
                if not resolved.is_relative_to(self.base):
                    self.issue("PATH_ESCAPE", path + "/location/path", "Resolved relative path escapes the project base directory.")
                    return
        elif kind == "local_absolute":
            windows = PureWindowsPath(text).is_absolute()
            posix = PurePosixPath(text).is_absolute()
            if not (windows or posix) or ".." in PureWindowsPath(text).parts or ".." in PurePosixPath(text).parts:
                self.issue("INVALID_ABSOLUTE_PATH", path + "/location/path", "Declare a full Windows/UNC or POSIX absolute path without traversal.")
                return
            foreign = (windows and os.name != "nt") or (posix and not windows and os.name == "nt")
            if foreign and self.check_files and asset["availability"] == "present":
                self.issue("FOREIGN_PATH_UNVERIFIABLE", path + "/location/path", "This host cannot verify the declared absolute path flavour; record it remotely or validate on its host.")
                return
            if not foreign:
                resolved = Path(text)
        else:
            uri = urlsplit(text)
            if uri.scheme not in {"https", "s3", "gs", "ssh"} or not uri.netloc or uri.username or uri.password or not uri.path:
                self.issue("INVALID_REMOTE_PATH", path + "/location/path", "Use an explicit https/s3/gs/ssh resource URI without embedded credentials.")
                return
            if asset["availability"] == "present":
                self.issue("REMOTE_AVAILABILITY", path + "/availability", "Use remote/not_checked/missing for a remote locator, not a locally present claim.")
            if self.check_files:
                self.file_checks["remote_assets_not_fetched"] += 1
                self.issue("REMOTE_NOT_FETCHED", path + "/location", "Remote content was not fetched or verified.", "warning")
            return
        if asset["availability"] == "remote":
            self.issue("LOCAL_AVAILABILITY", path + "/availability", "A local locator cannot use remote availability.")
        if not self.check_files or asset["availability"] != "present":
            return
        if resolved is None:
            self.issue("BASE_DIR_REQUIRED", path + "/location", "A project base directory is needed to check a relative file.")
            return
        if not resolved.is_file():
            self.issue("FILE_MISSING", path + "/location", "A file declared present is not accessible at the resolved location.")
            return
        self.file_checks["local_files_checked"] += 1
        stat = resolved.stat()
        if asset["bytes"] is not None and stat.st_size != asset["bytes"]:
            self.issue("SIZE_MISMATCH", path + "/bytes", "Actual file size differs from the declared size.")
        if self.verify_hashes:
            if asset["sha256"] is None:
                self.issue("HASH_REQUIRED", path + "/sha256", "A SHA256 is required for requested content verification.")
            else:
                h = hashlib.sha256()
                with resolved.open("rb") as f:
                    for block in iter(lambda: f.read(1024 * 1024), b""):
                        h.update(block)
                self.file_checks["hashes_checked"] += 1
                if h.hexdigest() != asset["sha256"]:
                    self.issue("HASH_MISMATCH", path + "/sha256", "Actual file content differs from its declared SHA256.")
        end = resolved.stat()
        if (stat.st_size, stat.st_mtime_ns) != (end.st_size, end.st_mtime_ns):
            self.issue("FILE_CHANGED_DURING_CHECK", path + "/location", "The file changed while being verified; no stable verification is claimed.")

    def result_state(self, result, path):
        complete = result["status"] == "completed"
        if not complete:
            if result["reason"] is None:
                self.issue("MISSING_REASON", path + "/reason", "Unavailable/unverified results require an explicit reason.")
            if any(result[x] is not None for x in ["row_count", "tested_count", "significant_count"]) or result["findings"] not in {"not_assessed", "not_applicable"}:
                self.issue("MISSING_IS_NOT_ZERO", path, "Missing/failed/omitted/not-estimable results use null counts and no finding claim.")
            return
        if result["reason"] is not None or result["row_count"] is None or not result["output_asset_ids"]:
            self.issue("INCOMPLETE_RESULT", path, "Completed results require output assets, a known row count and no missing-result reason.")
        tested = result["tested_count"]
        significant = result["significant_count"]
        if significant is not None and (tested is None or significant > tested):
            self.issue("INVALID_FINDING_COUNTS", path, "A significant count requires a tested count and cannot exceed it.")
        findings = result["findings"]
        if findings in {"zero_findings", "has_findings"}:
            if tested is None or tested <= 0 or significant is None or not result["family_id"] or not result["tested_universe_asset_id"]:
                self.issue("FINDINGS_REQUIRE_TESTED_UNIVERSE", path, "Finding claims require a nonempty tested universe, count and original family identity.")
            if findings == "zero_findings" and significant != 0:
                self.issue("ZERO_FINDINGS_COUNT", path, "Zero findings requires significant_count=0 after a completed assessed family.")
            if findings == "has_findings" and (significant is None or significant <= 0):
                self.issue("POSITIVE_FINDINGS_COUNT", path, "Positive findings requires significant_count>0.")
        elif significant is not None or tested is not None:
            self.issue("UNASSESSED_COUNTS", path, "Unassessed/not-applicable finding states do not carry tested or significant counts.")

    def run_state(self, run, path):
        if run["status"] in {"blocked", "failed", "cancelled", "omitted"} and run["reason"] is None:
            self.issue("MISSING_REASON", path + "/reason", "Blocked/failed/cancelled/omitted runs require a reason.")
        if run["execution_kind"] == "new_computation" and run["status"] == "completed":
            evidence = run["execution_evidence"]
            if evidence["status"] != "verified" or not evidence["receipt_asset_ids"]:
                self.issue("COMPLETION_EVIDENCE_REQUIRED", path, "Completed new computation requires traceable completion receipt assets; scheduler submission alone is insufficient.")
            software = run["software"]
            if not software["source_commit"] and not software["container_digest"]:
                self.issue("SOURCE_IDENTITY_REQUIRED", path + "/software", "Record the executing source commit or pinned container identity.")
        if run["execution_kind"] == "import" and run["status"] == "completed":
            provenance = run["provenance"]
            if not all(provenance[x] for x in ["source_schema", "source_version", "source_artifact_id"]):
                self.issue("IMPORT_PROVENANCE_REQUIRED", path + "/provenance", "A completed import preserves source schema/version/artifact identity; no new native job is required.")

    def project(self, doc):
        groups = [("samples", "sample_id"), ("contrasts", "contrast_id"), ("references", "reference_id"),
                  ("assets", "asset_id"), ("methods", "method_id"), ("runs", "run_id"), ("results", "result_id")]
        indexes = {name: self.index(doc[name], key, "/" + name) for name, key in groups}
        samples, contrasts, references, assets, methods, runs = [indexes[x] for x in ["samples", "contrasts", "references", "assets", "methods", "runs"]]
        single = doc["assay_mode"] != "bulk_rna"
        if (single and doc["analysis_level"] == "bulk") or (not single and doc["analysis_level"] != "bulk"):
            self.issue("ASSAY_LEVEL_MISMATCH", "/analysis_level", "Preserve bulk versus cell/pseudobulk origin; pseudobulk is not relabelled raw bulk.")
        if single and doc["input_kind"] == "bulk_fastq":
            self.issue("RAW_SINGLE_CELL_UNSUPPORTED", "/input_kind", "Raw scRNA/snRNA FASTQ upstream is outside this bulk pipeline; use qualified processed inputs.")
        for i, sample in enumerate(doc["samples"]):
            path = f"/samples/{i}"
            if sample["source_assay"] != doc["assay_mode"]:
                self.issue("SOURCE_ASSAY_MISMATCH", path + "/source_assay", "Keep the original assay identity; do not relabel single-cell origin as bulk.")
            if doc["analysis_level"] == "pseudobulk":
                aggregation = sample.get("aggregation")
                if not sample["donor_id"] or not sample["cell_type"] or not aggregation:
                    self.issue("PSEUDOBULK_LINEAGE_REQUIRED", path, "Pseudobulk requires donor, cell type and declared aggregation lineage.")
                elif aggregation["source_assay"] != sample["source_assay"] or not {"donor_id", "cell_type"}.issubset(aggregation["group_by"]):
                    self.issue("PSEUDOBULK_GROUPING", path + "/aggregation", "Keep source assay and donor-by-cell-type aggregation keys.")
            if sample.get("aggregation"):
                self.refs(sample["aggregation"]["source_asset_ids"], assets, path + "/aggregation/source_asset_ids", "asset")
        for i, contrast in enumerate(doc["contrasts"]):
            path = f"/contrasts/{i}"
            self.refs(contrast["sample_ids"], samples, path + "/sample_ids", "sample")
            if contrast["kind"] == "group_difference":
                num, den = contrast["numerator_sample_ids"], contrast["denominator_sample_ids"]
                self.refs(num + den, samples, path, "sample")
                if set(num) & set(den):
                    self.issue("CONTRAST_ARM_OVERLAP", path, "Numerator and denominator samples must be disjoint.")
                if set(num) | set(den) != set(contrast["sample_ids"]):
                    self.issue("CONTRAST_SCOPE_MISMATCH", path, "Declared comparison scope must contain exactly both arms.")
                if contrast["effect_direction"] != "numerator_minus_denominator":
                    self.issue("CONTRAST_DIRECTION", path, "Group differences explicitly use numerator minus denominator; no alphabetic fallback.")
                if contrast["paired"]:
                    pairs = contrast.get("pair_ids", {})
                    if set(pairs) != set(contrast["sample_ids"]) or Counter(pairs.get(x) for x in num) != Counter(pairs.get(x) for x in den) or any(n != 1 for n in Counter(pairs.get(x) for x in num).values()):
                        self.issue("PAIR_MAPPING_REQUIRED", path, "Paired comparisons require an explicit one-to-one pair map covering both arms.")
            elif contrast["effect_direction"] != "declared_linear_combination" or not any(contrast["weights"].values()):
                self.issue("LINEAR_CONTRAST_DEFINITION", path, "Record a nonzero coefficient combination and its declared effect direction.")
        for i, reference in enumerate(doc["references"]):
            self.refs(reference["asset_ids"], assets, f"/references/{i}/asset_ids", "asset")
            if reference["identity_status"] == "known" and (reference["taxon_id"] is None or not reference["genome_build"]):
                self.issue("REFERENCE_IDENTITY", f"/references/{i}", "Known reference identity requires taxon and genome build; unresolved identity remains explicit.")
        for i, asset in enumerate(doc["assets"]):
            path = f"/assets/{i}"
            self.refs(asset["sample_ids"], samples, path + "/sample_ids", "sample")
            self.refs([asset["reference_id"]] if asset["reference_id"] else [], references, path + "/reference_id", "reference")
            self.asset(asset, path)
        for i, run in enumerate(doc["runs"]):
            path = f"/runs/{i}"
            self.run_state(run, path)
            for field, known, name in [("sample_ids", samples, "sample"), ("contrast_ids", contrasts, "contrast"), ("input_asset_ids", assets, "asset"), ("reference_ids", references, "reference"), ("parent_run_ids", runs, "run")]:
                self.refs(run[field], known, path + "/" + field, name)
            self.refs(run["execution_evidence"]["receipt_asset_ids"], assets, path + "/execution_evidence/receipt_asset_ids", "asset")
            self.refs([run["method_id"]], methods, path + "/method_id", "method")
            if run["execution_kind"] == "import" and run["provenance"]["source_artifact_id"]:
                self.refs([run["provenance"]["source_artifact_id"]], assets, path + "/provenance/source_artifact_id", "asset")
            method = methods.get(run["method_id"])
            compatible = {"import": {"import"}, "render": {"render"}, "new_computation": {"compute", "annotate", "diagnostic"}}
            if method and method["operation"] not in compatible[run["execution_kind"]]:
                self.issue("METHOD_EXECUTION_KIND", path + "/execution_kind", "Method operation describes this adapter/worker operation and must agree with import/new-computation/render intent.")
            if run["execution_kind"] == "new_computation" and run["status"] == "completed":
                for aid in run["execution_evidence"]["receipt_asset_ids"]:
                    if aid in assets and (assets[aid]["role"] != "evidence" or assets[aid]["availability"] not in {"present", "remote"} or assets[aid]["sha256"] is None or assets[aid]["bytes"] is None):
                        self.issue("COMPLETION_RECEIPT_IDENTITY", path + "/execution_evidence", "A verified completion receipt requires an available evidence asset with size and hash, not a missing or unidentifiable placeholder.")
            if run["run_id"] in run["parent_run_ids"]:
                self.issue("RUN_CYCLE", path + "/parent_run_ids", "A run cannot depend on itself.")
            for cid in run["contrast_ids"]:
                if cid in contrasts and not set(contrasts[cid]["sample_ids"]).issubset(run["sample_ids"]):
                    self.issue("RUN_CONTRAST_SAMPLES", path, "The model/run sample scope must include its comparison samples.")
            if doc["intent"] in {"import", "report"} and run["execution_kind"] == "new_computation":
                self.issue("INTENT_COMPUTE_CONFLICT", path + "/execution_kind", "Import/report intent does not authorize new scientific computation.")
            active = run["execution_kind"] == "new_computation" and run["status"] in {"running", "completed"}
            if active:
                method = methods.get(run["method_id"])
                if method and method["source_sha256"] is None:
                    self.issue("METHOD_SOURCE_REQUIRED", path + "/method_id", "Active/new completed computation requires the actual method source fingerprint.")
                if single and method and (method["inference_unit"] == "biological_sample" or (method["inference_unit"] == "donor" and any(not samples[s]["donor_id"] for s in run["sample_ids"] if s in samples))):
                    self.issue("DONOR_INFERENCE_REQUIRED", path, "Preserve explicit donor mapping for declared donor inference; scRNA/snRNA rows cannot be relabelled bulk biological samples. Other declared native/test units retain their identities and require method-specific qualification later.")
                if any(references[r]["identity_status"] != "known" for r in run["reference_ids"] if r in references):
                    self.issue("UNRESOLVED_COMPUTE_REFERENCE", path, "Resolve each requested physical reference before active computation.")
                if any(assets[a]["availability"] not in {"present", "remote"} for a in run["input_asset_ids"] if a in assets):
                    self.issue("UNAVAILABLE_COMPUTE_INPUT", path, "Active computation cannot claim missing/not-checked inputs as ready.")
        self._run_cycles(runs)
        for i, result in enumerate(doc["results"]):
            path = f"/results/{i}"
            self.result_state(result, path)
            self.refs([result["run_id"]], runs, path + "/run_id", "run")
            self.refs(result["output_asset_ids"], assets, path + "/output_asset_ids", "asset")
            self.refs([result["tested_universe_asset_id"]] if result["tested_universe_asset_id"] else [], assets, path + "/tested_universe_asset_id", "asset")
            self.refs([result["contrast_id"]] if result["contrast_id"] else [], contrasts, path + "/contrast_id", "contrast")
            run = runs.get(result["run_id"])
            if run is None:
                continue
            if result["contrast_id"] and result["contrast_id"] not in run["contrast_ids"]:
                self.issue("RESULT_CONTRAST_SCOPE", path, "Result comparison must belong to its run.")
            method = methods.get(run["method_id"])
            if method:
                for field in ["measurement_unit", "inference_unit", "statistics_origin"]:
                    if result[field] != method[field]:
                        self.issue("RESULT_METHOD_UNIT_OR_ORIGIN", path + "/" + field, "Result units/statistical origin must match the declared method; a cell P value cannot be relabelled donor inference.")
            if result["status"] == "completed":
                if run["status"] != "completed":
                    self.issue("RESULT_RUN_NOT_COMPLETE", path + "/run_id", "Completed results require their producing/import run to be completed.")
                if run["execution_kind"] == "import" and run["provenance"]["source_completion"] != "completed":
                    self.issue("SOURCE_COMPLETION_UNKNOWN", path, "An import can succeed while source science is unknown/failed; mark its result unverified instead of scientific completion.")
                universe = assets.get(result["tested_universe_asset_id"])
                if result["findings"] in {"zero_findings", "has_findings"} and universe and (universe["availability"] not in {"present", "remote"} or universe["sha256"] is None or universe["bytes"] is None):
                    self.issue("TESTED_UNIVERSE_IDENTITY", path + "/tested_universe_asset_id", "Assessed findings require the retained tested-universe asset identity; a missing placeholder cannot establish a tested family.")
                for aid in result["output_asset_ids"]:
                    if aid in assets and assets[aid]["measurement_unit"] is not None and assets[aid]["measurement_unit"] != result["measurement_unit"]:
                        self.issue("ASSET_RESULT_UNIT", path + "/measurement_unit", "Output asset and result measurement units disagree; preserve the declared native scale.")
                    if aid in assets and (assets[aid]["availability"] not in {"present", "remote"} or assets[aid]["sha256"] is None or assets[aid]["bytes"] is None):
                        self.issue("COMPLETED_ASSET_IDENTITY", path + "/output_asset_ids", "Completed output assets require declared availability, size and content hash; actual verification is a separate read-only option.")

    def _run_cycles(self, runs):
        from collections import deque, defaultdict
        indegree = {key: 0 for key in runs}
        children = defaultdict(list)
        for key, run in runs.items():
            for parent in run["parent_run_ids"]:
                if parent in runs:
                    indegree[key] += 1
                    children[parent].append(key)
        queue = deque(key for key, degree in indegree.items() if degree == 0)
        seen = 0
        while queue:
            key = queue.popleft()
            seen += 1
            for child in children[key]:
                indegree[child] -= 1
                if indegree[child] == 0:
                    queue.append(child)
        if seen != len(runs):
            self.issue("RUN_CYCLE", "/runs", "Run lineage must be acyclic; shared models are referenced without recursive dependencies.")


def validate_document(document: Any, *, base_dir=None, check_files=False, verify_hashes=False) -> dict:
    """Validate metadata and optional local bytes; never mutate the document or run tools."""
    check = _Check(base_dir, check_files, verify_hashes)
    kind = document.get("kind") if isinstance(document, dict) else None
    if not isinstance(kind, str) or kind not in KINDS or document.get("schema_version") != CONTRACT_VERSION:
        check.issue("UNSUPPORTED_CONTRACT", "/", "Expected an explicit Katana exchange contract 1.0.0 and project/asset/run/result kind. Historical platform schema 2.0 requires an explicit adapter; no automatic migration is performed.")
    else:
        check.finite(document)
        validator = Draft202012Validator(_schema(kind))
        for error in sorted(validator.iter_errors(document), key=lambda e: str(list(e.absolute_path))):
            check.issue("SCHEMA", "/" + "/".join(str(p) for p in error.absolute_path), error.message)
        if not check.issues:
            if kind == "project":
                check.project(document)
            elif kind == "asset":
                check.asset(document, "")
            elif kind == "run":
                check.run_state(document, "")
            else:
                check.result_state(document, "")
            if kind != "project":
                check.issue("PROJECT_CONTEXT_NOT_CHECKED", "/", "Standalone structure/state checked; cross-record IDs and provenance require the containing project.", "warning")
    return {"contract_version": CONTRACT_VERSION, "document_kind": kind,
            "valid": not any(i["severity"] == "error" for i in check.issues),
            "scope": "metadata_and_local_hashes" if verify_hashes else "metadata_and_local_files" if check_files else "metadata_only",
            "execution_authorized": False, "scientific_tools_executed": 0,
            "recorded_states": {name: dict(Counter(r[field] for r in (document.get(name) if isinstance(document.get(name), list) else []) if isinstance(r, dict) and isinstance(r.get(field), str)))
                                for name, field in [("runs", "status"), ("results", "status"), ("assets", "availability")]} if kind == "project" else {},
            "issues": check.issues, "file_verification": check.file_checks}


def validate_file(path: str | Path, *, check_files=False, verify_hashes=False) -> dict:
    path = Path(path)
    try:
        document = load_document(path)
        return validate_document(document, base_dir=path.parent, check_files=check_files, verify_hashes=verify_hashes)
    except (ValueError, OSError) as error:
        return {"contract_version": CONTRACT_VERSION, "document_kind": None, "valid": False,
                "scope": "document_load", "execution_authorized": False, "scientific_tools_executed": 0,
                "issues": [{"code": "DOCUMENT_LOAD", "path": "/", "severity": "error", "message": str(error)}],
                "file_verification": {"local_files_checked": 0, "hashes_checked": 0, "remote_assets_not_fetched": 0}}
