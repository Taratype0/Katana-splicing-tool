"""Explicit presentation-only recovery of a proven MotifTrace render failure.

This never changes the failed task or promotes its status to completed. A new
task executes validation, file transport, presentation and artifact registration.
The original comparison's statistical tables are preserved byte for byte.
"""
from pathlib import Path

from .runtime import canonical_hash, file_hash, file_record, read_json, verify_files


RENDER_FAILURE = "Portable event report exceeds 64MiB; complete tables retained, no truncation"


def prepare_presentation_recovery(path, *, plan, task, base_sha, manifest_sha):
    from katana_next.reports.motif_native_index import SCIENTIFIC_FILES
    failure_path = Path(path).expanduser().resolve()
    failure = read_json(failure_path)
    expected = {
        "schema": "katana.flow.task_receipt/v1", "status": "failed",
        "error": RENDER_FAILURE, "plan_id": plan["plan_id"],
        "task_id": task["id"], "task_contract_sha256": canonical_hash(task),
        "workflow_id": "prioritize.motiftrace_consequences",
        "source_project_manifest_sha256": manifest_sha,
        "base_request_sha256": base_sha,
    }
    if (task.get("workflow_id") != expected["workflow_id"] or
        any(failure.get(key) != value for key, value in expected.items()) or
        not failure.get("finished_at") or
        failure.get("receipt_sha256") != canonical_hash({k: v for k, v in failure.items() if k != "receipt_sha256"})):
        raise ValueError("Presentation recovery requires the exact signed render failure of this unchanged task and its dependencies")
    roots = list((failure_path.parent / "results" / "artifacts" /
                  "prioritize.motiftrace_consequences").glob("*/event_consequences/RESULTS.json"))
    if len(roots) != 1:
        raise ValueError("Presentation recovery requires one retained scientific result family")
    source = roots[0].parent.resolve()
    if (source / "COMPLETE.json").exists():
        raise ValueError("An already completed family must use normal result reuse, not failed-presentation recovery")
    records, members = [file_record(failure_path)], {}
    for name in sorted(SCIENTIFIC_FILES):
        member = source / name
        if member.is_symlink() or not member.resolve().is_relative_to(source):
            raise ValueError("Retained scientific family must contain owned regular files")
        record = file_record(member)
        records.append(record)
        members[name] = {"bytes": record["size_bytes"], "sha256": record["sha256"]}
    evidence = {
        "failure_receipt": str(failure_path), "failure_receipt_sha256": file_hash(failure_path),
        "task_id": task["id"], "task_contract_sha256": canonical_hash(task),
        "scientific_plan_sha256": members["ANALYSIS_PLAN_BEFORE_SCORING.json"]["sha256"],
        "members": members,
    }
    return PresentationRecovery(source, evidence, records)


class PresentationRecovery:
    def __init__(self, source, evidence, inputs):
        self.source, self.evidence, self.inputs = source, evidence, inputs
        self.calls = 0
        self.input_validation = None

    def run(self, complete, selected_ids, models, fasta, destination, **kwargs):
        from katana_next.reports.motif_native_index import (
            validate_recorded_consequences_inputs, finalize_recorded_consequences,
        )
        if self.calls:
            raise ValueError("One presentation recovery may finalize only one complete comparison")
        verify_files(self.inputs)
        checked = {key: value for key, value in kwargs.items() if key != "check_cancelled"}
        self.input_validation = validate_recorded_consequences_inputs(
            self.source, complete, selected_ids, models, fasta,
            failure_receipt=self.evidence["failure_receipt"], **checked,
        )
        kwargs.get("check_cancelled", lambda: None)()
        result = finalize_recorded_consequences(self.source, destination, self.evidence)
        verify_files(self.inputs)
        self.calls += 1
        return result

    def execution_evidence(self):
        if self.calls != 1:
            raise ValueError("The new task did not execute its declared presentation recovery")
        return {
            "presentation_recovery": self.evidence,
            "scientific_input_validation": self.input_validation,
            "scientific_computation_performed": False,
            "preserved_scientific_family": str(self.source),
            "original_failed_attempt_preserved": True,
        }
