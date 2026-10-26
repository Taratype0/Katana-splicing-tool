"""Reuse immutable observed networks and distribute missing bootstrap resamples.

The observed full-data partition is never fitted by this adapter. Each worker
computes exactly one predeclared biological-unit resample; reduction replays its
verified labels through the unchanged consensus function, including kME and the
continuous membership evidence. Scientific modules are never edited.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
from threading import RLock

import numpy as np
import pandas as pd

_REPLAY_LOCK = RLock()


def digest(path: str | Path) -> str:
    value = sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def _json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as target:
        json.dump(value, target, indent=2, sort_keys=True)
        target.write("\n")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def kernel_identity() -> dict:
    import katana_cosplice
    root = Path(katana_cosplice.__file__).resolve().parent
    return {p.relative_to(root).as_posix(): digest(p)
            for p in sorted(root.rglob("*.py"))}


def _config_checks(config) -> None:
    if config.network.engine != "tom_hierarchical" or config.network.module_detection != "wgcna_blockwise":
        raise ValueError("Adoption currently supports the exact WGCNA TOM backbone only")
    if config.consensus.clustering_method != "backbone_confidence_only":
        raise ValueError("Adoption must preserve the full observed backbone")
    if config.uncertainty.enabled and config.uncertainty.network_mode != "diagnostic_only":
        raise ValueError("Posterior/reliability bootstrap requires an explicitly supported input contract")
    if config.gene_prior.enabled or config.multiobjective_selection.enabled:
        raise ValueError("Prior and tuned bootstrap adoption require an additional explicit contract")
    if config.consensus.apply_threshold_labels:
        raise ValueError("The approved default retains continuous metrics without threshold labels")
    if not 1 <= config.consensus.bootstrap_iterations <= np.iinfo(np.uint16).max:
        raise ValueError("Unsupported bootstrap iteration count")


def prepare_plan(*, config: str | Path, matrix: str | Path,
                 assignments: str | Path, samples: str | Path,
                 destination: str | Path, source_acceptance: str | Path | None = None,
                 expected_source_hashes: dict | None = None) -> Path:
    """Freeze paths, full input hashes, biological draws and frozen kernel identity."""
    from katana_cosplice.config import load_config
    from katana_cosplice.resampling import bootstrap_plan

    target = Path(destination).resolve()
    if target.exists():
        raise FileExistsError("Bootstrap plan destination already exists")
    cfg = load_config(config)
    _config_checks(cfg)
    inputs = {name: Path(path).resolve() for name, path in
              {"config": config, "matrix": matrix, "assignments": assignments,
               "samples": samples}.items()}
    identities = {name: {"path": str(p), "bytes": p.stat().st_size, "sha256": digest(p)}
                  for name, p in inputs.items()}
    if expected_source_hashes:
        for name, expected in expected_source_hashes.items():
            if name not in identities or identities[name]["sha256"] != expected:
                raise ValueError("Observed source identity mismatch: " + name)
    psi = pd.read_csv(inputs["matrix"], sep="\t", index_col=0)
    assigned = pd.read_csv(inputs["assignments"], sep="\t", low_memory=False)
    metadata = pd.read_csv(inputs["samples"], sep="\t", low_memory=False)
    if psi.empty or not psi.index.is_unique or not psi.columns.is_unique:
        raise ValueError("Observed matrix identities must be nonempty and unique")
    if not np.isfinite(psi.to_numpy(dtype=float)).all():
        raise ValueError("Observed matrix contains nonfinite values")
    if assigned.event_key.duplicated().any() or set(assigned.event_key) != set(psi.index):
        raise ValueError("Observed assignments do not exactly cover matrix events")
    resampling = bootstrap_plan(metadata, psi.columns)
    resampling.require_resampleable()
    rng = np.random.default_rng(cfg.consensus.random_seed)
    draws = []
    for iteration in range(1, cfg.consensus.bootstrap_iterations + 1):
        indices = resampling.draw(rng)
        draws.append({"iteration": iteration, "indices": indices.tolist(),
                      "draw_sha256": resampling.draw_hash(indices)})
    proof = None
    if source_acceptance is not None:
        p = Path(source_acceptance).resolve()
        accepted = json.loads(p.read_text(encoding="utf-8"))
        if accepted.get("passed") is not True:
            raise ValueError("Observed network source has not passed its existing-results acceptance")
        proven = {row["sha256"] for row in accepted.get("files", [])}
        if not {identities[k]["sha256"] for k in ("matrix", "assignments", "samples")} <= proven:
            raise ValueError("Observed acceptance does not prove the exact bootstrap inputs")
        proof = {"path": str(p), "bytes": p.stat().st_size, "sha256": digest(p)}
    value = {"schema": "katana.observed_network_bootstrap_plan/v1", "created_utc": _now(),
             "inputs": identities, "observed_acceptance": proof,
             "events": len(psi), "samples": len(psi.columns),
             "event_order_sha256": sha256(json.dumps(list(map(str, psi.index))).encode()).hexdigest(),
             "sample_order": list(map(str, psi.columns)), "bootstrap_iterations": len(draws),
             "resampling": resampling.audit(), "draws": draws, "kernel": kernel_identity(),
             "full_data_network_refit": False,
             "numerical_parameters_changed": False,
             "interpretation": "Internal biological-unit resampling confidence on the immutable observed partition"}
    target.mkdir(parents=True, exist_ok=False)
    plan = target / "PLAN.json"
    _json(plan, value)
    return plan


def load_plan(path: str | Path) -> tuple[dict, object, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    from katana_cosplice.config import load_config
    plan = json.loads(Path(path).read_text(encoding="utf-8"))
    if plan.get("schema") != "katana.observed_network_bootstrap_plan/v1":
        raise ValueError("Unsupported bootstrap adoption plan")
    for key, identity in plan["inputs"].items():
        p = Path(identity["path"])
        if not p.is_file() or p.stat().st_size != identity["bytes"] or digest(p) != identity["sha256"]:
            raise ValueError("Immutable bootstrap input changed: " + key)
    proof = plan.get("observed_acceptance")
    if proof and (not Path(proof["path"]).is_file() or digest(proof["path"]) != proof["sha256"]):
        raise ValueError("Observed acceptance receipt changed")
    if kernel_identity() != plan["kernel"]:
        raise ValueError("Frozen network kernel differs from the bootstrap plan")
    cfg = load_config(plan["inputs"]["config"]["path"])
    _config_checks(cfg)
    if cfg.consensus.bootstrap_iterations != plan["bootstrap_iterations"]:
        raise ValueError("Bootstrap count changed")
    psi = pd.read_csv(plan["inputs"]["matrix"]["path"], sep="\t", index_col=0)
    assigned = pd.read_csv(plan["inputs"]["assignments"]["path"], sep="\t", low_memory=False)
    samples = pd.read_csv(plan["inputs"]["samples"]["path"], sep="\t", low_memory=False)
    if len(psi) != plan["events"] or list(map(str, psi.columns)) != plan["sample_order"]:
        raise ValueError("Observed matrix shape/order changed")
    return plan, cfg, psi, assigned, samples


def partition_settings(config):
    return replace(config.network, consensus_bootstrap_iterations=0,
                   community_consensus_runs=1, bootstrap_iterations=0)


def draw_matrix(plan: dict, matrix: pd.DataFrame, iteration: int) -> pd.DataFrame:
    if not 1 <= iteration <= plan["bootstrap_iterations"]:
        raise ValueError("Bootstrap iteration outside the frozen plan")
    draw = plan["draws"][iteration - 1]
    if draw["iteration"] != iteration:
        raise ValueError("Bootstrap draw identity mismatch")
    selected = matrix.iloc[:, draw["indices"]].copy()
    selected.columns = [f"coassign_{iteration:03d}_{i:03d}" for i in range(len(draw["indices"]))]
    return selected


def verify_replicate(plan_path: str | Path, iteration: int, directory: str | Path,
                     *, require_native: bool = False) -> tuple[dict, np.ndarray]:
    plan = json.loads(Path(plan_path).read_text(encoding="utf-8"))
    root = Path(directory)
    receipt = json.loads((root / "replicate.receipt.json").read_text(encoding="utf-8"))
    if (receipt.get("schema") != "katana.network_bootstrap_replicate/v1"
            or receipt.get("status") != "completed" or receipt.get("iteration") != iteration
            or receipt.get("plan_sha256") != digest(plan_path)
            or receipt.get("draw_sha256") != plan["draws"][iteration - 1]["draw_sha256"]
            or receipt.get("kernel") != plan["kernel"]):
        raise ValueError("Bootstrap replicate proof is missing, failed, or belongs to a different plan")
    p = root / "labels.tsv"
    identity = receipt["output"]
    if p.stat().st_size != identity["bytes"] or digest(p) != identity["sha256"]:
        raise ValueError("Bootstrap replicate output changed")
    labels = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    order = sha256(json.dumps(labels.event_key.tolist()).encode()).hexdigest()
    if len(labels) != plan["events"] or order != plan["event_order_sha256"] or labels.module.eq("").any():
        raise ValueError("Bootstrap labels changed event identity or order")
    if require_native:
        native = receipt.get("native_job_id")
        if not native:
            raise ValueError("Bootstrap replicate lacks native execution identity")
        p = subprocess.run(["sacct", "-nP", "-j", str(native),
                            "--format=JobID,State,ExitCode"],
                           capture_output=True, text=True, timeout=30, check=True)
        rows = [line.split("|") for line in p.stdout.splitlines() if line.strip()]
        matches = [row for row in rows if row[0] == str(native)]
        if len(matches) != 1 or matches[0][1:3] != ["COMPLETED", "0:0"]:
            raise ValueError("Bootstrap replicate native job has not actually completed successfully")
        receipt = {**receipt, "verified_sacct": p.stdout}
    return receipt, labels.module.to_numpy(dtype=object)


def run_replicate(plan_path: str | Path, iteration: int, destination: str | Path) -> dict:
    from katana_cosplice.network import discover_module_labels
    plan, cfg, matrix, _assignments, _samples = load_plan(plan_path)
    target = Path(destination)
    if target.exists():
        receipt, _ = verify_replicate(plan_path, iteration, target)
        return {**receipt, "reused_without_recomputation": True}
    target.mkdir(parents=True, exist_ok=False)
    _json(target / "execution.claim.json", {"plan_sha256": digest(plan_path),
           "iteration": iteration, "started_utc": _now(), "pid": os.getpid()})
    receipt = {"schema": "katana.network_bootstrap_replicate/v1", "iteration": iteration,
               "plan_sha256": digest(plan_path), "draw_sha256": plan["draws"][iteration - 1]["draw_sha256"],
               "kernel": plan["kernel"], "started_utc": _now(),
               "native_job_id": (os.environ.get("SLURM_ARRAY_JOB_ID", "") + "_" + os.environ["SLURM_ARRAY_TASK_ID"]
                                 if os.environ.get("SLURM_ARRAY_TASK_ID") else os.environ.get("SLURM_JOB_ID")),
               "full_data_network_refit": False}
    try:
        bootstrap = draw_matrix(plan, matrix, iteration)
        labels = discover_module_labels(bootstrap, partition_settings(cfg),
                    event_reliability=None, reliability_floor=cfg.uncertainty.reliability_floor,
                    reliability_power=cfg.uncertainty.reliability_power,
                    gene_prior=None, gene_prior_weight=cfg.gene_prior.weight)
        if len(labels) != len(matrix):
            raise ValueError("Bootstrap discovery did not preserve the complete event universe")
        output = target / "labels.tsv"
        pd.DataFrame({"event_key": matrix.index.astype(str), "module": labels}).to_csv(output, sep="\t", index=False)
        receipt.update(status="completed", output={"path": "labels.tsv", "bytes": output.stat().st_size, "sha256": digest(output)})
    except Exception as error:
        receipt.update(status="failed", error=repr(error))
        raise
    finally:
        receipt["ended_utc"] = _now()
        _json(target / "replicate.receipt.json", receipt)
    return receipt


@contextmanager
def replay_discoveries(plan: dict, matrix: pd.DataFrame, config, labels: list[np.ndarray]):
    """Scope only the original consensus function's expensive discovery callback."""
    import katana_cosplice.consensus as consensus
    if len(labels) != plan["bootstrap_iterations"]:
        raise ValueError("Bootstrap replay must include every declared replicate")
    calls = []
    expected_settings = partition_settings(config)
    original = consensus.discover_module_labels
    def replay(bootstrap, settings, **kwargs):
        iteration = len(calls) + 1
        if iteration > len(labels) or settings != expected_settings:
            raise ValueError("Unexpected discovery request while reducing bootstrap results")
        pd.testing.assert_frame_equal(bootstrap, draw_matrix(plan, matrix, iteration), check_exact=True)
        if kwargs.get("event_reliability") is not None or kwargs.get("gene_prior") is not None:
            raise ValueError("Bootstrap replay cannot replace unmatched uncertainty or gene-prior computation")
        calls.append(iteration)
        return labels[iteration - 1].copy()
    with _REPLAY_LOCK:
        consensus.discover_module_labels = replay
        try:
            yield
            if calls != list(range(1, len(labels) + 1)):
                raise ValueError("Original consensus did not consume every planned replicate")
        finally:
            consensus.discover_module_labels = original


def reduce_replicates(plan_path: str | Path, directories: list[str | Path], *,
                      require_native: bool = True):
    from katana_cosplice.consensus import build_bootstrap_coassignment_consensus
    plan, cfg, matrix, assignments, samples = load_plan(plan_path)
    if len(directories) != plan["bootstrap_iterations"]:
        raise ValueError("Incomplete bootstrap replicate collection")
    evidence, labels = [], []
    for iteration, directory in enumerate(directories, start=1):
        receipt, result = verify_replicate(plan_path, iteration, directory, require_native=require_native)
        evidence.append(receipt); labels.append(result)
    with replay_discoveries(plan, matrix, cfg, labels):
        result = build_bootstrap_coassignment_consensus(matrix, assignments, samples,
                    cfg.network, cfg.consensus, event_reliability=None,
                    reliability_floor=cfg.uncertainty.reliability_floor,
                    reliability_power=cfg.uncertainty.reliability_power,
                    gene_prior=None, gene_prior_weight=cfg.gene_prior.weight)
    return result, {"schema": "katana.network_bootstrap_reduction/v1", "plan_sha256": digest(plan_path),
                    "full_data_network_refit": False, "replicates": evidence,
                    "kernel": plan["kernel"], "completed_utc": _now(),
                    "native_success_required": require_native}
