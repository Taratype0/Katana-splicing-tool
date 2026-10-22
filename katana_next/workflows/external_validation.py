"""Explicit external validation of frozen modules; never discovers modules."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.config import load_config
from katana_cosplice.formal_preservation import run_formal_preservation
from katana_cosplice.utils import bh_adjust
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.bulk import _plot
from katana_next.workflows.reliability import diagnostic_sources, read_matrix
from katana_next.infrastructure.fingerprints import fingerprint_file
from katana_next.workflows.sources import validate_source


def external_config(context):
    assets = [a for a in context.manifest.assets if a.kind == "cosplice_config"]
    if len(assets) != 1:
        raise ValueError("Declare exactly one CoSplice configuration with external_cohorts")
    config = load_config(assets[0].path)
    if not config.external_cohorts:
        raise ValueError("No external cohorts configured; add PSI and sample metadata before validation")
    if not config.annotation_version or config.genome_build != context.manifest.annotation.genome_build:
        raise ValueError("Declare the reference annotation version and matching genome build in CoSplice configuration")
    for cohort in config.external_cohorts:
        if cohort.genome_build != config.genome_build or cohort.annotation_version != config.annotation_version:
            raise ValueError(f"{cohort.name}: genome/annotation differs or is undeclared. No automatic coordinate conversion is performed.")
        for path in (cohort.psi_matrix, cohort.sample_metadata):
            if not path.is_file():
                raise ValueError(f"{cohort.name}: missing external input {path}")
    rscript = str(context.manifest.metadata.get("rscript_path") or config.network.wgcna_rscript or "")
    if not rscript or not Path(rscript).is_file():
        raise ValueError("Configure an existing Rscript with WGCNA and NetRep installed")
    return assets[0], config, rscript


def align_external(reference, target, samples, reference_samples, labels):
    if not reference.index.is_unique or not reference.columns.is_unique:
        raise ValueError("Reference matrix axes must be unique")
    if not {"event_key", "consensus_module"} <= set(labels) or labels.event_key.duplicated().any():
        raise ValueError("Frozen membership requires one unique event_key per row")
    if not set(labels.event_key) <= set(reference.index):
        raise ValueError("Frozen membership contains events absent from reference")
    if not target.index.is_unique or not target.columns.is_unique:
        raise ValueError("External matrix axes must be unique")
    if "sample_id" not in samples or samples.sample_id.duplicated().any():
        raise ValueError("External metadata requires unique sample_id")
    if set(target.columns) != set(samples.sample_id):
        raise ValueError("External PSI columns and sample metadata do not match exactly")
    if len(target.columns) < 4:
        raise ValueError("At least four external samples are required; this minimum is not proof of adequate power")
    for name, matrix, metadata in (("Reference", reference, reference_samples), ("External", target, samples)):
        if "sample_id" not in metadata or metadata.sample_id.duplicated().any() or set(matrix.columns) != set(metadata.sample_id):
            raise ValueError(f"{name} PSI columns and sample metadata do not match exactly")
        if "biological_sample_id" not in metadata or metadata.biological_sample_id.isna().any() or metadata.biological_sample_id.astype(str).str.strip().eq("").any():
            raise ValueError(f"{name} validation requires declared biological_sample_id; cells/technical replicates are not independent samples")
        if metadata.biological_sample_id.astype(str).str.strip().duplicated().any():
            raise ValueError(f"{name}: aggregate technical/cell measurements to biological units before external validation")
        values = matrix.to_numpy(float)
        if np.isinf(values).any() or np.any(np.isfinite(values) & ((values < 0) | (values > 1))):
            raise ValueError(f"{name} PSI must be within [0,1], with missing values left missing")
    if len(reference.columns) < 4:
        raise ValueError("At least four reference samples are required")
    if set(samples.biological_sample_id.astype(str).str.strip()) & set(reference_samples.biological_sample_id.astype(str).str.strip()):
        raise ValueError("Reference and target share biological identifiers. Confirm independence and use cohort-qualified identities.")
    usable_ref = reference.notna().sum(axis=1).ge(3) & reference.std(axis=1).gt(1e-8)
    usable_target = target.notna().sum(axis=1).ge(3) & target.std(axis=1).gt(1e-8)
    shared = reference.index[usable_ref & reference.index.isin(target.index[usable_target])]
    audit = labels.copy()
    audit["present_in_external"] = audit.event_key.isin(target.index)
    audit["usable_in_both"] = audit.event_key.isin(shared)
    audit["reference_nonmissing"] = audit.event_key.map(reference.notna().sum(axis=1))
    audit["external_nonmissing"] = audit.event_key.map(target.notna().sum(axis=1))
    audit["reference_variable"] = audit.event_key.map(reference.std(axis=1).gt(1e-8))
    audit["external_variable"] = audit.event_key.map(target.std(axis=1).gt(1e-8)).eq(True)
    audit["exclusion_reason"] = audit.apply(lambda row: ";".join(reason for test, reason in (
        (not row.present_in_external, "absent_external_event"),
        (row.reference_nonmissing < 3, "reference_low_observation_count"),
        (not row.reference_variable, "reference_no_variation"),
        (row.present_in_external and row.external_nonmissing < 3, "external_low_observation_count"),
        (row.present_in_external and not row.external_variable, "external_no_variation"),
    ) if test) or "eligible", axis=1)
    return reference.loc[shared], target.loc[shared], audit


def frozen_power(context, sources, parameters):
    explicit = int(parameters.get("soft_power", 0))
    if explicit:
        return explicit
    artifacts = context.catalog.list_artifacts(context.manifest.project_id, analysis_id="programs.cosplice")
    config = next((a for a in artifacts if a.run_fingerprint == sources["matrix__psi"].run_fingerprint
        and a.name == "Resolved CoSplice production configuration"), None)
    if config is None:
        raise ValueError("Frozen discovery configuration is unavailable. Supply its actual soft power explicitly; no value is guessed.")
    validate_source(context, config)
    return load_config(config.path).network.soft_power


def cohort_checkpoint(directory, run_fingerprint, keys):
    try:
        marker = json.loads((directory/"complete.json").read_text())
        if marker.get("run_fingerprint") != run_fingerprint:
            return None
        for key in keys:
            path = directory/f"{key}.parquet"
            if marker["files"][key] != fingerprint_file(path, force_full=True).model_dump(mode="json"):
                return None
        return {key:pd.read_parquet(directory/f"{key}.parquet") for key in keys}
    except (OSError, ValueError, KeyError, TypeError):
        return None


def frozen_membership(reference, observed, refined=None):
    """The observed network fixes the universe; refined exports may omit grey."""
    if "event_key" not in observed or observed.event_key.duplicated().any() or not set(observed.event_key) <= set(reference.index):
        raise ValueError("Observed frozen network membership is not a unique reference subset")
    selected = observed if refined is None else refined
    column = next((c for c in ("analysis_module", "module", "consensus_module") if c in selected), None)
    if column is None or selected.event_key.duplicated().any() or not set(selected.event_key) <= set(observed.event_key):
        raise ValueError("Refined membership is not a unique subset of the frozen network")
    if selected[column].isna().any() or selected[column].astype(str).str.strip().eq("").any():
        raise ValueError("Frozen module labels must be explicit; use grey for unassigned events")
    labels = selected.set_index("event_key")[[column]].rename(columns={column:"consensus_module"})
    labels = labels.reindex(observed.event_key).fillna({"consensus_module":"grey"}).rename_axis("event_key").reset_index()
    return reference.reindex(observed.event_key), labels


class ExternalPreservationWorkflow(WorkflowBase):
    version = "1.1.1"
    definition = WorkflowDefinition(
        "programs.external_validation", "programs", "Validate frozen modules externally",
        "Are these fixed modules preserved in the declared independent cohorts?",
        "Official WGCNA modulePreservation and NetRep on identical exact-event common backgrounds; module membership is frozen.",
        "Module coverage, WGCNA Z statistics, NetRep statistics and batch-wide FDR, with per-cohort checkpoints.",
        "fixed module in an external biological cohort", required_asset_kinds=("cosplice_config",),
        parameters=(
            ParameterDefinition("matrix_artifact_id", "Frozen discovery run", "artifact", "", "Uses cached discovery only; no upstream run.", options=("programs.cosplice", "2. Build PSI matrix: psi")),
            ParameterDefinition("permutations", "Permutations", "number", 1000, "Low values are exploratory; p-value resolution is limited by the permutation count.", minimum=10, maximum=10000),
            ParameterDefinition("soft_power", "NetRep signed-adjacency power", "number", 0, "Zero reads the frozen discovery configuration. An explicit value is recorded as an override, not silently inferred.", minimum=0, maximum=30),
            ParameterDefinition("min_module_size", "Minimum shared events per module", "number", 5, "Smaller modules are not tested, not classified as non-preserved.", minimum=5),
            ParameterDefinition("memory_budget_gb", "Dense-matrix safety budget (GiB)", "number", 4.0, "Conservative working-array estimate; not a guaranteed peak memory limit. No silent event cap.", minimum=.1),
            ParameterDefinition("seed", "Permutation seed", "number", 43, "Frozen before all cohort tests.", minimum=0),
        ), limitations="Exact identities and matching declared genome/annotation only. Missing cells are row-mean imputed after coverage filtering by the existing formal comparator. Independence requires trustworthy metadata. Preservation does not demonstrate phenotype or functional replication.")

    def input_paths(self, context, parameters):
        asset, config, rscript = external_config(context)
        sources = diagnostic_sources(context, parameters)
        frozen_power(context, sources, parameters)
        extras = [a.path for a in context.catalog.list_artifacts(context.manifest.project_id, analysis_id="programs.cosplice")
            if a.run_fingerprint == sources["matrix__psi"].run_fingerprint and a.name == "Resolved CoSplice production configuration"]
        return [asset.path, Path(rscript), *extras, *[p for c in config.external_cohorts for p in (c.psi_matrix, c.sample_metadata)]]

    def parent_artifacts(self, context, parameters):
        external_config(context)
        return list(diagnostic_sources(context, parameters).values())

    def run(self, context, job, run_fingerprint, parameters):
        _, config, rscript = external_config(context)
        bound = parameters.get("source_artifact_ids", {}).get("programs.cosplice::2. Build PSI matrix: psi")
        sources = diagnostic_sources(context, {**parameters, "matrix_artifact_id": bound or parameters.get("matrix_artifact_id", "")})
        power = frozen_power(context, sources, parameters)
        reference = read_matrix(sources["matrix__psi"])
        reference_samples = pd.read_parquet(sources["matrix__samples"].path)
        raw_matrix_events = len(reference)
        # The null background is the frozen network universe, including grey,
        # not every unfiltered rMATS event that never entered this network.
        observed = pd.read_parquet(sources["network_observed__assignments"].path)
        refined = (pd.read_parquet(sources["consensus_modules__assignments"].path)
                   if "consensus_modules__assignments" in sources else None)
        reference, labels = frozen_membership(reference, observed, refined)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"source_discovery_run": sources["matrix__psi"].run_fingerprint,
                             "method": self.definition.method, "limitations": self.definition.limitations})
        parents = tuple(a.artifact_id for a in sources.values())
        summaries = {"wgcna": [], "netrep": [], "audit": [], "coverage": []}
        for number, cohort in enumerate(config.external_cohorts):
            job.check_cancelled()
            job.report(.05+.8*number/len(config.external_cohorts), f"Validating frozen modules in {cohort.name}; discovery is not rerun")
            checkpoint = writer.directory / f"cohort_{number:03d}"
            checkpoint.mkdir(parents=True, exist_ok=True)
            marker = checkpoint / "complete.json"
            cached = cohort_checkpoint(checkpoint, run_fingerprint, summaries)
            if cached is not None:
                job.report(.05+.8*number/len(config.external_cohorts), f"Checkpoint hit: {cohort.name}")
                for key in summaries:
                    summaries[key].append(cached[key])
                continue
            target = pd.read_csv(cohort.psi_matrix, sep="\t", index_col=0)
            target.index = target.index.astype(str)
            samples = pd.read_csv(cohort.sample_metadata, sep="\t", dtype=str)
            ref, test, coverage = align_external(reference, target, samples, reference_samples, labels)
            estimate = 16*8*len(ref)**2 / 1024**3
            if estimate > float(parameters["memory_budget_gb"]):
                raise ValueError(f"{cohort.name}: dense-array safety estimate {estimate:.2f} GiB exceeds budget. No events were silently removed.")
            # The full frozen denominator stays in coverage; the formal kernel
            # receives only the aligned universe, with module labels unchanged.
            aligned_labels = labels[labels.event_key.isin(ref.index)].reset_index(drop=True)
            fit = run_formal_preservation(ref, test, aligned_labels, cohort=cohort.name, arm="Frozen CoSplice", rscript=rscript,
                permutations=int(parameters["permutations"]), soft_power=int(power),
                min_module_size=int(parameters["min_module_size"]), random_seed=int(parameters["seed"]),
                temp_root=checkpoint / "runtime")
            fit.audit["dense_array_estimate_gib"] = estimate
            fit.audit["reference_samples"] = len(reference.columns)
            fit.audit["source_raw_matrix_events"] = raw_matrix_events
            fit.audit["frozen_network_background_events"] = len(reference)
            fit.audit["external_biological_units"] = samples.biological_sample_id.nunique()
            fit.audit["small_cohort_warning"] = min(len(reference.columns), len(test.columns)) < 15
            fit.audit["low_permutation_warning"] = int(parameters["permutations"]) < 1000
            fit.audit["reference_imputed_cells"] = int(ref.isna().sum().sum())
            fit.audit["external_imputed_cells"] = int(test.isna().sum().sum())
            fit.wgcna["cohort"] = cohort.name
            fit.netrep["cohort"] = cohort.name
            coverage["cohort"] = cohort.name
            current = {"wgcna": fit.wgcna, "netrep": fit.netrep, "audit": fit.audit, "coverage": coverage}
            for key, frame in current.items():
                frame.to_parquet(checkpoint / f"{key}.parquet", index=False)
                summaries[key].append(frame)
            context.results.write_plot_spec({"run_fingerprint": run_fingerprint, "completed": True,
                "files": {key:fingerprint_file(checkpoint/f"{key}.parquet",force_full=True).model_dump(mode="json") for key in current}}, marker)
        combined = {key: pd.concat(frames, ignore_index=True) for key, frames in summaries.items()}
        netrep = combined["netrep"]
        pcols = [c for c in netrep if c.startswith("pvalue_")]
        if pcols:
            adjusted = bh_adjust(netrep[pcols].to_numpy(float).ravel()).reshape(len(netrep),len(pcols))
            for i, column in enumerate(pcols):
                netrep[column.replace("pvalue_", "batch_FDR_")] = adjusted[:,i]
        for key, frame in combined.items():
            writer.table(f"External preservation {key}", f"{key}.parquet", frame,
                entity_type="module_validation", parents=parents)
        coverage = combined["coverage"].groupby(["cohort", "consensus_module"]).agg(
            reference_events=("event_key", "size"), usable_events=("usable_in_both", "sum")).reset_index()
        coverage["usable_fraction"] = coverage.usable_events/coverage.reference_events
        writer.table("External module coverage denominators", "module_coverage.parquet", coverage, entity_type="module_validation", parents=parents)
        wgcna = combined["wgcna"]
        zcolumn = next((c for c in wgcna if c.startswith("Zsummary")), None)
        if zcolumn:
            display = wgcna[["cohort", "module", zcolumn]].rename(columns={zcolumn: "preservation_zsummary"})
            writer.plot("Official WGCNA module preservation", "wgcna_preservation.vl.json", _plot(display.to_dict("records"), "rect",
                {"x": {"field": "cohort", "type": "nominal"}, "y": {"field": "module", "type": "nominal"},
                 "color": {"field": "preservation_zsummary", "type": "quantitative", "title": "WGCNA preservation Zsummary"},
                 "tooltip": [{"field": c} for c in display]}, title="Are fixed modules preserved externally?",
                description="Official permutation Zsummary, not a p-value. Interpret together with module coverage, sample size and NetRep tests."), entity_type="module_validation", parents=parents)
        writer.plot("External module event coverage", "coverage.vl.json", _plot(coverage.to_dict("records"), "rect",
            {"x": {"field": "cohort", "type": "nominal"}, "y": {"field": "consensus_module", "type": "nominal"},
             "color": {"field": "usable_fraction", "type": "quantitative", "scale": {"domain": [0,1]}},
             "tooltip": [{"field": c} for c in coverage]}, title="How much of each frozen module could be tested?",
            description="Coverage is not preservation significance; inspect WGCNA and NetRep separately."), entity_type="module_validation", parents=parents)
        return writer.finish()
