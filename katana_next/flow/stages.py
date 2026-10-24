"""Headless adapters for computations required by the fixed study report.

The report renderer never calls these functions.  The flow scheduler gives each
adapter its own output directory and only accepts its verified output ledger.
Original algorithms, source projects and previous results are read-only.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from threading import RLock
from typing import Any
from unittest.mock import patch


ADAPTERS = frozenset({"reference.bed12", "qc.genebody", "qc.sequencing", "as.input_projection", "network.fit",
                      "network.bootstrap", "network.associate", "network.annotate",
                      "mapp.prepare", "mapp.maei"})
_PLOT_LOCK = RLock()
_MAPP_SCRIPTS = (
    "modules/MAEI/scripts/optimize-MARA-model-parameters.R",
    "modules/MAEI/scripts/model.cpp",
    "modules/MAEI/scripts/optim.c",
    "modules/MAEI/scripts/zeroin.c",
    "modules/MAEI/scripts/calculate-Zscores.py",
    "modules/MAEI/scripts/extend-tables.py",
    "modules/MAEI/scripts/calculate_Zscores_statistical_significancy.R",
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _executable(value: Any, default: str) -> str:
    from .container_runtime import require_image_path
    name = str(value or default)
    executable = shutil.which(name)
    if executable:
        # Preserve a venv's lexical bin/python path: resolving its symlink
        # before invocation would select the base environment instead.
        return require_image_path(Path(executable).absolute(), "executable")
    candidate = Path(name).expanduser()
    if candidate.is_file() and (os.name == "nt" or os.access(candidate, os.X_OK)):
        return require_image_path(candidate.absolute(), "executable")
    raise ValueError(f"Executable is unavailable: {name}")


def _r_environment(rscript, parameters):
    from .container_runtime import sanitize_container_environment, validate_runtime_settings, require_image_path
    sanitize_container_environment()
    require_image_path(rscript, "Rscript")
    validate_runtime_settings(parameters)
    from katana_cosplice.benchmark import _r_subprocess_env
    env = _r_subprocess_env(rscript)
    libraries = parameters.get("r_libraries")
    if libraries:
        if not isinstance(libraries, list) or not all(isinstance(x,str) and Path(x).is_dir() for x in libraries):
            raise ValueError("r_libraries must list existing R package library directories")
        env["R_LIBS_USER"] = os.pathsep.join(str(Path(x).resolve()) for x in libraries)
    elif parameters.get("r_library"):
        if not Path(parameters["r_library"]).is_dir():
            raise ValueError("Declared R package library does not exist")
        env["R_LIBS_USER"] = str(Path(parameters["r_library"]).resolve())
    return env


@lru_cache(maxsize=16)
def _r_packages(rscript, packages, environment):
    # Cached within one preflight process; all tasks sharing the same actual R
    # executable and library environment need only one non-mutating probe.
    quoted = ",".join("'"+p+"'" for p in packages)
    expression = ("p<-c("+quoted+");missing<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];"
                  "cat('KATANA_MISSING_R_PACKAGES=',paste(missing,collapse=','),sep='')")
    result = subprocess.run([rscript,"--vanilla","-e",expression], capture_output=True,
                            text=True, env=dict(environment), timeout=30, check=False)
    marker = "KATANA_MISSING_R_PACKAGES="
    if result.returncode or marker not in result.stdout:
        raise ValueError("R dependency probe failed: " + result.stderr.strip())
    missing = result.stdout.split(marker,1)[1].strip()
    if missing:
        raise ValueError("Required R packages are unavailable in the declared runtime: " + missing)


def _require_r_packages(rscript, packages, parameters):
    env = _r_environment(rscript, parameters)
    _r_packages(rscript, tuple(packages), tuple(sorted(env.items())))


@contextmanager
def _network_r_libraries(parameters):
    """One isolated headless process; preserve the caller's R library setting."""
    from .container_runtime import validate_runtime_settings
    validate_runtime_settings(parameters)
    if not parameters.get("r_library") and not parameters.get("r_libraries"):
        yield
        return
    previous = os.environ.get("R_LIBS_USER")
    libraries = parameters.get("r_libraries") or [parameters["r_library"]]
    with _PLOT_LOCK:
        os.environ["R_LIBS_USER"] = os.pathsep.join(str(Path(p).resolve()) for p in libraries)
        try:
            yield
        finally:
            if previous is None:
                os.environ.pop("R_LIBS_USER",None)
            else:
                os.environ["R_LIBS_USER"] = previous


def _adapter(spec: dict) -> str:
    return str(spec.get("adapter") or spec.get("workflow_id", "").removeprefix("flow."))


def readiness(spec: dict) -> list[str]:
    """Validate executable prerequisites without creating output or submitting work.

    Computed dependency paths can be declared later by the scheduler.  Scientific
    table contents are validated again immediately before native execution.
    """
    from .bundled_runtime import resolve_bundled_stage
    from .container_runtime import in_container, validate_task_runtime_boundaries
    try:
        project_path = spec.get("inputs", {}).get("project")
        project = json.loads(Path(project_path).read_text(encoding="utf-8-sig")) if in_container() and project_path else None
        validate_task_runtime_boundaries(spec, project=project, project_path=project_path)
        spec = resolve_bundled_stage(spec)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return [str(exc)]
    adapter = _adapter(spec)
    issues: list[str] = []
    inputs = spec.get("inputs") or {}
    parameters = spec.get("parameters") or {}
    deferred = set(spec.get("readiness_deferred_inputs") or [])
    deferred_parameters = set(spec.get("readiness_deferred_parameters") or [])
    if adapter == "qc.sequencing":
        if "project" in deferred:
            return []
        from .sequencing import readiness as sequencing_readiness
        return sequencing_readiness(spec)
    if adapter == "mapp.maei" and inputs.get("applicability_proof"):
        from .mapp_inputs import validate_inapplicability
        try:
            validate_inapplicability(inputs["applicability_proof"], sha=_sha, scope=spec.get("scope"))
        except (ValueError, OSError, KeyError, TypeError) as exc:
            return [str(exc)]
        return []
    from .dependency_checks import python_dependency_issues
    issues.extend(python_dependency_issues(spec, resolve_executable=_executable))
    if issues:
        # In particular, mapp_inputs imports Bio via isoform_consequences.
        # Report a broken runtime before importing that adapter or touching
        # scientific inputs; do not fall through to a later missing-Bio crash.
        return issues
    required = {
        "as.input_projection": ("source", "project"),
        "reference.bed12": ("gtf",), "qc.genebody": ("bam", "reference_bed"),
        "network.fit": ("config",), "network.bootstrap": ("config",),
        "network.associate": ("config",), "network.annotate": ("config",),
        "mapp.prepare": ("native_root", "events", "gtf", "genome_fasta"),
        "mapp.maei": ("native_root", "inclusion_table", "sitecount_matrix", "design"),
    }
    if adapter not in required:
        return [f"No scientific compute adapter is implemented for {adapter!r}"]
    if adapter.startswith("network.") and not inputs.get("config"):
        required[adapter] = ("project",)
    for name in required[adapter]:
        if name in deferred:
            continue
        value = inputs.get(name)
        if not value:
            issues.append(f"{adapter}: declare inputs.{name}")
        elif name == "native_root":
            from .mapp_inputs import COUNTER
            for relative in ((COUNTER,) if adapter == "mapp.prepare" else _MAPP_SCRIPTS):
                if not (Path(value) / relative).is_file():
                    issues.append(f"Missing official MAPP source: {relative}")
        elif not Path(value).is_file():
            issues.append(f"{adapter}: missing {name}: {value}")
    try:
        if adapter == "network.fit" and not issues:
            rscript = parameters.get("rscript")
            cfg = None
            if inputs.get("config"):
                from katana_cosplice.config import load_config
                cfg = load_config(inputs["config"])
                cfg.validate()
                rscript = rscript or cfg.network.wgcna_rscript
            if inputs.get("project") and "project" not in deferred and "network_scope" in parameters:
                from katana_next.infrastructure.manifest_store import ManifestStore
                from .network_inputs import network_scope_identity, resolve_network_rmats_mode, validate_external_network_config
                manifest = ManifestStore().load(inputs["project"])[1]
                identity = network_scope_identity(manifest, parameters)
                mode = resolve_network_rmats_mode(manifest, parameters)
                if cfg is not None:
                    if cfg.rmats_mode != mode:
                        raise ValueError("External network configuration mode differs from the resolved task mode")
                    validate_external_network_config(cfg, manifest, identity, _sha)
            if not rscript and inputs.get("project"):
                if "project" in deferred:
                    # A raw resource template deliberately has no fabricated
                    # samples or outputs before its genuine upstream producer.
                    raw = json.loads(Path(inputs["project"]).read_text(encoding="utf-8"))
                    rscript = raw.get("metadata", {}).get("rscript_path")
                else:
                    from katana_next.infrastructure.manifest_store import ManifestStore
                    rscript = ManifestStore().load(inputs["project"])[1].metadata.get("rscript_path")
            rscript = _executable(rscript, "Rscript")
            _require_r_packages(rscript,("WGCNA",),parameters)
        if adapter == "qc.genebody":
            _executable(parameters.get("rseqc_executable"), "geneBody_coverage.py")
            if not str(parameters.get("sample_id") or spec.get("scope", {}).get("sample_id", "")).strip():
                issues.append("qc.genebody requires an explicit sample_id")
        if adapter == "mapp.maei":
            rscript = _executable(parameters.get("rscript"), "Rscript")
            if not issues:
                _require_r_packages(rscript, ("Rcpp", "optparse"), parameters)
            _executable(parameters.get("python"), sys.executable)
            if not inputs.get("source_sha256"):
                issues.append("mapp.maei requires inputs.source_sha256 for all seven official source files")
            elif inputs.get("native_root"):
                if set(inputs["source_sha256"]) != set(_MAPP_SCRIPTS):
                    issues.append("Bind exactly all seven official MAEI source files")
                else:
                    for relative, digest in inputs["source_sha256"].items():
                        path = Path(inputs["native_root"]) / relative
                        if path.is_file() and _sha(path) != digest:
                            issues.append(f"Official MAPP source fingerprint changed: {relative}")
            if not {"inclusion_table", "sitecount_matrix", "design", "quantity_provenance"} & deferred:
                if all(inputs.get(k) and Path(inputs[k]).is_file() for k in ("inclusion_table", "sitecount_matrix", "design")):
                    _maei_tables(inputs)
        if adapter == "mapp.prepare":
            from .mapp_inputs import window_id
            _executable(parameters.get("python"), sys.executable)
            if not inputs.get("counter_sha256"):
                issues.append("mapp.prepare requires inputs.counter_sha256")
            elif inputs.get("native_root"):
                from .mapp_inputs import COUNTER
                counter = Path(inputs["native_root"]) / COUNTER
                if counter.is_file() and _sha(counter) != inputs["counter_sha256"]:
                    issues.append("Official MAPP kmer counter source fingerprint changed")
            if "sample_ids" not in deferred_parameters and not parameters.get("sample_ids"):
                issues.append("mapp.prepare requires explicit parameters.sample_ids")
            if "site" not in deferred_parameters and parameters.get("site") not in {"3ss", "5ss"}:
                issues.append("mapp.prepare requires explicit parameters.site (3ss/5ss)")
            if {"window_start","window_end"} & deferred_parameters:
                pass
            elif "window_start" not in parameters or "window_end" not in parameters:
                issues.append("mapp.prepare requires the predeclared window_start/window_end")
            else:
                window_id(int(parameters["window_start"]), int(parameters["window_end"]))
            if "quantification" not in deferred and not inputs.get("quantification"):
                issues.append("mapp.prepare requires sample-mapped existing transcript quantification")
            for sample, value in ({} if "quantification" in deferred else inputs.get("quantification", {})).items():
                if not Path(value).is_file():
                    issues.append(f"{sample}: missing existing transcript quantification {value}")
            if "genome_fasta" not in deferred and inputs.get("genome_fasta") and not Path(str(inputs["genome_fasta"])+".fai").is_file():
                issues.append("mapp.prepare requires the matched genome FASTA index")
    except (ValueError, OSError, KeyError, TypeError, subprocess.TimeoutExpired) as exc:
        issues.append(str(exc))
    return issues


def _run(command: list[str], out: Path, name: str, *, env: dict | None = None) -> None:
    """A failed or killed native command never receives a successful stage receipt."""
    with (out / f"{name}.stdout.log").open("xb") as stdout, (out / f"{name}.stderr.log").open("xb") as stderr:
        result = subprocess.run(command, cwd=out, stdout=stdout, stderr=stderr,
                                env=env, check=False)
    (out / f"{name}.command.json").write_text(json.dumps({
        "argv": command, "returncode": result.returncode,
    }, indent=2), encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"{name} exited {result.returncode}; see {out / (name + '.stderr.log')}")


def _ledger(paths: list[Path], adapter: str) -> list[dict]:
    records = []
    for path in dict.fromkeys(paths):
        path = path.resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Required {adapter} output is missing or empty: {path}")
        metadata = _output_metadata(path, adapter)
        records.append({"path": str(path), "size_bytes": path.stat().st_size,
                        "sha256": _sha(path), "name": path.name,
                        "entity_type": "flow_compute_output", "analysis_id": "flow." + adapter,
                        **metadata})
    if not records:
        raise ValueError(f"{adapter} produced no verified output")
    return records


def _output_metadata(path: Path, adapter: str) -> dict:
    if adapter.startswith("network."):
        from katana_cosplice.stage_registry import stage_spec
        from katana_next.workflows.programs import _entity_type
        keys = ("matrix", "stable_qc", "network_observed", "network_residual",
                "consensus_modules", "associate_observed", "associate_residual", "annotate_modules")
        scope = "full_network_with_bootstrap" if adapter == "network.bootstrap" else "observed_network_baseline_only"
        extra = {"completion_scope": scope}
        if path.suffix.lower() == ".png":
            prefix = path.parent.parent
            if "dendrogram" in path.name:
                data = prefix / "event_linkage.npy"
            elif "program_score_heatmap" in path.name:
                data = prefix / "module_program_scores.tsv"
            elif path.name == "bootstrap_coassignment_heatmap.png":
                data = prefix / "bootstrap_coassignment_matrix.npz"
            elif path.name == "consensus_module_sizes.png":
                data = prefix / "consensus_module_summary.tsv"
            elif "consensus" in str(prefix):
                data = prefix / "consensus_event_module_assignments.tsv"
            elif "association" in str(prefix):
                data = prefix / "module_trait_association.tsv"
            elif "annotation" in str(prefix):
                data = prefix / "module_evidence_summary.tsv"
            else:
                data = prefix / "event_module_assignments.tsv"
            if data.is_file():
                extra["data_links"] = [{"role": "plot_data", "label": "Complete numerical plot source",
                                        "path": str(data), "contrast_id": "all", "policy": None,
                                        "unit": "network event/module"}]
        metadata = {"analysis_id": "programs.cosplice", "entity_type": "cosplice_provenance",
                    "metadata": extra}
        normalized = path.as_posix()
        for key in keys:
            stage = stage_spec(key)
            for item in stage.artifacts:
                if normalized.endswith("/" + Path(item.relative_path).as_posix()):
                    return {**metadata, "name": f"{stage.title}: {item.key.replace('_', ' ')}",
                            "entity_type": _entity_type(item.key), "stage": key,
                            "artifact_key": item.key}
        if path.name == "event_linkage.npy":
            return {**metadata, "name": "Observed network: full numerical event linkage",
                    "entity_type": "network_linkage"}
        if path.name in {"module_hub_events.tsv", "module_hub_genes.tsv"}:
            return {**metadata, "name": "Observed network: " + path.name.removesuffix(".tsv").replace("_", " "),
                    "entity_type": "network_hub_events" if path.name == "module_hub_events.tsv" else "network_hub_genes"}
        return metadata
    if adapter == "qc.genebody":
        metadata = {}
        if path.suffix == ".png":
            metadata["data_links"] = [{"role": "plot_data", "label": "Measured gene-body coverage profile",
                                       "path": str(path.with_name("gene_body_coverage.tsv")),
                                       "contrast_id": "all", "policy": None, "unit": "sample/percentile"}]
        return {"analysis_id": "quality.read_coverage_bias", "entity_type": "coverage_profile", "metadata": metadata,
                "name": {"gene_body_coverage.tsv": "Measured RSeQC gene-body coverage profiles",
                         "gene_body_coverage.png": "Measured RSeQC 5-prime to 3-prime coverage"}.get(path.name, path.name)}
    if adapter == "mapp.maei":
        return {"analysis_id": "prioritize.mapp_native", "entity_type": "mapp_native",
                "name": "Official MAPP MAEI: " + path.name}
    return {}


def execute_stage(spec: dict, output_dir: str | Path) -> dict:
    """Execute one independent native stage and return its complete file ledger.

    The caller owns cancellation and scheduling.  Existing output directories are
    rejected: resume is performed using a verified scheduler receipt, never by
    assuming that partially written files imply success.
    """
    from .bundled_runtime import resolve_bundled_stage
    spec = resolve_bundled_stage(spec)
    if spec.get("readiness_deferred_inputs") or spec.get("readiness_deferred_parameters"):
        raise ValueError("Resolve dependency inputs before executing a scientific stage")
    if _adapter(spec) == "qc.sequencing":
        from .sequencing import execute_stage as sequence_stage
        return sequence_stage(spec,output_dir)
    problems = readiness(spec)
    if problems:
        raise ValueError("; ".join(problems))
    out = Path(output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Reconcile the existing stage output before retrying: {out}")
    out.mkdir(parents=True, exist_ok=True)
    adapter = _adapter(spec)
    if adapter == "reference.bed12":
        evidence, paths = _gtf_bed12(spec, out)
    elif adapter == "as.input_projection":
        from .eligibility_projection import project_existing_as
        evidence = project_existing_as(spec['inputs']['source'], spec['inputs']['project'],
            out / 'projection', contrast_ids=spec.get('scope', {}).get('contrast_ids', []))
        paths = sorted((out / 'projection').iterdir())
    elif adapter == "qc.genebody":
        evidence, paths = _gene_body(spec, out)
    elif adapter.startswith("network."):
        evidence, paths = _network(spec, out)
    elif adapter == "mapp.prepare":
        from .mapp_inputs import prepare_inputs
        evidence, paths = prepare_inputs(spec, out, run=_run, sha=_sha, executable=_executable)
    elif adapter == "mapp.maei":
        if spec.get("inputs", {}).get("applicability_proof"):
            from .mapp_inputs import validate_inapplicability
            proof, source_audit = validate_inapplicability(spec["inputs"]["applicability_proof"], sha=_sha, scope=spec.get("scope"))
            audit = out / "all_event_eligibility.tsv"
            shutil.copyfile(source_audit, audit)
            source_proof = Path(spec["inputs"]["applicability_proof"])
            proof = {**proof, "eligibility_table": str(audit),
                     "source_applicability_proof": {"path": str(source_proof.resolve()), "sha256": _sha(source_proof)}}
            proof_path = out / "native_maei_applicability.json"
            proof_path.write_text(json.dumps(proof, indent=2), encoding="utf-8")
            evidence, paths = {"status": "not_applicable", "native_fitting_performed": False,
                               "applicability_proof": str(proof_path)}, [audit, proof_path]
        else:
            evidence, paths = _maei(spec, out)
    else:  # readiness has already rejected unsupported adapters.
        raise ValueError(adapter)
    records = _ledger(paths, adapter)
    if adapter.startswith("network."):
        for record in records:
            record.setdefault("metadata", {}).update({
                key: evidence[key] for key in ("bootstrap_enabled", "stability_status", "hub_metric") if key in evidence})
    identity = evidence.get("network_identity") if adapter.startswith("network.") else None
    if identity:
        for record in records:
            metadata = record.setdefault("metadata", {})
            metadata.update(network_identity=identity, interpretation_status="minimum_met")
            if identity["network_scope"] == "contrast":
                metadata["data_links"] = [{**link, "contrast_id": identity["contrast_ids"][0]}
                                          for link in metadata.get("data_links", [])]
    if adapter == 'as.input_projection':
        for record in records:
            if Path(record['path']).name == 'measured_as_input_projection.parquet':
                record.update(analysis_id='flow.as_input_projection', kind='table',
                    name='Measured AS input with current RS eligibility', entity_type='splice_event',
                    metadata={'read_support_policy': evidence['read_support_policy'],
                        'scientific_model_fitted': False, 'native_statistics_reestimated': False,
                        'source_result_modified': False})
    receipt = {"schema": "katana.flow.stage/v1", "adapter": adapter,
               "status": "not_applicable" if evidence.get("status") == "not_applicable" else "complete",
               "outputs": records, "evidence": evidence}
    (out / "stage.receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    return receipt


def _gtf_bed12(spec: dict, out: Path) -> tuple[dict, list[Path]]:
    """Use matched annotation exon coordinates; GTF 1-based -> BED 0-based."""
    import gzip
    import re
    source = Path(spec["inputs"]["gtf"])
    transcripts: dict[str, dict] = {}
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt", encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF row {number}")
            if fields[2] != "exon":
                continue
            attributes = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
            identity = attributes.get("transcript_id")
            if not identity or fields[6] not in {"+", "-"}:
                raise ValueError(f"Exon row {number} requires transcript_id and strand")
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid exon interval at GTF row {number}")
            record = transcripts.setdefault(identity, {"chrom": fields[0], "strand": fields[6], "exons": set()})
            if (record["chrom"], record["strand"]) != (fields[0], fields[6]):
                raise ValueError(f"Transcript {identity} has conflicting chromosomes or strands")
            record["exons"].add((start, end))
    if not transcripts:
        raise ValueError("Matched annotation contains no strand-oriented transcript exons")
    destination = out / "reference_transcripts.bed"
    with destination.open("x", encoding="utf-8", newline="\n") as stream:
        for identity, item in sorted(transcripts.items()):
            exons = sorted(item["exons"])
            if any(right[0] < left[1] for left, right in zip(exons, exons[1:])):
                raise ValueError(f"Overlapping exons in transcript {identity}; annotation needs explicit resolution")
            start, end = exons[0][0], exons[-1][1]
            fields = [item["chrom"], start, end, identity, 0, item["strand"], start, end, 0,
                      len(exons), ",".join(str(b-a) for a,b in exons),
                      ",".join(str(a-start) for a,b in exons)]
            stream.write("\t".join(map(str, fields)) + "\n")
    return {"gtf_sha256": _sha(source), "transcripts": len(transcripts),
            "coordinate_conversion": "GTF 1-based closed -> BED 0-based half-open"}, [destination]


def _gene_body(spec: dict, out: Path) -> tuple[dict, list[Path]]:
    import pandas as pd
    from katana_next.science.coverage_bias import read_gene_body_coverage
    inputs, parameters = spec["inputs"], spec.get("parameters") or {}
    sample = str(parameters.get("sample_id") or spec.get("scope", {}).get("sample_id", ""))
    # Keep the owned staged BAM path: its adjacent index may be a separate
    # symlink to another immutable upstream task directory. Resolving the BAM
    # would silently lose that valid adjacent index before native RSeQC opens it.
    bam, bed = Path(inputs["bam"]).expanduser().absolute(), Path(inputs["reference_bed"]).resolve()
    bed_chroms: set[str] = set()
    rows = 0
    with bed.open(encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 12 or fields[5] not in {"+", "-"}:
                raise ValueError(f"Strand-oriented BED12 is required at row {number}")
            start, end, count = int(fields[1]), int(fields[2]), int(fields[9])
            sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
            offsets = [int(x) for x in fields[11].rstrip(",").split(",")]
            if (start < 0 or end <= start or count < 1 or len(sizes) != count or len(offsets) != count
                    or min(sizes) < 1 or min(offsets) < 0
                    or any(offset + size > end-start for offset,size in zip(offsets,sizes))):
                raise ValueError(f"Invalid BED12 transcript at row {number}")
            bed_chroms.add(fields[0]); rows += 1
    if not rows:
        raise ValueError("Reference BED12 is empty")
    try:
        import pysam
        handle = pysam.AlignmentFile(str(bam), "rb")
        header = handle.header.to_dict()
    except ImportError:
        import bamnostic
        index = next((p for p in (Path(str(bam) + ".bai"), bam.with_suffix(".bai"),
                                  Path(str(bam) + ".csi")) if p.is_file()), None)
        if index is None:
            raise ValueError("Gene-body QC requires a valid BAM index")
        handle = bamnostic.AlignmentFile(str(bam), "rb", filepath_index=str(index), require_index=True)
        header = dict(handle.header)
        if isinstance(header.get("HD"), list):
            header["HD"] = header["HD"][0]
    try:
        if header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError("Gene-body QC requires a coordinate-sorted BAM")
        if not handle.check_index():
            raise ValueError("Gene-body QC requires a valid BAM index")
        missing = bed_chroms - set(handle.references)
        if missing:
            raise ValueError(f"BED/BAM reference chromosomes differ: {sorted(missing)[:8]}")
    finally:
        handle.close()
    minimum = int(parameters.get("min_transcript_length", 100))
    if minimum < 100:
        raise ValueError("RSeQC requires min_transcript_length of at least 100 bases")
    command = [_executable(parameters.get("rseqc_executable"), "geneBody_coverage.py"),
               "-i", str(bam), "-r", str(bed), "-o", str(out / "gene_body"), "-l", str(minimum)]
    _run(command, out, "RSeQC_gene_body")
    profile = out / "gene_body.geneBodyCoverage.txt"
    frame = read_gene_body_coverage(profile)
    if frame.sample_id.nunique() != 1:
        raise ValueError("One per-sample RSeQC task must emit exactly one profile")
    native_sample = str(frame.sample_id.iloc[0])
    frame["sample_id"] = sample
    if frame.coverage.sum() <= 0:
        raise ValueError(f"No measurable gene-body coverage for {sample}")
    destination = out / "gene_body_coverage.tsv"
    frame.to_csv(destination, sep="\t", index=False)
    # Render from the measured native profile, independent of RSeQC's R graphics.
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(frame.percentile, frame.coverage / frame.coverage.mean(), color="#2B64A5")
    ax.set(xlabel="Gene body: 5′ to 3′ percentile", ylabel="Coverage / sample mean", title=sample)
    fig.tight_layout()
    plot = out / "gene_body_coverage.png"
    fig.savefig(plot, dpi=160); plt.close(fig)
    return {"sample_id": sample, "native_sample_label": native_sample,
            "bam": str(bam), "bam_resolved_source": str(bam.resolve()),
            "reference_bed_sha256": _sha(bed),
            "reference_transcripts": rows, "min_transcript_length": minimum,
            "measurement": "RSeQC gene-body coverage; not FastQC per-base quality",
            "project_asset": {"kind": "gene_body_coverage", "path": str(destination),
                              "metadata": {"sample_id": sample}}}, [
                destination, profile, plot, out / "RSeQC_gene_body.command.json",
            ]


def _iterative_tree(linkage, output: Path, name: str) -> Path:
    """Plot exact linkage geometry without recursive traversal of a deep tree."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    matrix = np.asarray(linkage, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != 4 or len(matrix) == 0 or not np.isfinite(matrix).all():
        raise ValueError("Expected a finite scipy linkage matrix with n-1 rows")
    leaves = len(matrix) + 1
    root = 2 * leaves - 2
    order: list[int] = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node < leaves:
            order.append(node)
        else:
            left, right = map(int, matrix[node-leaves, :2])
            if left < 0 or right < 0 or left >= node or right >= node:
                raise ValueError("Linkage child does not precede its parent")
            stack.extend([right, left])
    if sorted(order) != list(range(leaves)):
        raise ValueError("Linkage does not cover every event exactly once")
    x = {node: 5.0 + 10.0 * i for i, node in enumerate(order)}
    heights = {node: 0.0 for node in order}
    segments = []
    for index, row in enumerate(matrix):
        node = leaves + index
        left, right = map(int, row[:2])
        height = float(row[2])
        segments.append([(x[left], heights[left]), (x[left], height),
                         (x[right], height), (x[right], heights[right])])
        x[node] = (x[left] + x[right]) / 2; heights[node] = height
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.add_collection(LineCollection(segments, colors="#2B64A5", linewidths=0.5))
    ax.set(xlim=(0, 10*leaves), ylim=(0, max(heights.values())*1.05 or 1),
           xlabel="Events in saved linkage leaf order", ylabel="Correlation distance",
           title=f"{name.title()} event dendrogram", xticks=[])
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180); plt.close(fig)
    return output


def _ordinal_index_compat(rankdata):
    """Restore integer index transport; do not change any rank or tie order.

    SciPy 1.18 changed rankdata to always return floating point, including
    ordinal ranks. dynamicTreeCut 0.1.1 indexes NumPy arrays with these ranks.
    https://docs.scipy.org/doc/scipy/release/1.18.0-notes.html
    """
    import numpy as np
    def compatible(values, method="average", **kwargs):
        ranked = rankdata(values,method=method,**kwargs)
        if method == "ordinal":
            array = np.asarray(ranked)
            if not np.isfinite(array).all() or not np.equal(array,np.floor(array)).all():
                raise ValueError("Dynamic Tree Cut ordinal indices must remain finite exact integers")
            return array.astype(np.intp,copy=False)
        return ranked
    return compatible


@contextmanager
def persist_network_trees():
    """Save numerical linkage before graphics; never change the fit or modules.

    Scoped to one headless task process.  A process-level lock also protects users
    who call this wrapper from threads.  The original renderer handles the other
    figures unchanged; the dendrogram uses the exact saved linkage iteratively.
    """
    import numpy as np
    import pandas as pd
    import katana_cosplice.pipeline as pipeline
    import importlib
    dynamic = importlib.import_module("dynamicTreeCut.dynamicTreeCut")
    original = []

    def render(result, samples, output_dir, network_name):
        output_dir = Path(output_dir)
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        tree = result.linkage_matrix
        if tree is not None:
            np.save(output_dir.parent / "event_linkage.npy", tree, allow_pickle=False)
            pd.DataFrame({"leaf_id": range(len(result.selected_matrix)),
                          "event_key": result.selected_matrix.index.astype(str)}).to_csv(
                              output_dir.parent / "event_linkage_leaf_ids.tsv", sep="\t", index=False)
        plots = original[0](replace(result, linkage_matrix=None), samples, output_dir, network_name)
        if tree is not None:
            plots.append(_iterative_tree(tree, output_dir / f"{network_name}_event_dendrogram.png", network_name))
        return plots

    with _PLOT_LOCK:
        original.append(pipeline.render_network_plots)
        with patch.object(pipeline, "render_network_plots", render), \
                patch.object(dynamic,"rankdata",_ordinal_index_compat(dynamic.rankdata)):
            yield


def _network(spec: dict, out: Path) -> tuple[dict, list[Path]]:
    from katana_cosplice.config import load_config
    from katana_cosplice.pipeline import CoSplicePipeline
    from katana_next.workflows.programs import write_production_configuration
    inputs, parameters = spec["inputs"], spec.get("parameters") or {}
    adapter = _adapter(spec)
    projection = None
    manifest = None
    identity = None
    bootstrap_enabled = parameters.get("bootstrap_enabled", True)
    if not isinstance(bootstrap_enabled, bool):
        raise ValueError("bootstrap_enabled must be true or false")
    if adapter == "network.bootstrap" and not bootstrap_enabled:
        raise ValueError("Network bootstrap is disabled by this explicit task policy")
    if "hub_metric" in parameters and parameters["hub_metric"] not in {"abs_kME", "legacy_composite"}:
        raise ValueError("hub_metric must be abs_kME or legacy_composite")
    # Only new explicitly scoped flow tasks use the six-sample comparison
    # contract. Old sealed tasks and standalone CoSplice retain their contract.
    if "network_scope" in parameters:
        from .network_inputs import network_scope_identity, verify_network_engine_scope
        manifest = _network_project(spec)
        identity = network_scope_identity(manifest, parameters)
        for key, value in identity.items():
            if key in spec.get("scope", {}) and spec["scope"][key] != value:
                raise ValueError("Network task scope and reproducible parameters disagree")
    if inputs.get("config"):
        config = load_config(inputs["config"])
        if identity and adapter == "network.fit":
            from .network_inputs import validate_external_network_config
            validate_external_network_config(config, manifest, identity, _sha)
    elif adapter == "network.fit" and inputs.get("project"):
        from .network_inputs import prepare_project_network
        config,projection,manifest = prepare_project_network(inputs["project"],out / "prepared_inputs",parameters,_sha)
    else:
        raise ValueError("Declare the verified project manifest or effective network configuration")
    # A YAML output_dir is not producer evidence. Only a scheduler-verified
    # dependency may supply an engine_root for staged reuse.
    source_root = Path(inputs["engine_root"]).resolve() if inputs.get("engine_root") else None
    if adapter != "network.fit" and source_root is None:
        raise ValueError("A verified network.fit dependency must supply inputs.engine_root")
    if source_root is not None and identity:
        verify_network_engine_scope(source_root, identity)
    engine = out / "engine"
    if source_root is not None and (source_root == engine or engine in source_root.parents or source_root in engine.parents):
        raise ValueError("Network source and isolated stage output overlap")
    # Copy only required caches, never the old state/claim or old target stage.
    # A successful producer receipt is checked by the scheduler before this call.
    for prefix in ("02_matrix", "03_stable_as_qc", "04_network"):
        source = source_root / prefix if source_root else None
        if source is not None and source.is_dir():
            shutil.copytree(source, engine / prefix)
    config.output_dir = engine
    if "hub_metric" in parameters:
        config.interpretation.hub_metric = parameters["hub_metric"]
    if not bootstrap_enabled:
        config.network.bootstrap_iterations = 0
        config.network.consensus_bootstrap_iterations = 0
    # Execution limits must propagate to fit AND bootstrap, including plans
    # created from an existing project rather than an external config file.
    for key, target, attribute in (
        ("tom_dense_max_events", config.network, "tom_dense_max_events"),
        ("consensus_dense_max_events", config.consensus, "dense_max_events"),
        ("wgcna_timeout_seconds", config.network, "wgcna_timeout_seconds"),
    ):
        if key in parameters:
            value = parameters[key]
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{key} must be an explicit positive integer")
            setattr(target, attribute, value)
    if identity:
        config.network.minimum_network_samples = identity["minimum_network_samples"]
    if adapter == "network.fit":
        from katana_next.workflows.programs import (
            configure_production_membership, configure_production_modules,
            configure_production_uncertainty,
        )
        configure_production_modules(config, parameters)
        config.network.engine = "tom_hierarchical"
        config.network.module_detection = "wgcna_blockwise"
        configure_production_membership(config, parameters)
        configure_production_uncertainty(config, parameters)
        config.gene_prior.enabled = False
        config.multiobjective_selection.enabled = False
        config.network.max_events = int(parameters.get("max_events", 0)) or None
        if bootstrap_enabled:
            iterations = parameters.get("bootstrap_iterations", 50)
            if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 2:
                raise ValueError("Enabled bootstrap_iterations must be an integer of at least 2")
            config.network.bootstrap_iterations = iterations
            config.consensus.bootstrap_iterations = iterations
        if parameters.get("rscript"):
            config.network.wgcna_rscript = str(parameters["rscript"])
        if inputs.get("project"):
            manifest = manifest or _network_project(spec)
            from katana_next.workflows.programs import configure_interpretation_resources
            resources = out / "network_resources"; resources.mkdir()
            configure_interpretation_resources(config,manifest,resources)
            config.gene_exclusion = dict(manifest.metadata.get("gene_filter") or {})
    if adapter in {"network.fit", "network.bootstrap"}:
        # Persist the same PATH resolution used by readiness. The native engine
        # receives a concrete executable, not a bare command interpreted as a file.
        config.network.wgcna_rscript = _executable(
            parameters.get("rscript") or config.network.wgcna_rscript, "Rscript")
    config.validate()
    resolved = write_production_configuration(config, out / "resolved_network.yaml")
    pipeline = CoSplicePipeline(config)
    ran: list[str] = []
    name = str(parameters.get("network_name", "observed"))
    with persist_network_trees(), _network_r_libraries(parameters):
        if adapter == "network.fit":
            prerequisites = [("matrix", engine / "02_matrix/sample_metadata_resolved.tsv"),
                             ("stable_qc", engine / "03_stable_as_qc/stable_psi_matrix.tsv")]
            if not prerequisites[0][1].is_file():
                pipeline.run_stage("audit"); ran.append("audit")
            for stage, required in prerequisites:
                if not required.is_file():
                    if stage == "matrix" and projection:
                        from .network_inputs import matrix_display_labels
                        with matrix_display_labels(projection["condition_aliases"]):
                            pipeline.run_stage(stage)
                    else:
                        pipeline.run_stage(stage)
                    ran.append(stage)
            if identity:
                marker = engine / "network_scope.json"
                marker.write_text(json.dumps(identity, indent=2), encoding="utf-8")
                verify_network_engine_scope(engine, identity)
            outputs = pipeline.fit_network(name); ran.append(f"network_{name}")
            target = engine / "04_network" / name
            if not (target / "event_linkage.npy").is_file():
                raise ValueError("Fixed report requires numerical network linkage; configured algorithm did not produce one")
            paths = [*outputs.values(), target / "event_linkage.npy", target / "event_linkage_leaf_ids.tsv",
                     *sorted((target / "plots").glob("*.png"))]
        elif adapter == "network.bootstrap":
            outputs = pipeline.consensus_modules(); ran.append("consensus_modules")
            paths = [*outputs.values(), *sorted((engine / "04c_consensus_modules/plots").glob("*.png"))]
        elif adapter == "network.associate":
            precomputed, paths = _matched_network_associations(spec, engine, name)
            ran.append(f"associate_{name}_declared_contrasts")
        else:
            outputs = pipeline.annotate_modules(); ran.append("annotate_modules")
            paths = [*outputs.values(), *sorted((engine / "06b_module_annotation/plots").glob("*.png"))]
            precomputed, enrichment_paths = _network_enrichment(spec, engine)
            paths.extend(enrichment_paths)
            if "hub_metric" in parameters:
                hub_paths = _network_hubs(spec, engine, name)
                precomputed.update(hub_paths)
                paths.extend(map(Path, hub_paths.values()))
                ran.append(f"hub_membership_{name}_cached")
    paths.append(resolved)
    if identity:
        marker = engine / "network_scope.json"
        if not marker.is_file():
            marker.write_text(json.dumps(identity, indent=2), encoding="utf-8")
        verify_network_engine_scope(engine, identity)
        paths.append(marker)
    evidence = {"engine_root": str(engine), "resolved_config": str(resolved),
            "source_engine_root": str(source_root) if source_root else None, "executed_stages": ran,
            "fit_algorithm_unchanged": True, "original_output_modified": False,
            "dependency_compatibility": "dynamicTreeCut0.1.1 ordinal ranks retain exact values with integer index dtype",
            "analysis_id": "programs.cosplice"}
    if "bootstrap_enabled" in parameters:
        evidence.update(bootstrap_enabled=bootstrap_enabled,
                        stability_status="enabled" if bootstrap_enabled else "not_requested")
    if "hub_metric" in parameters:
        evidence["hub_metric"] = parameters["hub_metric"]
    if identity:
        evidence["network_identity"] = identity
        evidence["interpretation_status"] = "minimum_met"
    if adapter in {"network.associate", "network.annotate"}:
        evidence["precomputed_views"] = precomputed
    if projection:
        evidence["project_projection"] = projection
    return evidence, list(map(Path, paths))


def _network_project(spec):
    from katana_next.infrastructure.manifest_store import ManifestStore
    path = spec.get("inputs", {}).get("project")
    if not path:
        raise ValueError("Network interpretation requires its verified project manifest")
    return ManifestStore().load(path)[1]


def _matched_network_associations(spec, engine, name):
    import pandas as pd
    from katana_next.reports.network_views import matched_module_associations
    manifest = _network_project(spec)
    requested = list(spec.get("scope", {}).get("contrast_ids") or spec.get("parameters", {}).get("contrast_ids") or [])
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("Declare each selected network contrast exactly once; no automatic pairwise comparisons")
    known = {c.contrast_id: c for c in manifest.contrasts}
    if not set(requested) <= set(known):
        raise ValueError("A requested network contrast is absent from the project")
    contrasts = [known[cid].model_dump(mode="json") for cid in requested]
    scores = pd.read_csv(engine / f"04_network/{name}/module_program_scores.tsv",sep="\t",index_col=0)
    samples = pd.read_csv(engine / "02_matrix/sample_metadata_resolved.tsv",sep="\t")
    full, selected, audit = matched_module_associations(scores,samples,contrasts)
    rejected = [a for a in audit if a.get("status") != "completed"]
    if rejected:
        raise ValueError("A declared matched network contrast is not estimable: " + json.dumps(rejected))
    out = engine / "05_association" / name; out.mkdir(parents=True,exist_ok=True)
    paths = {"full_tests": out / "trait_full_family.tsv", "selected_tests": out / "trait_selected.tsv",
             "trait_audit": out / "trait_audit.json"}
    full.to_csv(paths["full_tests"],sep="\t",index=False)
    selected.to_csv(paths["selected_tests"],sep="\t",index=False)
    paths["trait_audit"].write_text(json.dumps(audit,indent=2),encoding="utf-8")
    design = out / "declared_contrasts.json"
    design.write_text(json.dumps(contrasts,indent=2),encoding="utf-8")
    return {k: str(v) for k,v in paths.items()}, [*paths.values(),design]


def _network_enrichment(spec, engine):
    import pandas as pd
    from katana_next.workflows.module_functions import enrich_modules
    manifest = _network_project(spec)
    libraries = dict(spec.get("inputs", {}).get("libraries") or {})
    for asset in manifest.assets:
        if asset.kind != "gene_sets_gmt":
            continue
        if asset.metadata.get("sha256") and _sha(asset.path) != asset.metadata["sha256"]:
            raise ValueError(f"Declared gene-set resource changed: {asset.asset_id}")
        key = str(asset.metadata.get("library") or asset.asset_id)
        if key in libraries and Path(libraries[key]).resolve() != asset.path.resolve():
            raise ValueError(f"Conflicting module enrichment library identity: {key}")
        libraries[key] = str(asset.path)
    if not libraries or any(not Path(value).is_file() for value in libraries.values()):
        raise ValueError("Module GO/KEGG enrichment requires explicit matched GMT libraries")
    network_name = str(spec.get("parameters", {}).get("network_name", "observed"))
    assignments = pd.read_csv(engine / "04_network" / network_name / "event_module_assignments.tsv",sep="\t")
    annotation = pd.read_csv(engine / "02_matrix/event_metadata.tsv",sep="\t")
    tested, background = enrich_modules(assignments,annotation,libraries)
    symbol = next(key for key in ("gene_symbol", "gene_name") if key in annotation)
    members = assignments[["event_key","module"]].merge(annotation[["event_key",symbol]],on="event_key",how="left",validate="one_to_one")
    members[symbol] = members[symbol].fillna("").astype(str)
    genes = members.loc[members[symbol].ne(""),["module",symbol]].drop_duplicates().rename(columns={symbol:"gene_symbol"}).sort_values(["module","gene_symbol"])
    out = engine / "06b_module_annotation"; out.mkdir(parents=True,exist_ok=True)
    paths = {"tested": out / "module_enrichment.tsv.gz", "background": out / "module_background.tsv",
             "input_genes": out / "module_input_genes.tsv"}
    tested.to_csv(paths["tested"],sep="\t",index=False,compression="gzip")
    background.to_csv(paths["background"],sep="\t",index=False)
    genes.to_csv(paths["input_genes"],sep="\t",index=False)
    resources = out / "module_enrichment_resources.json"
    resources.write_text(json.dumps({"libraries": {key:{"path":str(value),"sha256":_sha(Path(value))}
                         for key,value in libraries.items()},"background":"all unique measured network genes including grey",
                         "correction":"BH across every module and term within each explicit library"},indent=2),encoding="utf-8")
    return {k:str(v) for k,v in paths.items()}, [*paths.values(),resources]


def _network_hubs(spec, engine, name="observed"):
    """Export membership from the verified cached fit without any network refit."""
    import pandas as pd
    from katana_cosplice.module_interpretation import module_hub_tables
    metric = spec.get("parameters", {}).get("hub_metric", "abs_kME")
    source = Path(engine) / "04_network" / name
    files = {"assignments": source / "event_module_assignments.tsv",
             "psi": source / "selected_event_matrix.tsv",
             "activity": source / "module_program_scores.tsv",
             "annotation": Path(engine) / "02_matrix/event_metadata.tsv"}
    assignments = pd.read_csv(files["assignments"], sep="\t", low_memory=False)
    psi = pd.read_csv(files["psi"], sep="\t", index_col=0)
    activity = pd.read_csv(files["activity"], sep="\t", index_col=0)
    annotation = pd.read_csv(files["annotation"], sep="\t", low_memory=False)
    events, genes = module_hub_tables(assignments, psi, activity, annotation, hub_metric=metric)
    out = Path(engine) / "06e_module_interpretation"
    out.mkdir(parents=True, exist_ok=True)
    paths = {"hub_events": out / "module_hub_events.tsv", "hub_genes": out / "module_hub_genes.tsv",
             "hub_parameters": out / "hub_parameters.json"}
    events.to_csv(paths["hub_events"], sep="\t", index=False)
    genes.to_csv(paths["hub_genes"], sep="\t", index=False)
    paths["hub_parameters"].write_text(json.dumps({
        "schema": "katana.network.hub_membership/v1", "hub_metric": metric,
        "signed_kME_retained": True,
        "event_score": ("absolute Pearson correlation with cached module activity" if metric == "abs_kME" else
                        "mean available module kME/confidence/connectivity percentiles"),
        "gene_score": ("maximum absolute kME among that gene's events within each module" if metric == "abs_kME" else
                       "mean module percentiles of maximum kME, maximum event score and event count"),
        "tie_rule": ("ascending event_key for events; ascending gene_key for genes" if metric == "abs_kME" else
                     "historical shared minimum ranks"),
        "module_source": name, "network_refitted": False, "bootstrap_performed": False,
        "stability_claim": "kME measures membership, not resampling stability",
        "sources": {key: {"path": str(path), "sha256": _sha(path)} for key, path in files.items()},
    }, indent=2), encoding="utf-8")
    return {key: str(value) for key, value in paths.items()}


def _maei_tables(inputs):
    import numpy as np
    import pandas as pd
    inclusion = pd.read_csv(inputs["inclusion_table"], sep="\t", index_col=0)
    counts = pd.read_csv(inputs["sitecount_matrix"], sep="\t", index_col=0)
    design = pd.read_csv(inputs["design"], sep="\t", dtype=str)
    if inclusion.empty or counts.empty or inclusion.index.has_duplicates or counts.index.has_duplicates:
        raise ValueError("MAEI inputs need nonempty matrices with unique exact exon identities")
    if not inclusion.index.equals(counts.index):
        raise ValueError("MAEI inclusion and sitecount matrices require the same exact exon order")
    if "sample" not in design or design["sample"].duplicated().any():
        raise ValueError("Official MAEI design requires unique sample identities in column sample")
    samples = design["sample"].tolist()
    n = len(samples)
    if n < 4 or inclusion.shape[1] != 2*n:
        raise ValueError("MAEI inclusion table needs included and total TPM columns for every declared sample")
    if list(inclusion.columns[:n]) != [s + "_included" for s in samples]:
        raise ValueError("MAEI included TPM columns/order do not match the explicit design")
    if list(inclusion.columns[n:]) != [s + "_total" for s in samples]:
        raise ValueError("MAEI total TPM columns/order do not match the explicit design")
    numeric = inclusion.to_numpy(dtype=float)
    motifs = counts.to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or not np.isfinite(motifs).all() or (numeric < 0).any() or (motifs < 0).any():
        raise ValueError("Native MAEI TPM and motif count inputs must be finite and nonnegative")
    excess = numeric[:, :n] - numeric[:, n:]
    tolerance = 8 * np.finfo(float).eps * np.maximum(1.0, numeric[:, n:])
    if (excess > tolerance).any() or (numeric[:, n:] <= 0).any():
        raise ValueError("MAEI total TPM must be positive and at least included TPM")
    # Independently summed transcript subsets can differ by a few binary ULPs.
    # Normalize only those rounding excursions in a private fitting input.
    rounded = excess > 0
    inclusion.attrs["tpm_roundoff_cells"] = int(rounded.sum())
    inclusion.attrs["tpm_roundoff_max"] = float(excess[rounded].max()) if rounded.any() else 0.0
    if rounded.any():
        inclusion.iloc[:, n:] = np.maximum(numeric[:, :n], numeric[:, n:])
    if not inputs.get("quantity_provenance"):
        raise ValueError("MAEI requires explicit TPM/coordinate-bridge provenance, not raw event PSI")
    return inclusion, counts, samples


def _maei(spec: dict, out: Path) -> tuple[dict, list[Path]]:
    """Run the official standalone MAEI model on explicitly prepared TPM/counts.

    This does not substitute rMATS PSI for transcript quantities and does not
    claim to run MAPP's other eight subworkflows or KAPAC/polyadenylation.
    """
    import pandas as pd
    inputs, parameters = spec["inputs"], spec.get("parameters") or {}
    native = Path(inputs["native_root"]).resolve()
    bindings = inputs["source_sha256"]
    if set(bindings) != set(_MAPP_SCRIPTS):
        raise ValueError("Bind exactly all seven official MAEI source files by SHA256")
    for relative in _MAPP_SCRIPTS:
        if _sha(native / relative) != bindings[relative]:
            raise ValueError(f"Official MAPP source fingerprint changed: {relative}")
    inclusion, counts, samples = _maei_tables(inputs)
    fit_inclusion = Path(inputs["inclusion_table"]).resolve()
    roundoff_cells = inclusion.attrs.get("tpm_roundoff_cells", 0)
    if roundoff_cells:
        fit_inclusion = out / "INCLUSION_TABLE_roundoff_normalized.tsv"
        inclusion.to_csv(fit_inclusion, sep="\t")
    fraction = float(parameters.get("min_motif_fraction", 0.01))
    if not 0 <= fraction <= 1:
        raise ValueError("min_motif_fraction must be between zero and one")
    rscript = _executable(parameters.get("rscript"), "Rscript")
    python = _executable(parameters.get("python"), sys.executable)
    # Compile and execute byte-identical private copies, so Rcpp and Python
    # cannot leave generated files beside the user's original source checkout.
    for relative in _MAPP_SCRIPTS:
        copied = out / "native_source" / relative
        copied.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(native / relative, copied)
        if _sha(copied) != bindings[relative]:
            raise ValueError(f"Private MAEI source copy changed: {relative}")
    script_root = out / "native_source/modules/MAEI/scripts"
    activity, nuisance = out / "AUTHOR_NATIVE_ACTIVITY.tsv", out / "AUTHOR_NATIVE_NUISANCE_C.tsv"
    zscore, extended, result = (out / name for name in (
        "AUTHOR_NATIVE_ZSCORES.tsv", "AUTHOR_NATIVE_ZSCORES_EXTENDED.tsv",
        "AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv"))
    env = _r_environment(rscript, parameters)
    _run([rscript, "--vanilla", str(script_root / "optimize-MARA-model-parameters.R"),
          "--cpp_model", str(script_root / "model.cpp"), "--inclusion_table", str(fit_inclusion),
          "--sitecount_matrix", str(Path(inputs["sitecount_matrix"]).resolve()), "--out1", str(activity),
          "--out2", str(nuisance), "--min_motif_fraction", str(fraction), "--average_expressions", "FALSE"],
         out, "MAPP_MAEI_fit", env=env)
    _run([python, "-B", str(script_root / "calculate-Zscores.py"), "--activity-table", str(activity),
          "--design-table", str(Path(inputs["design"]).resolve()), "--outfile", str(zscore)], out, "MAPP_MAEI_Z", env=env)
    _run([python, "-B", str(script_root / "extend-tables.py"), "--input-table", str(zscore),
          "--output-table", str(extended)], out, "MAPP_MAEI_extend", env=env)
    # The official final display sort passes a one-column data.frame to order;
    # R >= 4.3 rejects it.  Preserve likelihood, Z, mixture and P calculations.
    compat = out / "display_order_compat.R"
    compat.write_text('order <- function(x, ...) {\n'
                      '  if (is.data.frame(x)) {\n'
                      '    if (ncol(x) != 1L) stop("Only a single display sort column is supported")\n'
                      '    x <- x[[1L]]\n  }\n  base::order(x, ...)\n}\n'
                      'source(Sys.getenv("AUTHOR_MAEI_SIGNIFICANCE_SCRIPT"), local=.GlobalEnv)\n', encoding="utf-8")
    env["AUTHOR_MAEI_SIGNIFICANCE_SCRIPT"] = str(script_root / "calculate_Zscores_statistical_significancy.R")
    _run([rscript, "--vanilla", str(compat), "--input", str(extended),
          "--sorting-strategy", "avg", "--output", str(result)], out, "MAPP_MAEI_native_significance", env=env)
    table = pd.read_csv(result, sep="\t")
    if table.empty:
        raise ValueError("Native MAPP significance output is empty")
    provenance = out / "native_source_binding.json"
    provenance.write_text(json.dumps({"source_sha256": bindings, "quantity_provenance": inputs["quantity_provenance"],
        "reference": "https://github.com/gruber-sciencelab/MAPP", "native_commit": inputs.get("native_commit"),
        "tpm_roundoff": {"cells": roundoff_cells, "max_absolute_adjustment": inclusion.attrs.get("tpm_roundoff_max", 0.0),
                         "tolerance": "8 * float64 epsilon * max(1,total TPM)", "source_unchanged": True},
        "scope": "official standalone MAEI fit/Z/mixture statistics; prepared quantities and motif counts"}, indent=2), encoding="utf-8")
    return {"native_fitting_performed": True, "quantification_rerun": False,
            "full_nine_module_MAPP": False, "KAPAC_polyadenylation": False,
            "native_statistics_unchanged": True, "sample_ids": samples,
            "events": len(inclusion), "motifs": len(counts.columns),
            "min_motif_fraction": fraction, "source_binding": str(provenance)}, [
                activity, nuisance, zscore, extended, result, compat, provenance,
                *sorted(out.glob("*.command.json")), *sorted(out.glob("INCLUSION_TABLE_roundoff_normalized.tsv")),
            ]
