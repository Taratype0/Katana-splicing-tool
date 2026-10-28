"""Project-to-CoSplice transport, preserving native vector/sample identity."""
import json
from pathlib import Path
import shutil
from contextlib import contextmanager
from dataclasses import replace
from unittest.mock import patch

import pandas as pd

from katana_cosplice.config import CoSpliceConfig
from katana_cosplice.utils import normalize_sample_id
from katana_next.infrastructure.manifest_store import ManifestStore


def resolve_network_rmats_mode(manifest, parameters):
    """Bind JC/JCEC from an explicit request or one declared scoped source mode.

    Never infer a count universe from a filename or default a JC-only project
    to JCEC. Ambiguous dual-mode input requires an explicit user/project choice.
    """
    scope = parameters.get("network_scope", "cohort")
    selected = set(parameters.get("contrast_ids") or [])
    assets = [a for a in manifest.assets if a.kind == "rmats_result"
              and (scope != "contrast" or a.metadata.get("contrast_id") in selected)]
    if not assets:
        raise ValueError("No declared rMATS measurement assets exist for this network scope")
    modes = set()
    for asset in assets:
        mode = str(asset.metadata.get("mode") or "").upper()
        if mode not in {"JC", "JCEC"}:
            raise ValueError(f"rMATS asset {asset.asset_id} must declare an explicit JC/JCEC mode")
        modes.add(mode)
    request = parameters.get("rmats_mode")
    if request is None:
        request = manifest.metadata.get("rmats_mode")
    if request is None:
        if len(modes) != 1:
            raise ValueError("Network input contains both JC and JCEC; explicitly choose network.parameters.rmats_mode or project rmats_mode")
        return next(iter(modes))
    mode = str(request).upper()
    if mode not in {"JC", "JCEC"}:
        raise ValueError("Network rmats_mode must explicitly be JC or JCEC")
    if mode not in modes:
        raise ValueError(f"Requested network mode {mode} has no declared scoped rMATS assets; available: {sorted(modes)}")
    return mode


def network_scope_identity(manifest, parameters, *, preserve_recorded_status=True):
    """Bind a discovery to its exact biological samples, independently of labels.

    Missing scope is the historical cohort contract. New flow plans explicitly
    choose contrast; persisted plans and the standalone engine keep their rules.
    """
    scope = parameters.get("network_scope", "cohort")
    if scope not in {"contrast", "cohort"}:
        raise ValueError("network.scope must be contrast or cohort")
    known = {c.contrast_id: c for c in manifest.contrasts}
    contrasts = list(parameters.get("contrast_ids") or [])
    if len(set(contrasts)) != len(contrasts) or set(contrasts) - set(known):
        raise ValueError("Network scope needs unique declared contrast IDs")
    if scope == "contrast":
        if len(contrasts) != 1:
            raise ValueError("Each contrast network must declare exactly one comparison")
        contrast = known[contrasts[0]]
        wanted = set(contrast.numerator.sample_ids) | set(contrast.denominator.sample_ids)
        samples = [s for s in manifest.samples if s.sample_id in wanted]
    else:
        samples = list(manifest.samples)
    sample_ids = [s.sample_id for s in samples]
    biological_ids = [s.biological_sample_id for s in samples]
    repeated = len(set(biological_ids)) != len(biological_ids)
    declared = manifest.metadata.get("comparison_designs", {}).get(contrasts[0], {}) if scope == "contrast" else {}
    paired = declared.get("paired") is True
    if len(set(sample_ids)) != len(sample_ids) or (repeated and not paired):
        raise ValueError("Network discovery requires distinct independent biological samples; technical records cannot inflate n")
    if paired:
        arms = [[s.biological_sample_id for s in samples if s.sample_id in arm.sample_ids]
                for arm in (contrast.numerator, contrast.denominator)]
        if not repeated or len(set(arms[0])) != len(arms[0]) or len(set(arms[1])) != len(arms[1]) or set(arms[0]) != set(arms[1]):
            raise ValueError("Paired contrast network requires exactly one sample per subject in each arm")
        if len(set(biological_ids)) < 3:
            raise ValueError("Paired module association requires at least three complete subject pairs")
    minimum = parameters.get("minimum_network_samples", 6 if scope == "contrast" else 15)
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 3:
        raise ValueError("network.parameters.minimum_network_samples must be an integer of at least 3")
    if len(samples) < minimum:
        raise ValueError(
            f"{scope} network has {len(samples)} sample observations; "
            f"the configured gate requires at least {minimum}. Set network.parameters.minimum_network_samples "
            "to your intended minimum (at least 3), or supply more independent samples."
        )
    legacy_status = "small_sample_exploratory" if len(samples) < 15 else "standard_sample_count"
    recorded_status = parameters.get("sample_size_status")
    # An older sealed task must retain its exact identity for verified reuse.
    # New plans use a neutral eligibility value, with actual n recorded beside it.
    status = legacy_status if preserve_recorded_status and recorded_status == legacy_status else "minimum_met"
    identity = {"network_scope": scope, "contrast_ids": contrasts,
                "sample_ids": sample_ids, "biological_sample_ids": biological_ids,
                "n_samples": len(samples), "minimum_network_samples": minimum,
                "sample_size_status": status}
    if paired:
        identity.update(paired=True, n_independent_units=len(set(biological_ids)),
                        correlation_observations=len(samples),
                        association_model="subject-blocked frozen module scores",
                        discovery_inference="observed correlation and descriptive kME; observations include repeated subjects")
    for key, value in identity.items():
        if key == "sample_size_status" and not preserve_recorded_status and recorded_status == legacy_status:
            continue
        if key in parameters and parameters[key] != value:
            raise ValueError(f"Network {key} disagrees with the declared project scope")
    return identity


def validate_external_network_config(config, manifest, identity, sha):
    """Reject an external config that silently widens or changes a contrast.

    The same native rMATS vectors and their exact side order must be bound to
    the declared samples; an unrelated cohort YAML cannot masquerade as a fit.
    """
    from katana_cosplice.audit import load_sample_metadata, resolve_rmats_root
    from katana_cosplice.matrix import _candidate_sources, _group_samples
    expected_mode = resolve_network_rmats_mode(manifest, {**identity, "rmats_mode": config.rmats_mode})
    if config.rmats_mode != expected_mode:
        raise ValueError("External network configuration mode differs from the declared input mode")
    samples = load_sample_metadata(config.sample_metadata)
    if set(samples.sample_id) != set(identity["sample_ids"]):
        raise ValueError("External network configuration samples do not exactly match the task scope")
    if "biological_sample_id" not in samples:
        raise ValueError("Scoped external network configuration requires biological_sample_id metadata")
    expected_units = dict(zip(identity["sample_ids"], identity["biological_sample_ids"]))
    if samples.set_index("sample_id").biological_sample_id.astype(str).to_dict() != expected_units:
        raise ValueError("External network biological units differ from the declared project")
    selected = load_sample_metadata(config.sample_metadata, config.selected_conditions)
    if set(selected.sample_id) != set(identity["sample_ids"]):
        raise ValueError("External network selected_conditions changes the declared sample scope")
    if identity["network_scope"] != "contrast":
        return
    if config.matrix_input.input_mode != "rmats_events":
        raise ValueError("Contrast flow configuration currently requires declared rMATS sources")
    cid = identity["contrast_ids"][0]
    assets = [a for a in manifest.assets if a.kind == "rmats_result"
              and a.metadata.get("contrast_id") == cid
              and str(a.metadata.get("mode", config.rmats_mode)).upper() == config.rmats_mode]
    by_digest = {}
    for asset in assets:
        digest = sha(asset.path)
        expected = asset.metadata.get("upstream_sha256") or asset.metadata.get("sha256")
        if expected and expected != digest:
            raise ValueError("Declared contrast rMATS source fingerprint changed")
        by_digest.setdefault(digest, []).append(asset)
    root = resolve_rmats_root(config.rmats_root)
    groups = sorted(samples.condition.unique().tolist(), key=len, reverse=True)
    sources, counts = _candidate_sources(root, config.rmats_mode, groups)
    if len(sources) != len(groups) * 5 or any(n != 1 for n in counts.values()):
        raise ValueError("External contrast configuration needs one unambiguous native source per group and event class")
    for (group, event), (path, side) in sources.items():
        order = _group_samples(root / "bamlist" / f"{group}_bamlist.txt")
        if set(order) != set(samples.loc[samples.condition.eq(group), "sample_id"]):
            raise ValueError("External rMATS bamlist differs from the scoped sample metadata")
        matches = by_digest.get(sha(path), [])
        valid = False
        for asset in matches:
            arm_order = ("denominator", "numerator") if asset.metadata.get("reverse_direction", False) else ("numerator", "denominator")
            arm = arm_order[side - 1]
            if (str(asset.metadata.get("event_type") or asset.path.name.split(".")[0]) == event
                    and order == list(asset.metadata.get(arm + "_sample_ids") or [])):
                valid = True
        if not valid:
            raise ValueError("External network rMATS source or native sample order is outside this exact contrast")


def verify_network_engine_scope(engine, identity):
    """Followups cannot consume another contrast's fit, even with equal n."""
    marker = Path(engine) / "network_scope.json"
    if not marker.is_file() or json.loads(marker.read_text(encoding="utf-8")) != identity:
        raise ValueError("Network fit dependency scope does not match this exact task")
    samples = pd.read_csv(Path(engine) / "02_matrix/sample_metadata_resolved.tsv", sep="\t", dtype=str)
    if (not samples.sample_id.is_unique or set(samples.sample_id) != set(identity["sample_ids"])):
        raise ValueError("Network fit matrix includes missing or foreign samples")
    if "biological_sample_id" not in samples:
        raise ValueError("Scoped network fit matrix requires biological_sample_id metadata")
    units = samples.set_index("sample_id").biological_sample_id.astype(str).to_dict()
    if units != dict(zip(identity["sample_ids"], identity["biological_sample_ids"])):
        raise ValueError("Network fit matrix biological units differ from the task scope")


def prepare_project_network(project_path, out, parameters, sha):
    """Project existing complete rMATS files; no alignment or AS refitting."""
    _, manifest = ManifestStore().load(project_path)
    identity = network_scope_identity(manifest, parameters)
    mode = resolve_network_rmats_mode(manifest, parameters)
    selected_samples = [s for s in manifest.samples if s.sample_id in set(identity["sample_ids"])]
    if not manifest.samples:
        raise ValueError("Network discovery requires the declared biological-sample cohort")
    if any(normalize_sample_id(s.sample_id) != s.sample_id for s in manifest.samples):
        raise ValueError("Sample IDs must be stable identifiers; paths or BAM suffixes belong in sample assets")
    rows = []
    conditions = list(dict.fromkeys(s.condition for s in selected_samples))
    aliases = {condition: f"group{index:04d}" for index,condition in enumerate(conditions,1)}
    groups = {c:[s.sample_id for s in selected_samples if s.condition == c] for c in conditions}
    sample_condition = {s.sample_id:s.condition for s in manifest.samples}
    for sample in selected_samples:
        rows.append({**sample.attributes,"sample_id":sample.sample_id,
                     "biological_sample_id":sample.biological_sample_id,
                     "condition":aliases[sample.condition],"condition_label":sample.condition,
                     **({"batch":sample.batch} if sample.batch else {})})
    out.mkdir(parents=True,exist_ok=True)
    samples_path = out / "sample_metadata_engine.tsv"
    pd.DataFrame(rows).to_csv(samples_path,sep="\t",index=False)
    root = out / "rmats_projection"
    (root / "bamlist").mkdir(parents=True)
    for group,ids in groups.items():
        (root / "bamlist" / (aliases[group]+"_bamlist.txt")).write_text(",".join(ids)+"\n",encoding="utf-8")
    contrasts = {c.contrast_id:c for c in manifest.contrasts}
    mapping,slots = [],{}
    for asset in manifest.assets:
        if asset.kind != "rmats_result" or str(asset.metadata.get("mode",mode)).upper() != mode:
            continue
        metadata = asset.metadata
        cid = metadata.get("contrast_id")
        if identity["network_scope"] == "contrast" and cid not in identity["contrast_ids"]:
            continue
        if cid not in contrasts:
            raise ValueError(f"rMATS resource has an undeclared comparison: {asset.asset_id}")
        contrast = contrasts[cid]
        arms = {}
        for arm in ("numerator","denominator"):
            declared = list(getattr(contrast,arm).sample_ids)
            observed = list(metadata.get(arm+"_sample_ids") or [])
            if len(set(observed)) != len(observed) or set(observed) != set(declared):
                raise ValueError(f"{cid}: native {arm} sample vector order is missing or incompatible")
            labels = {sample_condition[s] for s in observed}
            if len(labels) != 1:
                raise ValueError(f"{cid}: one rMATS arm must represent one declared condition")
            condition = next(iter(labels))
            if set(groups[condition]) != set(observed):
                raise ValueError(f"{cid}: rMATS arms differ from the scoped network samples")
            arms[arm] = (condition,observed)
        order = ("denominator","numerator") if metadata.get("reverse_direction",False) else ("numerator","denominator")
        first,second = (arms[arm] for arm in order)
        event = str(metadata.get("event_type") or asset.path.name.split(".")[0])
        if event not in {"SE","A3SS","A5SS","MXE","RI"}:
            raise ValueError(f"Unknown rMATS event class: {event}")
        destination = root / (aliases[first[0]]+"-"+aliases[second[0]]) / "rmats_post" / f"{event}.MATS.{mode}.txt"
        source_sha = sha(asset.path)
        expected = metadata.get("upstream_sha256") or metadata.get("sha256")
        if expected and source_sha != expected:
            raise ValueError(f"rMATS source fingerprint changed: {asset.asset_id}")
        if destination in slots:
            if slots[destination] != source_sha:
                raise ValueError("Repeated comparison arms have conflicting native rMATS sources; declare an unambiguous measurement source")
            continue
        destination.parent.mkdir(parents=True,exist_ok=True)
        reorders = [([original.index(sid) for sid in groups[condition]]) for condition,original in (first,second)]
        changed = any(indices != list(range(len(indices))) for indices in reorders)
        if not changed:
            shutil.copyfile(asset.path,destination)
            if sha(destination) != source_sha:
                raise ValueError("Native rMATS projection bytes changed during copying")
        else:
            frame = pd.read_csv(asset.path,sep="\t")
            for side,indices in enumerate(reorders,1):
                for column in (f"IncLevel{side}",f"IJC_SAMPLE_{side}",f"SJC_SAMPLE_{side}"):
                    if column not in frame:
                        raise ValueError(f"rMATS input lacks required sample vector {column}")
                    def reorder(value):
                        tokens = str(value).split(",")
                        if len(tokens) != len(indices):
                            raise ValueError(f"{cid}/{event}: vector length differs from declared sample identities")
                        return ",".join(tokens[i] for i in indices)
                    frame[column] = frame[column].map(reorder)
            frame.to_csv(destination,sep="\t",index=False)
        slots[destination] = source_sha
        mapping.append({"contrast_id":cid,"event_type":event,"mode":mode,"source":str(asset.path),
                        "source_sha256":source_sha,"projected":str(destination),"projected_sha256":sha(destination),
                        "native_arm1_sample_ids":first[1],"native_arm2_sample_ids":second[1],
                        "matrix_arm1_sample_ids":groups[first[0]],"matrix_arm2_sample_ids":groups[second[0]],
                        "vector_order_transport_changed":changed,"statistics_refitted":False})
    if not mapping:
        raise ValueError("No compatible existing rMATS measurements are declared for the network")
    available = {row["event_type"] for row in mapping}
    if available != {"SE","A3SS","A5SS","MXE","RI"}:
        raise ValueError("The network measurement source must retain all five rMATS event-class schemas")
    mapping_path = out / "rmats_projection.mapping.json"
    mapping_path.write_text(json.dumps({"projection":mapping,"condition_aliases":aliases,"network_identity":identity,
                           "sample_identity_preserved":True,"source_results_modified":False},indent=2),encoding="utf-8")
    config = CoSpliceConfig(project_name=manifest.name,rmats_root=root,sample_metadata=samples_path,
        output_dir=out / "engine",gtf=manifest.annotation.gtf,genome_fasta=manifest.annotation.genome_fasta,
        transcript_fasta=manifest.annotation.transcript_fasta,genome_build=manifest.annotation.genome_build,
        annotation_version=str(manifest.metadata.get("annotation_version", "")),rmats_mode=mode,
        selected_conditions=list(aliases.values()),gene_exclusion=dict(manifest.metadata.get("gene_filter") or {}))
    rscript = parameters.get("rscript") or manifest.metadata.get("rscript_path")
    if rscript:
        config.network.wgcna_rscript = str(rscript)
    config.network.min_module_size = 20
    config.network.wgcna_merge_cut_height = .25
    config.network.minimum_network_samples = identity["minimum_network_samples"]
    return config,{"condition_aliases":aliases,"projection_mapping":str(mapping_path),
                   "transport_inputs":str(out),"statistics_refitted":False,"network_identity":identity},manifest


def restore_condition_labels(engine, aliases):
    """Only transport aliases change; biological groups/order remain identical."""
    inverse = {token:label for label,token in aliases.items()}
    target = engine / "02_matrix/sample_metadata_resolved.tsv"
    frame = pd.read_csv(target,sep="\t",dtype=str)
    if not set(frame.condition) <= set(inverse):
        raise ValueError("Projected network metadata contains an unknown condition alias")
    frame["condition"] = frame.condition.map(inverse)
    frame.to_csv(target,sep="\t",index=False)


@contextmanager
def matrix_display_labels(aliases):
    """Resolve transport labels before the normal writer commits its outputs.

    No count/PSI values, sample IDs, row order or biological grouping change.
    This avoids altering a matrix file after its official artifact hash is saved.
    """
    import katana_cosplice.pipeline as pipeline
    inverse = {token:label for label,token in aliases.items()}
    original = pipeline.build_sample_matrices
    def build(config):
        result = original(config)
        samples = result.samples.copy()
        if not set(samples.condition) <= set(inverse):
            raise ValueError("Projected network metadata contains an unknown condition alias")
        samples["condition"] = samples.condition.map(inverse)
        return replace(result,samples=samples)
    with patch.object(pipeline,"build_sample_matrices",build):
        yield
