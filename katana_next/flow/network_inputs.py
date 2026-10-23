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


def prepare_project_network(project_path, out, parameters, sha):
    """Project existing complete rMATS files; no alignment or AS refitting."""
    _, manifest = ManifestStore().load(project_path)
    if not manifest.samples:
        raise ValueError("Network discovery requires the declared biological-sample cohort")
    if any(normalize_sample_id(s.sample_id) != s.sample_id for s in manifest.samples):
        raise ValueError("Sample IDs must be stable identifiers; paths or BAM suffixes belong in sample assets")
    rows = []
    conditions = list(dict.fromkeys(s.condition for s in manifest.samples))
    aliases = {condition: f"group{index:04d}" for index,condition in enumerate(conditions,1)}
    groups = {c:[s.sample_id for s in manifest.samples if s.condition == c] for c in conditions}
    sample_condition = {s.sample_id:s.condition for s in manifest.samples}
    for sample in manifest.samples:
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
    mode = str(parameters.get("rmats_mode") or manifest.metadata.get("rmats_mode") or "JCEC").upper()
    if mode not in {"JC","JCEC"}:
        raise ValueError("Choose explicit JC or JCEC counts for the complete network cohort")
    contrasts = {c.contrast_id:c for c in manifest.contrasts}
    mapping,slots = [],{}
    for asset in manifest.assets:
        if asset.kind != "rmats_result" or str(asset.metadata.get("mode",mode)).upper() != mode:
            continue
        metadata = asset.metadata
        cid = metadata.get("contrast_id")
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
                raise ValueError(f"{cid}: partial-condition rMATS arms cannot define this full-cohort network matrix")
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
    mapping_path.write_text(json.dumps({"projection":mapping,"condition_aliases":aliases,
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
    return config,{"condition_aliases":aliases,"projection_mapping":str(mapping_path),
                   "transport_inputs":str(out),"statistics_refitted":False},manifest


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
