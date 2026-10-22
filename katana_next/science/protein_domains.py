"""Compare sequence-scanned isoforms, not transferred reference coordinates."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import uuid

import pandas as pd

from katana_next.science.isoform_consequences import normalized_id
from katana_next.science.motif_significance import _wsl_path


def resource_paths(context):
    databases = [a for a in context.manifest.assets if a.kind == "pfam_hmm"]
    runtimes = [a for a in context.manifest.assets if a.kind == "hmmer_runtime"]
    if not databases and not runtimes:
        return None
    if len(databases) != 1 or len(runtimes) != 1:
        raise ValueError("Exactly one Pfam database and one HMMER runtime are required")
    return databases[0],runtimes[0]


def scan_domains(features, database, runtime, directory, job, *, cache_root=None, previous_directories=()):
    directory.mkdir(parents=True,exist_ok=True)
    with database.path.open("rb") as handle:
        digest = hashlib.file_digest(handle,"sha256").hexdigest()
    if digest != database.metadata.get("sha256") or not all(database.metadata.get(k) for k in ("version","source_url")):
        raise ValueError("Pfam requires matching SHA256, database version and source URL")
    valid = features[features.protein_status.eq("validated_annotated_ORF")]
    records=valid[["transcript_id","protein_sequence"]].to_dict("records")
    config=json.loads(runtime.path.read_text())
    if not config.get("version"):
        raise ValueError("Declare the exact PyHMMER runtime version")
    proteins, result_path = directory/"proteins.json", directory/"domains.json"
    worker=Path(__file__).resolve().parents[2]/"scripts"/"pfam_worker.py"
    key=hashlib.sha256(json.dumps({"proteins":records,"database":digest,"runtime":config,
        "worker":hashlib.sha256(worker.read_bytes()).hexdigest()},sort_keys=True).encode()).hexdigest()
    marker=directory/"checkpoint.json"
    expected=dict(zip(valid.transcript_id,valid.protein_sha256))
    def validate(result):
        if result["engine_version"] != config["version"]:
            raise ValueError("Installed PyHMMER version differs from declared runtime")
        if {r["transcript_id"]:r["protein_sha256"] for r in result["proteins"]} != expected or result["profiles_scanned"] < 1:
            raise ValueError("Pfam scan omitted proteins or scanned no models")
    def publish(source, target):
        target.mkdir(parents=True,exist_ok=True)
        for name in ("domains.json", "checkpoint.json"):
            temporary=target/(name+"."+uuid.uuid4().hex+".tmp")
            shutil.copyfile(source/name,temporary)
            temporary.replace(target/name)
    shared=cache_root/key if cache_root else None
    candidates=[directory]+([shared] if shared else [])+list(previous_directories)
    for candidate in dict.fromkeys(candidates):
        job.check_cancelled()
        try:
            saved=json.loads((candidate/"checkpoint.json").read_text())
            if saved["key"] != key:
                continue
            raw=(candidate/"domains.json").read_bytes()
            if saved["result_sha256"] != hashlib.sha256(raw).hexdigest():
                continue
            result=json.loads(raw)
            validate(result)
        except (OSError,ValueError,KeyError,TypeError):
            continue
        if candidate != directory:
            publish(candidate,directory)
            (directory/"reuse.json").write_text(json.dumps({"source":str(candidate),"scan_fingerprint":key,
                "result_sha256":saved["result_sha256"],"scope":"Identical protein sequences, database, runtime and worker; checked complete protein coverage"},indent=2))
        if shared and shared != candidate:
            publish(directory,shared)
        job.report(.4,"Pfam cache hit: sequence, database, runtime and complete scan verified")
        return result,key
    proteins.write_text(json.dumps(records),encoding="utf-8")
    if not records:
        return {"proteins":[],"domains":[],"profiles_scanned":0,"status":"no_validated_proteins"},key
    paths=[worker,database.path,proteins,result_path]
    command=(["wsl","-d",config["wsl_distribution"],"--",config["python"],"-B",*map(_wsl_path,paths)]
             if config.get("wsl_distribution") else [config["python"],"-B",*map(str,paths)])
    (directory/"command.json").write_text(json.dumps({"command":command,"database":database.metadata,"runtime":config},indent=2))
    with (directory/"scan.log").open("w",encoding="utf-8") as log:
        process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
        try:
            while True:
                job.check_cancelled()
                try:
                    code=process.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    continue
        except BaseException:
            process.terminate(); process.wait(timeout=15)
            raise
    if code:
        raise RuntimeError(f"Pfam scan failed; inspect {directory/'scan.log'}")
    result=json.loads(result_path.read_text())
    validate(result)
    marker.write_text(json.dumps({"key":key,"result_sha256":hashlib.sha256(result_path.read_bytes()).hexdigest()}))
    if shared:
        publish(directory,shared)
    return result,key


def domain_consequences(switches, scan):
    scanned={r["transcript_id"] for r in scan["proteins"]}
    counts={tx:Counter() for tx in scanned}
    for hit in scan["domains"]:
        if hit["transcript_id"] not in scanned:
            raise ValueError("Domain hit has no scanned protein")
        counts[hit["transcript_id"]][hit["domain_accession"]] += 1
    result=switches.copy(); rows=[]
    for ix,row in result.iterrows():
        a,b=map(normalized_id,(row.numerator_dominant_transcript,row.denominator_dominant_transcript))
        if a not in scanned or b not in scanned or not row.annotation_gene_match:
            result.loc[ix,"domain_consequence"]="not_evaluable_unscanned_or_gene_mismatch"
            continue
        gains,losses = counts[a]-counts[b],counts[b]-counts[a]
        result.loc[ix,"domain_consequence"]=("predicted_domain_gain_and_loss" if gains and losses else "predicted_domain_gain" if gains else "predicted_domain_loss" if losses else "same_detected_domain_counts")
        for domain in sorted(set(counts[a]) | set(counts[b])):
            rows.append({"source_row_id":row.row_id,"comparison_id":row.comparison_id,"gene_id":row.gene_id,
                "numerator_transcript":a,"denominator_transcript":b,"domain_accession":domain,
                "numerator_domain_count":counts[a][domain],"denominator_domain_count":counts[b][domain],
                "delta_domain_count":counts[a][domain]-counts[b][domain],
                "scope":"GA-threshold sequence prediction; absence means not detected, not proof of biological domain loss"})
    return result,pd.DataFrame(rows,columns=["source_row_id","comparison_id","gene_id","numerator_transcript","denominator_transcript","domain_accession","numerator_domain_count","denominator_domain_count","delta_domain_count","scope"])
