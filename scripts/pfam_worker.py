"""Versioned Pfam/HMMER scanning of every supplied validated protein."""
import argparse
import hashlib
import json
from pathlib import Path

import pyhmmer


def scan(database, proteins, output, cpus=2):
    sequences = json.loads(Path(proteins).read_text())
    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
        pyhmmer.easel.TextSequence(name=row["transcript_id"].encode(),sequence=row["protein_sequence"]).digitize(alphabet)
        for row in sequences])
    domains, profiles = [], 0
    with pyhmmer.plan7.HMMFile(database) as models:
        for hits in pyhmmer.hmmer.hmmsearch(models, block, cpus=cpus, bit_cutoffs="gathering"):
            profiles += 1
            for hit in hits.included:
                for domain in hit.domains.included:
                    alignment = domain.alignment
                    domains.append({"transcript_id":hit.name.decode(),"domain_accession":hits.query.accession.decode(),
                        "domain_name":hits.query.name.decode(),"start_1based":alignment.target_from,
                        "end_inclusive":alignment.target_to,"hmm_start":alignment.hmm_from,"hmm_end":alignment.hmm_to,
                        "hmm_length":hits.query.M,"domain_score":domain.score,"domain_i_evalue":domain.i_evalue})
    result = {"engine":"PyHMMER/hmmsearch", "engine_version":pyhmmer.__version__, "threshold":"Pfam model-specific GA sequence and domain thresholds",
        "profiles_scanned":profiles,"proteins":[{"transcript_id":r["transcript_id"],"protein_sha256":hashlib.sha256(r["protein_sequence"].encode()).hexdigest()} for r in sequences],
        "domains":domains,"limitations":"Sequence-model predictions; no clan overlap arbitration, no experimental proof of domain function"}
    temp = Path(output).with_suffix(".partial")
    temp.write_text(json.dumps(result),encoding="utf-8")
    temp.replace(output)


if __name__ == "__main__":
    p=argparse.ArgumentParser(__doc__)
    p.add_argument("database"); p.add_argument("proteins"); p.add_argument("output")
    p.add_argument("--cpus",type=int,default=2)
    args=p.parse_args(); scan(args.database,args.proteins,args.output,args.cpus)
