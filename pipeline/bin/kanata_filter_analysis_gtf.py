#!/usr/bin/env python3
"""Apply the approved mitochondrial/pseudogene exclusions before native fitting.

Keep the full input GTF untouched. Classification matches Katana's gene-filter
rules: declared chromosome membership and a gene biotype containing pseudogene.
Nuclear genes are never removed merely because their symbols start with MT.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

MITOCHONDRIAL_CHROMOSOMES = {"m", "mt", "chrm", "chrmt", "mitochondrion_genome"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def records(path: Path):
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                yield line, None
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {number}: expected nine fields")
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
            if not attrs.get("gene_id"):
                raise ValueError(f"GTF line {number}: missing gene_id")
            yield line, (fields, attrs)


def filter_gtf(source: Path, output: Path) -> dict:
    if output.exists():
        raise ValueError("Output exists; use a new isolated task directory")
    genes = {}
    for _, record in records(source):
        if record is None:
            continue
        fields, attrs = record
        gene_id = attrs["gene_id"]
        entry = genes.setdefault(gene_id, {"gene_id": gene_id, "chromosome": fields[0], "gene_type": "", "gene_name": ""})
        if fields[0] != entry["chromosome"]:
            raise ValueError(f"Gene appears on multiple chromosomes: {gene_id}")
        kind = attrs.get("gene_type") or attrs.get("gene_biotype", "")
        if kind:
            if entry["gene_type"] and kind != entry["gene_type"]:
                raise ValueError(f"Conflicting gene biotypes for {gene_id}")
            entry["gene_type"] = kind
        if attrs.get("gene_name"):
            entry["gene_name"] = attrs["gene_name"]
    if not genes or not any(entry["gene_type"] for entry in genes.values()):
        raise ValueError("Pseudogene exclusion requires actual gene_type/gene_biotype annotation")
    for entry in genes.values():
        reasons = []
        if entry["chromosome"].casefold() in MITOCHONDRIAL_CHROMOSOMES:
            reasons.append("mitochondrial_encoded")
        if "pseudogene" in entry["gene_type"].casefold():
            reasons.append("pseudogene")
        entry["reasons"] = ";".join(reasons)
    excluded = {gene_id for gene_id, entry in genes.items() if entry["reasons"]}
    if len(excluded) == len(genes):
        raise ValueError("All genes are excluded; no analysis reference can be built")
    output.mkdir(parents=True)
    destination = output / "annotation.analysis.gtf"
    input_features = retained_features = 0
    with destination.open("w", encoding="utf-8", newline="") as handle:
        for line, record in records(source):
            if record is None:
                handle.write(line)
            else:
                input_features += 1
                if record[1]["gene_id"] not in excluded:
                    handle.write(line)
                    retained_features += 1
    with (output / "annotation_decisions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_id", "gene_name", "gene_type", "chromosome", "excluded", "reasons"], delimiter="\t")
        writer.writeheader()
        for entry in sorted(genes.values(), key=lambda item: item["gene_id"]):
            writer.writerow({**entry, "excluded": bool(entry["reasons"])})
    receipt = {"schema": "katana.upstream_gene_exclusions/v1", "status": "completed",
        "application_stage": "before_quantification_counts_and_differential_models",
        "exclude_pseudogenes": True, "exclude_mitochondrial": True,
        "mitochondrial_chromosomes": sorted(MITOCHONDRIAL_CHROMOSOMES),
        "pseudogene_rule": "case-insensitive gene_type/gene_biotype contains pseudogene",
        "unknown_biotype_policy": "retain_and_record", "input_gtf_sha256": sha256(source),
        "analysis_gtf_sha256": sha256(destination), "input_genes": len(genes),
        "excluded_genes": len(excluded), "retained_genes": len(genes) - len(excluded),
        "unknown_biotype_genes": sum(not entry["gene_type"] for entry in genes.values()),
        "input_features": input_features, "retained_features": retained_features}
    (output / "filtering.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(filter_gtf(args.gtf, args.outdir)))


if __name__ == "__main__":
    main()
