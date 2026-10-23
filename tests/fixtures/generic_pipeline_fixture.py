"""Generate isolated generic-pipeline installation fixtures, never study data.

Run this module in the target environment so all saved resource paths belong to
that environment. The caller adds generated BAM/QC inputs before sealing a plan.
No external files, public datasets, previously completed analyses or network
connections are used. Fixture statistics are intentionally generated inputs;
new MotifTrace statistics must be computed by the real registered workflow.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from katana_next.full_fixture import create_full_fixture


PATTERNS = ("AAA", "CCC", "GGG", "TTT", "ACG", "CGT", "GTA", "TAC", "AAC", "CCA", "GGT", "TTA")


def create_generic_fixture(directory):
    root = Path(directory).resolve()
    if root.exists():
        raise FileExistsError("A generated fixture must use a new directory")
    project_path = create_full_fixture(root)
    project = json.loads(project_path.read_text(encoding="utf-8"))
    project["project_id"] = "generic-pipeline-generated-acceptance"
    project["name"] = "Generated software acceptance fixture (not biological results)"
    project["samples"] = [s for s in project["samples"] if int(s["sample_id"].rsplit("-", 1)[1]) <= 3]
    sample_ids = {s["sample_id"] for s in project["samples"]}
    project["contrasts"] = project["contrasts"][:2]
    contrasts = {c["contrast_id"]: c for c in project["contrasts"]}
    for sample in project["samples"]:
        sample["biological_sample_id"] = sample["sample_id"]
        sample["batch"] = "synthetic-batch"
    for contrast in contrasts.values():
        for arm in ("numerator", "denominator"):
            contrast[arm]["sample_ids"] = [s for s in contrast[arm]["sample_ids"] if s in sample_ids]
    project["assets"] = [a for a in project["assets"]
                         if not a["kind"].startswith("single_cell")
                         and (not a.get("metadata", {}).get("contrast_id")
                              or a["metadata"]["contrast_id"] in contrasts)]
    vectors = ("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel1", "IncLevel2")
    templates = []
    for asset in project["assets"]:
        if asset["kind"] != "rmats_result":
            continue
        path = Path(asset["path"])
        frame = pd.read_csv(path, sep="\t")
        for column in vectors:
            frame[column] = frame[column].map(lambda x: ",".join(str(x).split(",")[:3]))
        for name in ("numerator_sample_ids", "denominator_sample_ids"):
            asset["metadata"][name] = [s for s in asset["metadata"][name] if s in sample_ids]
        kind = asset["metadata"]["event_type"]
        # Ensure all generated A3SS/A5SS and RI coordinate contracts are valid
        # for the real local-RNA extractor, independent of network use.
        if kind in {"A3SS", "A5SS"}:
            for idx, row in frame.iterrows():
                suffix = (kind == "A5SS") == (row["strand"] == "+")
                start = int(row["longExonStart_0base"])
                frame.loc[idx, ["shortES", "shortEE"]] = [start, start + 95] if suffix else [start + 45, start + 140]
                frame.loc[idx, ["flankingES", "flankingEE"]] = [start + 300, start + 400] if suffix else [start - 300, start - 200]
        if kind == "RI":
            frame["riExonStart_0base"] = frame["upstreamES"]
            frame["riExonEnd"] = frame["downstreamEE"]
        if kind == "SE":
            extra = []
            for i in range(48):
                foreground = i < 24
                start = 60_000 + i * 700
                gene = 101 + i
                row = frame.iloc[0].to_dict()
                row.update(ID=1000 + i, GeneID=f"ENSMUSG{gene:011d}.1", geneSymbol=f"GENE{gene}",
                           chr="chr1", strand="+", exonStart_0base=start, exonEnd=start + 36,
                           upstreamES=start - 300, upstreamEE=start - 200,
                           downstreamES=start + 300, downstreamEE=start + 400)
                # Same denominators allow motif-blind matching between unique
                # foreground/background genes. Both comparisons use six samples.
                if foreground:
                    base = .30 + (i % 3 - 1) * .005
                    shift = .30 if asset["metadata"]["contrast_id"] == next(iter(contrasts)) else .20
                    one, two = [base + shift + d for d in (-.005, 0, .005)], [base + d for d in (-.005, 0, .005)]
                else:
                    base = .30 + (i % 3 - 1) * .005
                    one, two = [base + d for d in (-.005, 0, .005)], [base + d for d in (-.005, 0, .005)]
                row.update(IncLevel1=",".join(map(str, one)), IncLevel2=",".join(map(str, two)),
                           IncLevelDifference=float(np.mean(one) - np.mean(two)),
                           FDR=.001 if foreground else .8)
                for suffix, values in (("1", one), ("2", two)):
                    row["IJC_SAMPLE_" + suffix] = ",".join(str(round(x * 100)) for x in values)
                    row["SJC_SAMPLE_" + suffix] = ",".join(str(100 - round(x * 100)) for x in values)
                extra.append(row)
            frame = pd.concat([frame, pd.DataFrame(extra)], ignore_index=True)
            templates = extra
        frame.to_csv(path, sep="\t", index=False)

    inputs = root / "inputs"
    for condition in {s["condition"] for s in project["samples"]}:
        ids = [s["sample_id"] for s in project["samples"] if s["condition"] == condition]
        (inputs / "rmats/bamlist" / (condition + "_bamlist.txt")).write_text(",".join(ids), encoding="utf-8")
    pd.DataFrame([{k: s[k] for k in ("sample_id", "condition", "biological_sample_id")} for s in project["samples"]]).to_csv(inputs / "cosplice_samples.csv", index=False)

    genome = list("ACGT" * 26_000)
    for i, row in enumerate(templates):
        word = "AAACCCGGGTTT" * 3 if i < 24 else "ACGT" * 9
        genome[row["exonStart_0base"]:row["exonEnd"]] = word
    fasta = inputs / "fixture.genome.fa"
    sequence = "".join(genome)
    chunks, fai, offset = [], [], 0
    for chrom in ("chr1", "chr2", "chr3", "chr4"):
        header = ">" + chrom + "\n"
        offset += len(header)
        fai.append(f"{chrom}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n")
        chunks.append(header + sequence + "\n")
        offset += len(sequence) + 1
    fasta.write_text("".join(chunks), encoding="ascii", newline="\n")
    Path(str(fasta) + ".fai").write_text("".join(fai), encoding="ascii", newline="\n")

    # Append coherent reference transcripts for every added synthetic gene.
    gtf = inputs / "fixture.annotation.gtf"
    with gtf.open("a", encoding="utf-8", newline="\n") as stream:
        for row in templates:
            attrs = f'gene_id "{row["GeneID"]}"; transcript_id "TX{row["ID"]}"; gene_name "{row["geneSymbol"]}"; gene_type "protein_coding";'
            for left, right in ((row["upstreamES"], row["upstreamEE"]), (row["exonStart_0base"], row["exonEnd"]), (row["downstreamES"], row["downstreamEE"])):
                stream.write(f"chr1\tgenerated_fixture\texon\t{left + 1}\t{right}\t.\t+\t.\t{attrs}\n")
    with (inputs / "fixture.pathways.gmt").open("a", encoding="utf-8") as stream:
        stream.write("Generated_SE_foreground\tsynthetic fixture only\t" + "\t".join(f"GENE{i}" for i in range(101, 125)) + "\n")
        stream.write("Generated_SE_background\tsynthetic fixture only\t" + "\t".join(f"GENE{i}" for i in range(125, 149)) + "\n")

    def resource(identity, kind, path, **metadata):
        return {"asset_id": identity, "kind": kind, "path": str(path), "metadata": {
            "source_url": "urn:katana:generated-software-fixture", "version": "generated-v1",
            "species": "Mus_musculus", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "synthetic": True, **metadata}}

    motifs = inputs / "generated.consensus.tsv"
    pd.DataFrame([dict(motif_id=f"M{i:02}", RBP=f"FixtureRBP{i:02}", species="Mus_musculus",
                       evidence="user_defined", source="generated fixture only", motif=pattern)
                  for i, pattern in enumerate(PATTERNS)]).to_csv(motifs, sep="\t", index=False)
    pwm = inputs / "generated.pwm.tsv"
    rows = []
    for i, pattern in enumerate(PATTERNS):
        for position, base in enumerate(pattern, 1):
            rows.append({"motif_id": f"M{i:02}", "RBP": f"FixtureRBP{i:02}", "position": position,
                         "Organism": "Mus_musculus", "evidence_class": "user_defined",
                         "source": "generated fixture only", **{b: .97 if b == base else .01 for b in "ACGT"}})
    pd.DataFrame(rows).to_csv(pwm, sep="\t", index=False)
    runtime = inputs / "generated.fimo-runtime.json"
    runtime.write_text(json.dumps({"python": sys.executable}), encoding="utf-8")
    peaks = inputs / "generated.clip.bed"
    peaks.write_text("".join(f'chr1\t{row["exonStart_0base"] - 5}\t{row["exonEnd"] + 5}\tfixture_peak_{i}\t0\t+\n'
                             for i, row in enumerate(templates[:12])), encoding="utf-8")
    project["assets"].extend([
        resource("generated-pwm", "rbp_motif_database", pwm),
        resource("generated-fimo", "fimo_runtime", runtime),
        resource("generated-clip", "rbp_binding_peaks", peaks, genome_build="fixture-v1",
                 coordinate_system="0-based-half-open", assay="eCLIP", biosample="generated nonbiological fixture",
                 target="FixtureRBP00", accession="SYNTHETIC_NOT_EXPERIMENTAL", peak_calling="generated intervals; no experimental evidence"),
    ])
    project["metadata"].update(species="Mus_musculus", synthetic=True,
        purpose="Generated integration fixture only; no biological interpretation or validation claim.",
        synthetic_clip_notice="CLIP-shaped records are generated integration inputs, not observed binding.")
    project_path.write_text(json.dumps(project, indent=2), encoding="utf-8")
    suggestion = {"schema": "katana.flow.config/v1", "contrast_ids": list(contrasts), "design_factors": ["condition"],
                  "parameters": {"prioritize.motiftrace_consequences": {"custom_motif_path": str(motifs), "cluster_permutations": 99}},
                  "fixture_scope": "128 generated events, nine independent generated samples, two six-sample comparisons, twelve distinct motif models"}
    (root / "suggested-config.json").write_text(json.dumps(suggestion, indent=2), encoding="utf-8")
    return project_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    args = parser.parse_args()
    print(create_generic_fixture(args.outdir))
