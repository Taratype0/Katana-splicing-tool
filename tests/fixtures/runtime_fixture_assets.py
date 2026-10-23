"""Add generated BAM, sequencing-QC and exon-usage inputs before fixture sealing.

These records test software integration; QC and DEXSeq values are deliberately
synthetic input observations, never estimates from a real sequencing experiment.
No workflow is executed and no existing result is recomputed by this helper.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

import pandas as pd
import pysam


CHROMOSOMES = ("chr1", "chr2", "chr3", "chr4")
CHROMOSOME_LENGTH = 104_000


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _asset(identity, kind, path, **metadata):
    return {"asset_id": identity, "kind": kind, "path": str(path), "required": True,
            "metadata": {"synthetic": True, "source": "generated software fixture only",
                         "upstream_sha256": sha256(path), **metadata}}


def _bam(path, sample_id, fasta, sample_index):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": chrom, "LN": CHROMOSOME_LENGTH} for chrom in CHROMOSOMES],
              "RG": [{"ID": sample_id, "SM": sample_id, "LB": "synthetic", "PL": "ILLUMINA"}]}
    count = 0
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        for reference_id, chrom in enumerate(CHROMOSOMES):
            # All reference intervals have reads; this also covers annotated
            # chr1 exon bodies used by the geneBody_coverage runtime check.
            for start in range(sample_index % 5, CHROMOSOME_LENGTH - 75, 50):
                read = pysam.AlignedSegment(output.header)
                read.query_name = f"synthetic_{sample_id}_{count}"
                read.query_sequence = fasta.fetch(chrom, start, start + 75)
                read.flag = 0
                read.reference_id = reference_id
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = ((0, 75),)
                read.query_qualities = pysam.qualitystring_to_array("I" * 75)
                read.set_tags([("NH", 1), ("NM", 0), ("RG", sample_id)])
                output.write(read)
                count += 1
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path), "rb") as check:
        assert check.references == CHROMOSOMES
        assert check.lengths == (CHROMOSOME_LENGTH,) * 4
        assert check.has_index() and check.mapped == count
        depth = sum(sum(base) for base in check.count_coverage("chr1", 1200, 1299))
        assert depth > 0, "Generated chr1 annotation region must have nonzero coverage"
    return {"sample_id": sample_id, "read_count": count, "chr1_1200_1299_base_coverage": depth,
            "bam_sha256": sha256(path), "index_sha256": sha256(str(path) + ".bai")}


def _qc(directory, sample_id, read_count):
    directory.mkdir()
    assets = []
    for phase, mean in (("raw", 37.5), ("trimmed", 38.5)):
        path = directory / f"{sample_id}_{phase}_fastqc.zip"
        body = (
            "##FastQC\t0.12.1\n>>Basic Statistics\tpass\n#Measure\tValue\n"
            f"Filename\t{sample_id}_{phase}.synthetic.fastq.gz\nTotal Sequences\t{read_count}\n"
            "Sequence length\t75\n%GC\t50\n>>END_MODULE\n"
            ">>Per base sequence quality\tpass\n#Base\tMean\tMedian\n"
            + "".join(f"{position}\t{mean - position / 200}\t{mean}\n" for position in (1, 15, 30, 45, 60, 75))
            + ">>END_MODULE\n>>Sequence Length Distribution\tpass\n#Length\tCount\n"
            f"75\t{read_count}\n>>END_MODULE\n"
        )
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            member = zipfile.ZipInfo(f"{sample_id}_{phase}/fastqc_data.txt", date_time=(2026, 1, 1, 0, 0, 0))
            member.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(member, body)
        assets.append(_asset(f"generated-qc::{sample_id}::{phase}", "sequencing_qc", path,
                             sample_id=sample_id, qc_phase=phase))
    trim = directory / f"{sample_id}.fastq.gz_trimming_report.txt"
    trim.write_text(f"Total reads processed: {read_count:,}\nReads with adapters: 0 (0.0%)\n"
                    f"Reads written (passing filters): {read_count:,} (100.0%)\n", encoding="utf-8")
    assets.append(_asset(f"generated-qc::{sample_id}::trimming", "sequencing_qc", trim, sample_id=sample_id))
    star = directory / f"{sample_id}Log.final.out"
    star.write_text(f"Number of input reads | {read_count}\nUniquely mapped reads % | 100%\n"
                    "% of reads mapped to multiple loci | 0%\n% of reads mapped to too many loci | 0%\n"
                    "% of reads unmapped: too many mismatches | 0%\n% of reads unmapped: too short | 0%\n"
                    "% of reads unmapped: other | 0%\nAverage input read length | 75\n", encoding="utf-8")
    assets.append(_asset(f"generated-qc::{sample_id}::alignment", "sequencing_qc", star, sample_id=sample_id))
    return assets


def _exons(gtf):
    rows, seen, counts = [], set(), {}
    for line in Path(gtf).read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 9 or columns[2] != "exon":
            continue
        gene = re.search(r'gene_id "([^"]+)"', columns[8]).group(1)
        key = (gene, columns[0], columns[3], columns[4], columns[6])
        if key in seen:
            continue
        seen.add(key)
        counts[gene] = counts.get(gene, 0) + 1
        rows.append({"groupID": gene, "featureID": f"E{counts[gene]:03}",
                     "genomicData.seqnames": columns[0], "genomicData.start": int(columns[3]),
                     "genomicData.end": int(columns[4]), "genomicData.strand": columns[6]})
    return rows


def add_runtime_assets(project_path):
    path = Path(project_path).resolve(strict=True)
    project = json.loads(path.read_text(encoding="utf-8"))
    if project.get("project_id") != "generic-pipeline-generated-acceptance" or project.get("metadata", {}).get("synthetic") is not True:
        raise ValueError("Only the explicitly generated generic-pipeline fixture may be extended")
    if len(project["samples"]) != 9 or len(project["contrasts"]) != 2:
        raise ValueError("This fixture helper expects nine generated samples and two comparisons")
    destination = path.parent / "inputs" / "generated-runtime-assets"
    before = sha256(path)
    bam_records, added = [], []
    with pysam.FastaFile(project["annotation"]["genome_fasta"]) as fasta:
        if tuple(fasta.references) != CHROMOSOMES or tuple(fasta.lengths) != (CHROMOSOME_LENGTH,) * 4:
            raise ValueError("Generated reference must contain four 104000-base chromosomes")
        destination.mkdir(exist_ok=False)
        for number, sample in enumerate(project["samples"]):
            sample_id = sample["sample_id"]
            bam = destination / f"{sample_id}.bam"
            record = _bam(bam, sample_id, fasta, number)
            bam_records.append(record)
            sample.setdefault("assets", {}).update(bam=str(bam), bai=str(bam) + ".bai")
            sample.setdefault("attributes", {}).update(synthetic="true", source="fastq")
            added.extend(_qc(destination / f"qc-{sample_id}", sample_id, record["read_count"]))
    exon_rows = _exons(project["annotation"]["gtf"])
    dexseq_records = []
    for index, contrast in enumerate(project["contrasts"]):
        comparison = contrast["contrast_id"]
        frame = pd.DataFrame(exon_rows)
        frame["exonBaseMean"] = [100 + i % 29 for i in range(len(frame))]
        frame["dispersion"] = .1
        frame["stat"] = [12.0 if i < 24 else .2 for i in range(len(frame))]
        frame["pvalue"] = [.0002 * (i + 1) if i < 24 else .2 + (i % 20) / 40 for i in range(len(frame))]
        frame["padj"] = [.001 * (i + 1) if i < 24 else .5 + (i % 20) / 40 for i in range(len(frame))]
        frame["log2fold_control_treatment"] = [(-1 if i % 2 else 1) * (1 + index * .1 + (i % 7) / 10) for i in range(len(frame))]
        output = destination / f"{comparison}.dexseq.tsv"
        frame.to_csv(output, sep="\t", index=False)
        added.append(_asset(f"generated-dexseq::{comparison}", "dexseq_result", output,
                            contrast_id=comparison, native_fold_change_direction="control/treatment",
                            native_pvalue_column="pvalue", native_fdr_column="padj",
                            statistical_source="Generated test values; no differential fit or biological inference"))
        dexseq_records.append({"contrast_id": comparison, "rows": len(frame), "sha256": sha256(output)})
    project["assets"].extend(added)
    project["metadata"]["runtime_fixture_inputs"] = {
        "synthetic": True, "purpose": "Generated BAM/QC/exon-usage integration inputs only",
        "qc_statistics_source": "Generated FastQC/STAR-shaped records, not upstream tool execution",
        "dexseq_statistics_source": "Generated values for parser/report validation, not a new statistical fit"}
    from katana_next.domain.models import ProjectManifest
    ProjectManifest.model_validate(project)
    path.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    receipt = {"schema": "katana.generated_runtime_assets/v1", "status": "created", "synthetic": True,
               "scientific_workflow_started": False, "source_manifest_before_sha256": before,
               "source_manifest_after_sha256": sha256(path), "bam_records": bam_records,
               "sequencing_qc_assets": sum(a["kind"] == "sequencing_qc" for a in added),
               "dexseq_records": dexseq_records,
               "files": [{"path": p.relative_to(destination).as_posix(), "sha256": sha256(p), "bytes": p.stat().st_size}
                         for p in sorted(destination.rglob("*")) if p.is_file()]}
    receipt_path = destination / "generation.receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    args = parser.parse_args()
    print(add_runtime_assets(args.project))
