from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import tarfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/"
    "gencode.v41.annotation.gtf.gz"
)
GENOME_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/"
    "GRCh38.primary_assembly.genome.fa.gz"
)
STAR_URL = "https://codeload.github.com/alexdobin/STAR/tar.gz/refs/tags/2.7.11b"
PRIMARY_CONDITIONS = {"WT_DMSO", "KO_DMSO"}


def _sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _soft_samples(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                accession = line.split("=", 1)[1].strip()
                current = {"geo_accession": accession}
                records[accession] = current
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split("=", 1)[1].strip()
                if ":" in value:
                    key, item = value.split(":", 1)
                    current[key.strip().casefold()] = item.strip()
    return records


def _title_fields(title: str) -> tuple[str, str, str]:
    if title.startswith("786-O SETD2 "):
        genotype = "WT" if " SETD2 WT_" in title else "KO"
        treatment = "A1331852" if "_A1331852 " in title else "DMSO"
        return "786-O", genotype, treatment
    if title.startswith("A498_SETD2 addback"):
        return "A498", "WT_addback", "None"
    if title.startswith("A498_SETD2 deficient"):
        return "A498", "SETD2_deficient", "None"
    raise ValueError(f"Unrecognized GSE297044 sample title: {title}")


def build_manifest(source_root: Path, work_root: Path) -> pd.DataFrame:
    runinfo = pd.read_csv(source_root / "ena_runinfo.tsv", sep="\t", dtype=str)
    soft = _soft_samples(source_root / "GSE297044_family.soft.gz")
    rows: list[dict[str, object]] = []
    for row in runinfo.itertuples(index=False):
        title_match = re.search(r"GSM\d+: (.+?) Homo sapiens RNA-Seq", row.experiment_title)
        geo_match = re.search(r"(GSM\d+):", row.experiment_title)
        if not title_match or not geo_match:
            raise ValueError(f"Could not parse ENA experiment title: {row.experiment_title}")
        title = title_match.group(1)
        geo_accession = geo_match.group(1)
        cell_line, genotype, treatment = _title_fields(title)
        condition = f"{genotype}_{treatment}" if treatment != "None" else genotype
        role = (
            "primary_genotype"
            if cell_line == "786-O" and treatment == "DMSO"
            else "secondary_treatment"
            if cell_line == "786-O"
            else "secondary_addback"
        )
        urls = [f"https://{value}" for value in str(row.fastq_ftp).split(";")]
        byte_values = [int(value) for value in str(row.fastq_bytes).split(";")]
        record = soft.get(geo_accession, {})
        rows.append(
            {
                "sample_id": row.run_accession,
                "run_accession": row.run_accession,
                "geo_accession": geo_accession,
                "title": title,
                "cell_line": cell_line,
                "genotype": genotype,
                "treatment": treatment,
                "condition": condition,
                "batch": record.get("batch", ""),
                "validation_role": role,
                "library_layout": row.library_layout,
                "read1_url": urls[0],
                "read2_url": urls[1],
                "read1_bytes": byte_values[0],
                "read2_bytes": byte_values[1],
                "total_fastq_bytes": sum(byte_values),
                "genome_build": "GRCh38",
                "annotation": "GENCODE v41",
            }
        )
    manifest = pd.DataFrame(rows).sort_values(
        ["validation_role", "cell_line", "condition", "sample_id"]
    )
    inputs = work_root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(inputs / "sample_manifest.tsv", sep="\t", index=False)
    primary = manifest.loc[manifest["validation_role"].eq("primary_genotype")].copy()
    primary[["sample_id", "condition", "genotype", "treatment", "batch"]].to_csv(
        inputs / "sample_metadata_primary.tsv", sep="\t", index=False
    )
    audit = pd.DataFrame(
        [
            {"check": "organism", "value": "Homo sapiens", "passed": True},
            {"check": "primary_conditions", "value": ";".join(sorted(primary["condition"].unique())), "passed": set(primary["condition"]) == PRIMARY_CONDITIONS},
            {"check": "primary_replicates", "value": json.dumps(primary.groupby("condition")["sample_id"].nunique().to_dict(), sort_keys=True), "passed": bool(primary.groupby("condition")["sample_id"].nunique().eq(3).all())},
            {"check": "paired_end", "value": str(primary["library_layout"].unique().tolist()), "passed": bool(primary["library_layout"].eq("PAIRED").all())},
            {"check": "supplementary_count_semantics", "value": "featureCounts gene-assigned exon aggregation; not junction IJC/SJC", "passed": False},
            {"check": "primary_fastq_GB", "value": round(primary["total_fastq_bytes"].sum() / 1024**3, 3), "passed": True},
        ]
    )
    audit.to_csv(inputs / "input_audit.tsv", sep="\t", index=False)
    preregistration = {
        "dataset": "GSE297044 / PRJNA1262435",
        "primary_contrast": "786-O SETD2 KO DMSO vs 786-O SETD2 WT DMSO",
        "primary_samples": primary["sample_id"].tolist(),
        "replicates_per_group": 3,
        "primary_endpoint": "locked human ccRCC program projection and genotype association",
        "secondary_endpoints": [
            "A1331852 treatment response in 786-O WT and KO",
            "A498 SETD2 add-back direction sensitivity analysis",
        ],
        "event_evidence": "STAR BAM splice junction counts; rMATS annotation optional and separate",
        "external_validation_boundary": (
            "Six primary samples are not used for de novo parameter tuning or claims of "
            "standalone network-discovery superiority. SpliceSeq-to-BAM comparisons are "
            "reported as exact event only when coordinates can be harmonized; otherwise "
            "they are explicitly gene/program-level projections."
        ),
        "frozen_before_alignment": True,
    }
    (inputs / "preregistration.json").write_text(
        json.dumps(preregistration, indent=2), encoding="utf-8"
    )
    return manifest


def _download(url: str, target: Path, expected_bytes: int | None = None) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        if expected_bytes is not None and target.stat().st_size != expected_bytes:
            raise ValueError(
                f"Existing download has the wrong size: {target} "
                f"({target.stat().st_size} != {expected_bytes})"
            )
        return
    partial = target.with_suffix(target.suffix + ".part")
    subprocess.run(
        [
            "curl.exe",
            "--fail",
            "--location",
            "--retry",
            "8",
            "--retry-all-errors",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            url,
        ],
        check=True,
    )
    if expected_bytes is not None and partial.stat().st_size != expected_bytes:
        raise ValueError(
            f"Downloaded file has the wrong size: {partial} "
            f"({partial.stat().st_size} != {expected_bytes})"
        )
    partial.replace(target)


def _gunzip(source: Path, target: Path) -> None:
    if target.exists() and target.stat().st_size > 0:
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    with gzip.open(source, "rb") as reader, partial.open("wb") as writer:
        shutil.copyfileobj(reader, writer, length=8 * 1024 * 1024)
    partial.replace(target)


def prepare_references(work_root: Path) -> tuple[Path, Path, Path]:
    reference = work_root / "reference"
    gtf_gz = reference / "gencode.v41.annotation.gtf.gz"
    fasta_gz = reference / "GRCh38.primary_assembly.genome.fa.gz"
    _download(GTF_URL, gtf_gz)
    _download(GENOME_URL, fasta_gz)
    gtf = reference / "gencode.v41.annotation.gtf"
    fasta = reference / "GRCh38.primary_assembly.genome.fa"
    _gunzip(gtf_gz, gtf)
    _gunzip(fasta_gz, fasta)

    archive = work_root / "tools" / "STAR-2.7.11b.tar.gz"
    _download(STAR_URL, archive)
    tool_root = work_root / "tools" / "STAR-2.7.11b"
    star = tool_root / "bin" / "Linux_x86_64_static" / "STAR"
    if not star.exists():
        tool_root.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "r:gz") as handle:
            try:
                handle.extractall(tool_root.parent, filter="data")
            except TypeError:  # Python 3.11 builds without extraction filters.
                handle.extractall(tool_root.parent)
        extracted = tool_root.parent / "STAR-2.7.11b"
        if extracted != tool_root and extracted.exists():
            extracted.replace(tool_root)
    if not star.exists():
        raise FileNotFoundError(f"STAR binary was not found after extraction: {star}")
    return gtf, fasta, star


def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    return f"/mnt/{drive}/{resolved.as_posix().split(':', 1)[1].lstrip('/')}"


def _run_wsl(command: Iterable[str], log: Path | None = None) -> None:
    shell_command = " ".join(shlex.quote(str(value)) for value in command)
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        shell_command += f" > {shlex.quote(_wsl_path(log))} 2>&1"
    subprocess.run(["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", shell_command], check=True)


def build_star_index(work_root: Path, gtf: Path, fasta: Path, star: Path, threads: int) -> Path:
    index = work_root / "reference" / "STAR_index_v41"
    marker = index / ".katana_complete.json"
    required = [index / "Genome", index / "SA", index / "SAindex"]
    genome_sa_sparse_d = 2
    required_complete = all(path.exists() and path.stat().st_size > 0 for path in required)
    if required_complete:
        if not marker.exists():
            marker.write_text(
                json.dumps(
                    {
                        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                        "adopted_complete_index": True,
                        "STAR": str(star.resolve()),
                        "genome_fasta": str(fasta.resolve()),
                        "gtf": str(gtf.resolve()),
                        "sjdbOverhang": 149,
                        "genomeSAsparseD": genome_sa_sparse_d,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
        return index
    index.mkdir(parents=True, exist_ok=True)
    _run_wsl(
        [
            _wsl_path(star),
            "--runMode", "genomeGenerate",
            "--runThreadN", str(threads),
            "--genomeDir", _wsl_path(index),
            "--genomeFastaFiles", _wsl_path(fasta),
            "--sjdbGTFfile", _wsl_path(gtf),
            "--sjdbOverhang", "149",
            # Sparse suffix-array indexing keeps GRCh38 within common WSL
            # desktop memory limits without changing aligned read coordinates.
            "--genomeSAsparseD", str(genome_sa_sparse_d),
        ],
        work_root / "logs" / "STAR_genomeGenerate.log",
    )
    missing = [str(path) for path in required if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"STAR index finished without required files: {missing}")
    marker.write_text(
        json.dumps(
            {
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                "STAR": str(star.resolve()),
                "genome_fasta": str(fasta.resolve()),
                "gtf": str(gtf.resolve()),
                "sjdbOverhang": 149,
                "genomeSAsparseD": genome_sa_sparse_d,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return index


def download_fastq(
    manifest: pd.DataFrame,
    work_root: Path,
    subset: str,
    workers: int,
    sample_ids: list[str] | None = None,
) -> pd.DataFrame:
    selected = _select_subset(manifest, subset, sample_ids)
    jobs: list[tuple[str, Path, int]] = []
    for row in selected.itertuples(index=False):
        jobs.extend(
            [
                (row.read1_url, work_root / "fastq" / f"{row.sample_id}_1.fastq.gz", int(row.read1_bytes)),
                (row.read2_url, work_root / "fastq" / f"{row.sample_id}_2.fastq.gz", int(row.read2_bytes)),
            ]
        )
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = [
            executor.submit(_download, url, target, expected_bytes)
            for url, target, expected_bytes in jobs
        ]
        for future in futures:
            future.result()
    return selected


def _select_subset(
    manifest: pd.DataFrame,
    subset: str,
    sample_ids: list[str] | None = None,
) -> pd.DataFrame:
    if subset == "primary":
        selected = manifest.loc[manifest["validation_role"].eq("primary_genotype")].copy()
    elif subset == "786o":
        selected = manifest.loc[manifest["cell_line"].eq("786-O")].copy()
    elif subset == "all":
        selected = manifest.copy()
    else:
        raise ValueError(f"Unknown subset: {subset}")
    if sample_ids:
        missing = sorted(set(sample_ids).difference(selected["sample_id"].astype(str)))
        if missing:
            raise ValueError(f"Requested sample IDs are outside subset={subset}: {missing}")
        selected = selected.loc[selected["sample_id"].astype(str).isin(sample_ids)].copy()
    return selected


def align_samples(
    manifest: pd.DataFrame,
    work_root: Path,
    subset: str,
    index: Path,
    star: Path,
    threads: int,
    delete_fastq: bool,
    sample_ids: list[str] | None = None,
) -> None:
    for row in _select_subset(manifest, subset, sample_ids).itertuples(index=False):
        bam = work_root / "bam" / f"{row.sample_id}.bam"
        bai = Path(f"{bam}.bai")
        if bam.exists() and bai.exists():
            continue
        sample_dir = work_root / "alignment" / row.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        prefix = sample_dir / f"{row.sample_id}."
        read1 = work_root / "fastq" / f"{row.sample_id}_1.fastq.gz"
        read2 = work_root / "fastq" / f"{row.sample_id}_2.fastq.gz"
        if not read1.exists() or not read2.exists():
            raise FileNotFoundError(f"Missing FASTQ for {row.sample_id}")
        _run_wsl(
            [
                _wsl_path(star),
                "--runThreadN", str(threads),
                "--genomeDir", _wsl_path(index),
                "--readFilesIn", _wsl_path(read1), _wsl_path(read2),
                "--readFilesCommand", "zcat",
                "--twopassMode", "Basic",
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--limitBAMsortRAM", "4000000000",
                "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD", "XS",
                "--outFileNamePrefix", _wsl_path(prefix),
            ],
            work_root / "logs" / f"STAR_{row.sample_id}.log",
        )
        produced = Path(f"{prefix}Aligned.sortedByCoord.out.bam")
        bam.parent.mkdir(parents=True, exist_ok=True)
        produced.replace(bam)
        _run_wsl(["samtools", "index", "-@", str(max(1, threads // 2)), _wsl_path(bam)])
        if not bai.exists():
            raise FileNotFoundError(f"samtools index did not create {bai}")
        if delete_fastq:
            read1.unlink(missing_ok=True)
            read2.unlink(missing_ok=True)


def write_config(work_root: Path, gtf: Path, fasta: Path) -> Path:
    primary_samples = pd.read_csv(
        work_root / "inputs" / "sample_metadata_primary.tsv", sep="\t"
    )["sample_id"].astype(str)
    junction_root = work_root / "junction_counts"
    precomputed_complete = all(
        (junction_root / f"{sample_id}.junction_counts.tsv.gz").exists()
        for sample_id in primary_samples
    )
    config = {
        "project_name": "GSE297044_SETD2_ccRCC_primary",
        "rmats_root": None,
        "sample_metadata": str(work_root / "inputs" / "sample_metadata_primary.tsv"),
        "output_dir": str(work_root / "katana_cosplice_output"),
        "gtf": str(gtf),
        "genome_fasta": str(fasta),
        "bam_root": str(junction_root if precomputed_complete else work_root / "bam"),
        "genome_build": "GRCh38",
        "annotation_version": "GENCODE v41",
        "selected_conditions": ["WT_DMSO", "KO_DMSO"],
        "association_factors": ["genotype"],
        "batch_column": "batch",
        "matrix_input": {
            "input_mode": "bam_junctions",
            "backend": "precomputed_junctions" if precomputed_complete else "wsl_pysam",
            "min_mapping_quality": 1,
            "min_anchor_bp": 8,
            "max_intron_length": 500000,
            "min_total_junction_reads": 20,
            "min_samples_with_reads": 2,
            "min_cluster_fraction": 0.001,
            "exclude_duplicates": True,
            "max_parallel_samples": 2,
        },
        "stable_as": {
            "preset": "strict",
            "min_total_junction_reads": 20,
            "min_supported_replicates_per_condition": 3,
        },
        "network": {
            "random_seed": 20260901,
            "engine": "tom_hierarchical",
            "module_detection": "wgcna_blockwise",
        },
        "differential_network": {"enabled": False},
    }
    path = work_root / "cosplice_config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


def write_provenance(work_root: Path, manifest: pd.DataFrame) -> None:
    rows: list[dict[str, object]] = []
    provenance_files = [
        *sorted((work_root / "inputs").glob("*")),
        *sorted((work_root / "fastq").glob("*.fastq.gz")),
        *sorted((work_root / "reference").glob("*.gtf")),
        *sorted((work_root / "reference").glob("*.fa")),
        *sorted((work_root / "reference" / "hisat2_grch38").glob("**/*.ht2")),
        *sorted((work_root / "junction_counts").glob("*.junction_counts.tsv.gz")),
        *sorted((work_root / "junction_counts").glob("*.alignment_summary.json")),
    ]
    for path in provenance_files:
        if not path.is_file():
            continue
        rows.append({"path": str(path), "bytes": path.stat().st_size, "sha256": _sha256(path)})
    pd.DataFrame(rows).to_csv(work_root / "input_provenance.tsv", sep="\t", index=False)
    status = {
        "manifest_samples": len(manifest),
        "primary_samples": int(manifest["validation_role"].eq("primary_genotype").sum()),
        "primary_fastq_GB": float(
            manifest.loc[manifest["validation_role"].eq("primary_genotype"), "total_fastq_bytes"].sum()
            / 1024**3
        ),
        "reference": "GENCODE v41 / GRCh38 primary assembly",
    }
    (work_root / "preprocessing_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Prepare the preregistered GSE297044 SETD2 ccRCC external validation."
    )
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument(
        "--stage",
        choices=["manifest", "references", "index", "download", "align", "config", "all"],
        default="all",
    )
    parser.add_argument("--subset", choices=["primary", "786o", "all"], default="primary")
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--download-workers", type=int, default=4)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--delete-fastq-after-alignment", action="store_true")
    args = parser.parse_args(argv)

    source_root = args.source_root.resolve()
    work_root = args.work_root.resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    manifest_path = work_root / "inputs" / "sample_manifest.tsv"
    manifest = (
        pd.read_csv(manifest_path, sep="\t")
        if manifest_path.exists()
        else build_manifest(source_root, work_root)
    )
    if args.stage == "manifest":
        write_provenance(work_root, manifest)
        return 0

    gtf: Path | None = None
    fasta: Path | None = None
    star: Path | None = None
    if args.stage in {"references", "index", "align", "config", "all"}:
        gtf, fasta, star = prepare_references(work_root)
        if args.stage == "references":
            write_provenance(work_root, manifest)
            return 0
    if args.stage == "index":
        assert gtf is not None and fasta is not None and star is not None
        build_star_index(work_root, gtf, fasta, star, args.threads)
        write_provenance(work_root, manifest)
        return 0
    if args.stage in {"download", "all"}:
        download_fastq(
            manifest,
            work_root,
            args.subset,
            args.download_workers,
            args.sample_id,
        )
        if args.stage == "download":
            return 0
    if args.stage in {"align", "all"}:
        assert gtf is not None and fasta is not None and star is not None
        index = build_star_index(work_root, gtf, fasta, star, args.threads)
        align_samples(
            manifest,
            work_root,
            args.subset,
            index,
            star,
            args.threads,
            args.delete_fastq_after_alignment,
            args.sample_id,
        )
    if args.stage in {"config", "all"}:
        assert gtf is not None and fasta is not None
        write_config(work_root, gtf, fasta)
    write_provenance(work_root, manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
