from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import json
import os
from pathlib import Path
import subprocess
from typing import Callable, Iterable, Mapping

import numpy as np
import pandas as pd

from katana_cosplice.audit import load_sample_metadata, resolve_rmats_root
from katana_cosplice.config import CoSpliceConfig, MatrixInputConfig
from katana_cosplice.matrix import EVENT_TYPES, MatrixResult
from katana_cosplice.utils import file_signature, fingerprint, normalize_gene_id, normalize_sample_id, write_json


Progress = Callable[[str], None] | None
Junction = tuple[str, int, int]


@dataclass
class BamMatrixResult:
    matrix: MatrixResult
    junction_counts: pd.DataFrame
    clusters: pd.DataFrame
    rmats_mapping: pd.DataFrame


class _UnionFind:
    def __init__(self, values: Iterable[Junction]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: Junction) -> Junction:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: Junction, right: Junction) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def _emit(progress: Progress, message: str) -> None:
    if progress:
        progress(message)


def resolve_bam_sources(samples: pd.DataFrame, bam_root: Path) -> pd.DataFrame:
    if not bam_root.exists():
        raise FileNotFoundError(f"BAM root does not exist: {bam_root}")
    candidates = sorted(bam_root.glob("**/*.bam"))
    by_sample: dict[str, list[Path]] = defaultdict(list)
    for path in candidates:
        by_sample[normalize_sample_id(path.name)].append(path)

    rows: list[dict[str, object]] = []
    for row in samples.itertuples(index=False):
        sample_id = str(row.sample_id)
        matches = by_sample.get(sample_id, [])
        if not matches and hasattr(row, "genome_bam"):
            expected = normalize_sample_id(getattr(row, "genome_bam"))
            matches = by_sample.get(expected, [])
        if not matches:
            raise FileNotFoundError(f"No local BAM matched sample_id={sample_id} under {bam_root}")
        if len(matches) > 1:
            raise ValueError(f"Multiple BAM files matched sample_id={sample_id}: {matches}")
        bam = matches[0]
        bai_candidates = [Path(f"{bam}.bai"), bam.with_suffix(".bai")]
        bai = next((path for path in bai_candidates if path.exists()), None)
        if bai is None:
            raise FileNotFoundError(f"Missing BAM index for sample_id={sample_id}: {bam}")
        rows.append(
            {
                "sample_id": sample_id,
                "condition": str(row.condition),
                "bam_path": str(bam),
                "bai_path": str(bai),
                "bam_size": bam.stat().st_size,
            }
        )
    return pd.DataFrame(rows)


def _backend_name(requested: str) -> str:
    if requested in {"auto", "pysam"}:
        try:
            import pysam  # noqa: F401

            return "pysam"
        except ImportError:
            if requested == "pysam":
                raise RuntimeError("matrix_input.backend=pysam but pysam is not installed")
    if requested in {"auto", "wsl_pysam"} and os.name == "nt":
        check = subprocess.run(
            ["wsl.exe", "-d", "Ubuntu", "--", "python3", "-c", "import pysam"],
            capture_output=True,
            text=True,
            check=False,
        )
        if check.returncode == 0:
            return "wsl_pysam"
        if requested == "wsl_pysam":
            raise RuntimeError(
                "matrix_input.backend=wsl_pysam but Ubuntu WSL with python3-pysam is unavailable"
            )
    try:
        import bamnostic  # noqa: F401

        return "bamnostic"
    except ImportError as exc:
        raise RuntimeError("BAM junction mode requires pysam or bamnostic") from exc


def _record_junctions(record: object, settings: MatrixInputConfig) -> list[Junction]:
    if any(
        bool(getattr(record, attribute, False))
        for attribute in ("is_unmapped", "is_secondary", "is_supplementary", "is_qcfail")
    ):
        return []
    if settings.exclude_duplicates and bool(getattr(record, "is_duplicate", False)):
        return []
    mapq = int(getattr(record, "mapping_quality", getattr(record, "mapq", 0)) or 0)
    if mapq < settings.min_mapping_quality:
        return []
    chrom = str(getattr(record, "reference_name", "") or "")
    position = int(getattr(record, "reference_start", -1))
    operations = list(getattr(record, "cigartuples", None) or [])
    if not chrom or position < 0 or not operations:
        return []

    output: list[Junction] = []
    for index, (operation, length) in enumerate(operations):
        length = int(length)
        if operation == 3:
            left_anchor = 0
            for left_index in range(index - 1, -1, -1):
                left_operation, left_length = operations[left_index]
                if left_operation in {0, 7, 8}:
                    left_anchor = int(left_length)
                    break
                if left_operation in {3, 4, 5}:
                    break
            right_anchor = 0
            for right_index in range(index + 1, len(operations)):
                right_operation, right_length = operations[right_index]
                if right_operation in {0, 7, 8}:
                    right_anchor = int(right_length)
                    break
                if right_operation in {3, 4, 5}:
                    break
            if (
                min(left_anchor, right_anchor) >= settings.min_anchor_bp
                and 0 < length <= settings.max_intron_length
            ):
                output.append((chrom, position, position + length))
            position += length
        elif operation in {0, 2, 7, 8}:
            position += length
    return output


def _iter_bam_records(path: Path, backend: str) -> Iterable[object]:
    if backend == "pysam":
        import pysam

        with pysam.AlignmentFile(str(path), "rb") as handle:
            yield from handle.fetch(until_eof=True)
        return

    import bamnostic

    with bamnostic.AlignmentFile(str(path), "rb") as handle:
        yield from handle


def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"WSL backend requires a Windows drive path: {resolved}")
    suffix = resolved.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{suffix}"


def _count_with_wsl_pysam(
    bam: Path,
    output: Path,
    summary: Path,
    settings: MatrixInputConfig,
) -> dict[str, object]:
    helper = Path(__file__).with_name("wsl_count_junctions.py")
    command = [
        "wsl.exe",
        "-d",
        "Ubuntu",
        "--",
        "python3",
        _wsl_path(helper),
        "--bam",
        _wsl_path(bam),
        "--output",
        _wsl_path(output),
        "--summary",
        _wsl_path(summary),
        "--min-mapq",
        str(settings.min_mapping_quality),
        "--min-anchor",
        str(settings.min_anchor_bp),
        "--max-intron",
        str(settings.max_intron_length),
        "--scratch-copy",
    ]
    if not settings.exclude_duplicates:
        command.append("--include-duplicates")
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "WSL pysam junction counter failed: "
            f"return_code={completed.returncode}; stderr={completed.stderr.strip()}"
        )
    return json.loads(summary.read_text(encoding="utf-8"))


def _sample_cache_fingerprint(bam: Path, bai: Path, backend: str, settings: MatrixInputConfig) -> str:
    return fingerprint(
        {
            "schema": 1,
            "bam": file_signature(bam),
            "bai": file_signature(bai),
            "backend": backend,
            "settings": settings.__dict__,
        }
    )


def count_sample_junctions(
    bam: Path,
    bai: Path,
    settings: MatrixInputConfig,
    backend: str,
    cache_dir: Path,
    sample_id: str,
    progress: Progress = None,
) -> tuple[Counter[Junction], dict[str, object]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    counts_path = cache_dir / f"{sample_id}.junction_counts.tsv.gz"
    parameters_path = cache_dir / f"{sample_id}.parameters.json"
    cache_fingerprint = _sample_cache_fingerprint(bam, bai, backend, settings)
    if counts_path.exists() and parameters_path.exists():
        parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
        if parameters.get("fingerprint") == cache_fingerprint:
            cached = pd.read_csv(counts_path, sep="\t")
            counts = Counter(
                {
                    (str(row.chrom), int(row.intron_start_0based), int(row.intron_end_0based)): int(row.read_count)
                    for row in cached.itertuples(index=False)
                }
            )
            parameters["cache_status"] = "cache_hit"
            _emit(progress, f"BAM junction cache hit: {sample_id} ({len(counts):,} junctions)")
            return counts, parameters

    _emit(progress, f"Scanning BAM junctions: {sample_id} with {backend}")
    if backend == "wsl_pysam":
        summary_path = cache_dir / f"{sample_id}.wsl_summary.json"
        summary = _count_with_wsl_pysam(
            bam=bam,
            output=counts_path,
            summary=summary_path,
            settings=settings,
        )
        cached = pd.read_csv(counts_path, sep="\t")
        counts = Counter(
            {
                (str(row.chrom), int(row.intron_start_0based), int(row.intron_end_0based)): int(row.read_count)
                for row in cached.itertuples(index=False)
            }
        )
        summary_path.unlink(missing_ok=True)
        parameters = {
            "fingerprint": cache_fingerprint,
            "cache_status": "rebuilt",
            "sample_id": sample_id,
            "backend": backend,
            "bam_path": str(bam),
            "bai_path": str(bai),
            **summary,
        }
        write_json(parameters_path, parameters)
        _emit(
            progress,
            f"Finished BAM junction scan: {sample_id} "
            f"({len(counts):,} junctions; {sum(counts.values()):,} junction observations; "
            f"{float(summary.get('elapsed_seconds', 0)):.1f}s)",
        )
        return counts, parameters

    counts: Counter[Junction] = Counter()
    alignments_seen = 0
    spliced_alignments = 0
    for record in _iter_bam_records(bam, backend):
        alignments_seen += 1
        junctions = _record_junctions(record, settings)
        if junctions:
            spliced_alignments += 1
            counts.update(junctions)
    table = pd.DataFrame(
        [
            {
                "chrom": chrom,
                "intron_start_0based": start,
                "intron_end_0based": end,
                "read_count": count,
            }
            for (chrom, start, end), count in counts.items()
        ]
    )
    if table.empty:
        table = pd.DataFrame(
            columns=["chrom", "intron_start_0based", "intron_end_0based", "read_count"]
        )
    table.to_csv(counts_path, sep="\t", index=False, compression="gzip")
    parameters = {
        "fingerprint": cache_fingerprint,
        "cache_status": "rebuilt",
        "sample_id": sample_id,
        "backend": backend,
        "bam_path": str(bam),
        "bai_path": str(bai),
        "alignments_seen": alignments_seen,
        "spliced_alignments": spliced_alignments,
        "unique_junctions": len(counts),
        "junction_read_observations": int(sum(counts.values())),
    }
    write_json(parameters_path, parameters)
    _emit(
        progress,
        f"Finished BAM junction scan: {sample_id} "
        f"({len(counts):,} junctions; {sum(counts.values()):,} junction observations)",
    )
    return counts, parameters


def _parse_gtf_attributes(value: str) -> dict[str, str]:
    output: dict[str, str] = {}
    for item in value.split(";"):
        item = item.strip()
        if " " not in item:
            continue
        key, raw = item.split(" ", 1)
        output[key] = raw.strip().strip('"')
    return output


def annotate_junctions_with_gtf(junctions: set[Junction], gtf: Path) -> pd.DataFrame:
    if not gtf.exists():
        raise FileNotFoundError(f"GTF does not exist: {gtf}")
    starts_by_chrom: dict[str, set[int]] = defaultdict(set)
    ends_by_chrom: dict[str, set[int]] = defaultdict(set)
    for chrom, start, end in junctions:
        starts_by_chrom[chrom].add(start)
        ends_by_chrom[chrom].add(end)

    left_annotations: dict[tuple[str, int], list[tuple[str, str, str, str]]] = defaultdict(list)
    right_annotations: dict[tuple[str, int], list[tuple[str, str, str, str]]] = defaultdict(list)
    with gtf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom = fields[0]
            if chrom not in starts_by_chrom and chrom not in ends_by_chrom:
                continue
            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            use_left = end0 in starts_by_chrom.get(chrom, set())
            use_right = start0 in ends_by_chrom.get(chrom, set())
            if not use_left and not use_right:
                continue
            attributes = _parse_gtf_attributes(fields[8])
            annotation = (
                normalize_gene_id(attributes.get("gene_id", "")),
                attributes.get("gene_name", ""),
                attributes.get("transcript_id", ""),
                fields[6],
            )
            if use_left:
                left_annotations[(chrom, end0)].append(annotation)
            if use_right:
                right_annotations[(chrom, start0)].append(annotation)

    rows: list[dict[str, object]] = []
    for chrom, start, end in sorted(junctions):
        left = left_annotations.get((chrom, start), [])
        right = right_annotations.get((chrom, end), [])
        left_genes = {item[0] for item in left if item[0]}
        right_genes = {item[0] for item in right if item[0]}
        genes = sorted(left_genes.intersection(right_genes))
        symbols = sorted(
            {
                item[1]
                for item in [*left, *right]
                if item[0] in genes and item[1]
            }
        )
        left_transcripts = {item[2] for item in left if item[2]}
        right_transcripts = {item[2] for item in right if item[2]}
        direct_transcripts = sorted(left_transcripts.intersection(right_transcripts))
        strands = sorted(
            {
                item[3]
                for item in [*left, *right]
                if item[0] in genes and item[3] in {"+", "-"}
            }
        )
        rows.append(
            {
                "chrom": chrom,
                "intron_start_0based": start,
                "intron_end_0based": end,
                "gene_id": ";".join(genes),
                "gene_symbol": ";".join(symbols),
                "strand": strands[0] if len(strands) == 1 else ".",
                "annotated_splice_boundaries": bool(left and right),
                "annotated_direct_junction": bool(direct_transcripts),
                "annotated_transcript_count": len(direct_transcripts),
                "annotated_transcript_ids": ";".join(direct_transcripts),
            }
        )
    return pd.DataFrame(rows)


def _cluster_junctions(junctions: set[Junction]) -> dict[Junction, str]:
    union_find = _UnionFind(junctions)
    donor_lookup: dict[tuple[str, int], Junction] = {}
    acceptor_lookup: dict[tuple[str, int], Junction] = {}
    for junction in sorted(junctions):
        chrom, start, end = junction
        for lookup, key in (
            (donor_lookup, (chrom, start)),
            (acceptor_lookup, (chrom, end)),
        ):
            previous = lookup.get(key)
            if previous is None:
                lookup[key] = junction
            else:
                union_find.union(previous, junction)
    groups: dict[Junction, list[Junction]] = defaultdict(list)
    for junction in junctions:
        groups[union_find.find(junction)].append(junction)
    cluster_map: dict[Junction, str] = {}
    for index, members in enumerate(
        sorted(groups.values(), key=lambda values: min(values)),
        start=1,
    ):
        if len(members) < 2:
            continue
        cluster_id = f"JCL{index:07d}"
        for junction in members:
            cluster_map[junction] = cluster_id
    return cluster_map


def _filter_and_cluster(
    counts: pd.DataFrame,
    settings: MatrixInputConfig,
) -> tuple[pd.DataFrame, dict[Junction, str]]:
    total = counts.sum(axis=1)
    supported_samples = counts.gt(0).sum(axis=1)
    keep = total.ge(settings.min_total_junction_reads) & supported_samples.ge(
        settings.min_samples_with_reads
    )
    filtered = counts.loc[keep].copy()
    junctions = {
        (str(index[0]), int(index[1]), int(index[2]))
        for index in filtered.index
    }
    cluster_map = _cluster_junctions(junctions)
    filtered = filtered.loc[[index in cluster_map for index in filtered.index]]
    if filtered.empty:
        return filtered, {}

    cluster_series = pd.Series(
        [cluster_map[tuple(index)] for index in filtered.index],
        index=filtered.index,
        dtype="string",
    )
    junction_totals = filtered.sum(axis=1).astype(float)
    cluster_totals = junction_totals.groupby(cluster_series, observed=True).transform("sum")
    fractions = junction_totals.divide(cluster_totals.clip(lower=1.0))
    filtered = filtered.loc[fractions.ge(settings.min_cluster_fraction)]
    final_junctions = {
        (str(index[0]), int(index[1]), int(index[2]))
        for index in filtered.index
    }
    final_clusters = _cluster_junctions(final_junctions)
    filtered = filtered.loc[[index in final_clusters for index in filtered.index]]
    return filtered, final_clusters


def _rmats_row_junctions(
    row: Mapping[str, object],
    event_type: str,
) -> list[tuple[int, int, str]]:
    def value(name: str) -> int:
        return int(row[name])

    if event_type == "SE":
        return [
            (value("upstreamEE"), value("exonStart_0base"), "upstream_inclusion"),
            (value("exonEnd"), value("downstreamES"), "downstream_inclusion"),
            (value("upstreamEE"), value("downstreamES"), "skipping"),
        ]
    if event_type == "RI":
        return [(value("upstreamEE"), value("downstreamES"), "retained_intron")]
    if event_type == "MXE":
        return [
            (value("upstreamEE"), value("1stExonStart_0base"), "first_exon_upstream"),
            (value("1stExonEnd"), value("downstreamES"), "first_exon_downstream"),
            (value("upstreamEE"), value("2ndExonStart_0base"), "second_exon_upstream"),
            (value("2ndExonEnd"), value("downstreamES"), "second_exon_downstream"),
        ]
    intervals = [
        (value("longExonStart_0base"), value("longExonEnd"), "long"),
        (value("shortES"), value("shortEE"), "short"),
        (value("flankingES"), value("flankingEE"), "flanking"),
    ]
    output: list[tuple[int, int, str]] = []
    for left_index, left in enumerate(intervals):
        for right in intervals[left_index + 1 :]:
            if left[1] <= right[0]:
                output.append((left[1], right[0], f"{left[2]}_{right[2]}"))
            elif right[1] <= left[0]:
                output.append((right[1], left[0], f"{right[2]}_{left[2]}"))
    return output


def map_rmats_to_junctions(
    junctions: set[Junction],
    rmats_root: Path,
    mode: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    junction_lookup = set(junctions)
    for path in sorted(rmats_root.glob(f"*/rmats_post/*.MATS.{mode}.txt")):
        event_type = path.name.split(".", 1)[0]
        if event_type not in EVENT_TYPES:
            continue
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        comparison_id = path.parent.parent.name
        for row in frame.to_dict(orient="records"):
            chrom = str(row["chr"])
            row_junctions = _rmats_row_junctions(row, event_type)
            event_id = (
                f"{event_type}|{chrom}|{row['strand']}|"
                + "|".join(str(value) for value in row_junctions)
            )
            for start, end, role in row_junctions:
                junction = (chrom, start, end)
                if junction not in junction_lookup:
                    continue
                rows.append(
                    {
                        "event_key": f"JUNC|{chrom}|{start}|{end}",
                        "comparison_id": comparison_id,
                        "rmats_event_type": event_type,
                        "rmats_event_id": event_id,
                        "rmats_gene_id": (
                            normalize_gene_id(str(row.get("GeneID")))
                            if pd.notna(row.get("GeneID"))
                            else ""
                        ),
                        "rmats_gene_symbol": (
                            str(row.get("geneSymbol"))
                            if pd.notna(row.get("geneSymbol"))
                            else ""
                        ),
                        "junction_role": role,
                        "rMATS_FDR": pd.to_numeric(row.get("FDR", np.nan), errors="coerce"),
                        "raw_dPSI": pd.to_numeric(
                            row.get("IncLevelDifference", np.nan), errors="coerce"
                        ),
                    }
                )
    columns = [
        "event_key",
        "comparison_id",
        "rmats_event_type",
        "rmats_event_id",
        "rmats_gene_id",
        "rmats_gene_symbol",
        "junction_role",
        "rMATS_FDR",
        "raw_dPSI",
    ]
    return pd.DataFrame(rows, columns=columns).drop_duplicates()


def build_bam_junction_matrices(
    config: CoSpliceConfig,
    cache_dir: Path,
    progress: Progress = None,
) -> BamMatrixResult:
    if config.bam_root is None:
        raise ValueError("matrix_input.input_mode=bam_junctions requires bam_root")
    if config.gtf is None:
        raise ValueError("matrix_input.input_mode=bam_junctions requires gtf")

    samples = load_sample_metadata(config.sample_metadata, config.selected_conditions)
    sources = resolve_bam_sources(samples, config.bam_root)
    backend = _backend_name(config.matrix_input.backend)
    counters: dict[str, Counter[Junction]] = {}
    source_parameters: list[dict[str, object]] = []
    source_rows = list(sources.itertuples(index=False))

    def run_sample(row: object) -> tuple[str, Counter[Junction], dict[str, object]]:
        sample_id = str(row.sample_id)
        counts, parameters = count_sample_junctions(
            bam=Path(row.bam_path),
            bai=Path(row.bai_path),
            settings=config.matrix_input,
            backend=backend,
            cache_dir=cache_dir,
            sample_id=sample_id,
            progress=progress,
        )
        return sample_id, counts, parameters

    max_workers = min(config.matrix_input.max_parallel_samples, len(source_rows))
    if max_workers == 1:
        completed_samples = [run_sample(row) for row in source_rows]
    else:
        _emit(progress, f"Scanning BAMs with {max_workers} parallel sample workers")
        completed_samples = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(run_sample, row): row for row in source_rows}
            for future in as_completed(futures):
                completed_samples.append(future.result())
    for sample_id, counts, parameters in completed_samples:
        counters[sample_id] = counts
        source_parameters.append(parameters)

    sample_order = samples["sample_id"].tolist()
    count_matrix = pd.DataFrame(
        {
            sample_id: pd.Series(counter, dtype="float64")
            for sample_id, counter in counters.items()
        }
    )
    count_matrix = count_matrix.reindex(columns=sample_order).fillna(0).astype(np.int64)
    count_matrix.index = pd.MultiIndex.from_tuples(
        count_matrix.index.to_list(),
        names=["chrom", "intron_start_0based", "intron_end_0based"],
    )
    count_matrix = count_matrix.sort_index()
    count_matrix, cluster_map = _filter_and_cluster(count_matrix, config.matrix_input)
    if count_matrix.empty:
        raise ValueError(
            "No competing BAM junction clusters passed matrix_input filters; "
            "lower min_total_junction_reads/min_samples_with_reads or inspect BAM inputs"
        )

    junction_index = count_matrix.index.to_list()
    cluster_ids = pd.Series(
        [cluster_map[tuple(index)] for index in junction_index],
        index=count_matrix.index,
        dtype="string",
    )
    cluster_totals = count_matrix.groupby(cluster_ids).transform("sum")
    psi = count_matrix.divide(cluster_totals.where(cluster_totals.gt(0)))
    sjc = cluster_totals - count_matrix
    event_keys = [
        f"JUNC|{chrom}|{start}|{end}"
        for chrom, start, end in junction_index
    ]
    for frame in (psi, count_matrix, sjc):
        frame.index = event_keys
        frame.index.name = "event_key"

    retained_junctions = {
        (str(chrom), int(start), int(end))
        for chrom, start, end in junction_index
    }
    _emit(progress, f"Annotating {len(retained_junctions):,} retained junctions with GTF")
    annotation = annotate_junctions_with_gtf(retained_junctions, config.gtf)
    annotation["event_key"] = annotation.apply(
        lambda row: (
            f"JUNC|{row['chrom']}|{int(row['intron_start_0based'])}|"
            f"{int(row['intron_end_0based'])}"
        ),
        axis=1,
    )
    annotation["cluster_id"] = [
        cluster_map[(str(row.chrom), int(row.intron_start_0based), int(row.intron_end_0based))]
        for row in annotation.itertuples(index=False)
    ]
    annotation["event_type"] = "JUNCTION"
    annotation["event_coordinates"] = annotation.apply(
        lambda row: (
            f"{row['chrom']}:{int(row['intron_start_0based'])}-"
            f"{int(row['intron_end_0based'])}"
        ),
        axis=1,
    )

    rmats_root = resolve_rmats_root(config.rmats_root)
    _emit(progress, "Mapping retained BAM junctions to rMATS event annotations")
    rmats_mapping = map_rmats_to_junctions(retained_junctions, rmats_root, config.rmats_mode)
    if not rmats_mapping.empty:
        aggregate = (
            rmats_mapping.groupby("event_key", observed=True)
            .agg(
                rmats_event_types=("rmats_event_type", lambda values: ";".join(sorted(set(values)))),
                rmats_event_count=("rmats_event_id", "nunique"),
                rmats_comparison_count=("comparison_id", "nunique"),
                rmats_gene_ids=(
                    "rmats_gene_id",
                    lambda values: ";".join(
                        sorted({str(value) for value in values if pd.notna(value) and str(value)})
                    ),
                ),
                rmats_gene_symbols=(
                    "rmats_gene_symbol",
                    lambda values: ";".join(
                        sorted({str(value) for value in values if pd.notna(value) and str(value)})
                    ),
                ),
                rmats_junction_roles=(
                    "junction_role",
                    lambda values: ";".join(sorted(set(values))),
                ),
                best_rMATS_FDR=("rMATS_FDR", "min"),
                max_abs_rMATS_dPSI=("raw_dPSI", lambda values: pd.to_numeric(values, errors="coerce").abs().max()),
            )
            .reset_index()
        )
        annotation = annotation.merge(aggregate, on="event_key", how="left")
        annotation["event_type"] = annotation["rmats_event_types"].fillna("JUNCTION")
        annotation["gene_id"] = annotation["gene_id"].replace("", np.nan).fillna(
            annotation["rmats_gene_ids"]
        )
        annotation["gene_symbol"] = annotation["gene_symbol"].replace("", np.nan).fillna(
            annotation["rmats_gene_symbols"]
        )
    else:
        annotation["rmats_event_types"] = ""
        annotation["rmats_event_count"] = 0
        annotation["rmats_comparison_count"] = 0
        annotation["rmats_gene_ids"] = ""
        annotation["rmats_gene_symbols"] = ""
        annotation["rmats_junction_roles"] = ""
        annotation["best_rMATS_FDR"] = np.nan
        annotation["max_abs_rMATS_dPSI"] = np.nan
    annotation["source_kind"] = "bam_junction_cluster"
    annotation["usage_definition"] = "junction_reads / competing_cluster_reads"
    annotation = annotation.set_index("event_key").reindex(event_keys).reset_index()

    raw_counts = count_matrix.copy()
    raw_counts.index = event_keys
    raw_counts.index.name = "event_key"
    clusters = annotation[
        [
            "event_key",
            "cluster_id",
            "chrom",
            "intron_start_0based",
            "intron_end_0based",
            "strand",
            "gene_id",
            "gene_symbol",
            "annotated_splice_boundaries",
            "annotated_direct_junction",
        ]
    ].copy()
    sources = sources.copy()
    parameters_by_sample = {str(row["sample_id"]): row for row in source_parameters}
    sources["backend"] = backend
    sources["cache_status"] = sources["sample_id"].map(
        lambda sample: parameters_by_sample[str(sample)]["cache_status"]
    )
    sources["alignments_seen"] = sources["sample_id"].map(
        lambda sample: parameters_by_sample[str(sample)].get("alignments_seen", np.nan)
    )
    sources["unique_junctions"] = sources["sample_id"].map(
        lambda sample: parameters_by_sample[str(sample)].get("unique_junctions", np.nan)
    )
    sources["source_kind"] = "bam"
    matrix = MatrixResult(
        psi=psi,
        ijc=raw_counts,
        sjc=sjc,
        events=annotation,
        samples=samples,
        sources=sources,
    )
    return BamMatrixResult(
        matrix=matrix,
        junction_counts=raw_counts,
        clusters=clusters,
        rmats_mapping=rmats_mapping,
    )
