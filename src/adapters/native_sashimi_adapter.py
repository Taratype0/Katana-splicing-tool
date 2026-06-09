from __future__ import annotations

from dataclasses import dataclass
import gzip
import math
from pathlib import Path
import re
import shutil
import subprocess
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class ExonBlock:
    start: int
    end: int
    label: str
    role: str


@dataclass(frozen=True)
class JunctionDef:
    start: int
    end: int
    label: str
    role: str


@dataclass
class EventGeometry:
    chrom: str
    strand: str
    event_type: str
    exons: list[ExonBlock]
    junctions: list[JunctionDef]
    region_start: int
    region_end: int


class NativeSashimiAdapter:
    """Dependency-light Katana-native sashimi-like plotter.

    The adapter consumes the selected-event manifest already produced by Katana.
    It intentionally runs event-by-event so a failure never aborts the batch.
    """

    def __init__(self, samtools_path: str = "samtools") -> None:
        self.samtools_path = samtools_path

    def run_manifest(self, manifest: pd.DataFrame, output_root: Path, gtf_path: Path | None = None) -> tuple[str, pd.DataFrame]:
        logs: list[str] = []
        failures: list[dict[str, object]] = []
        for _, row in manifest.iterrows():
            try:
                message = self._run_one(row, output_root, gtf_path)
                logs.append(message)
            except Exception as exc:
                outdir = self._output_dir(row, output_root)
                failures.append(
                    {
                        "gene": str(row.get("geneSymbol", row.get("gene_symbol", ""))),
                        "event_type": str(row.get("event_type", "")),
                        "event_id": str(row.get("event_uid", row.get("event_id", ""))),
                        "comparison_id": str(row.get("comparison_id", "")),
                        "command": "Katana native sashimi",
                        "script_path": "src.adapters.native_sashimi_adapter.NativeSashimiAdapter",
                        "output_folder": str(outdir),
                        "return_code": "",
                        "stdout": "",
                        "stderr": "",
                        "error_message": str(exc),
                    }
                )
        return "\n".join(line for line in logs if line), pd.DataFrame(failures)

    def _run_one(self, row: pd.Series, output_root: Path, gtf_path: Path | None) -> str:
        samtools = self._samtools_executable()
        if samtools is None and not self._pysam_available() and not self._bamnostic_available():
            raise RuntimeError(
                "No BAM reader is available. Configure a samtools executable in Settings, install pysam, "
                "or install bamnostic in the Katana conda environment."
            )

        event_file = Path(str(row.get("event_file", "")))
        if not event_file.exists():
            raise RuntimeError(f"Missing event file: {event_file}")
        event_frame = pd.read_csv(event_file, sep="\t", low_memory=False)
        if event_frame.empty:
            raise RuntimeError(f"Event file is empty: {event_file}")
        event_row = event_frame.iloc[0]

        event_type = str(row.get("event_type", "") or event_row.get("event_type", "")).strip()
        if event_type not in {"SE", "RI", "A3SS", "A5SS", "MXE"}:
            raise RuntimeError(f"Unsupported native sashimi event type: {event_type or 'missing'}")

        geometry = self._event_geometry(event_row, event_type)
        gene_symbol = str(row.get("geneSymbol", row.get("gene_symbol", event_row.get("geneSymbol", ""))) or "").strip()
        gene_id = str(row.get("GeneID", row.get("gene_id", event_row.get("GeneID", ""))) or "").strip()
        if gtf_path is not None and gtf_path.exists():
            geometry = self._augment_geometry_with_gtf(geometry, gtf_path, gene_symbol, gene_id)

        b1_paths = self._read_bam_list(Path(str(row.get("b1_txt", ""))))
        b2_paths = self._read_bam_list(Path(str(row.get("b2_txt", ""))))
        if not b1_paths or not b2_paths:
            raise RuntimeError("Missing BAM inputs for one or both sashimi groups.")

        for bam in b1_paths + b2_paths:
            if not bam.exists():
                raise RuntimeError(f"Missing BAM: {bam}")
            if self._bam_index_path(bam) is None:
                raise RuntimeError(f"Missing BAM index for: {bam}")

        readable_mapper = _CompressedMapper.from_exons(geometry.exons, geometry.region_start, geometry.region_end)
        strict_mapper = _LinearMapper(geometry.region_start, geometry.region_end)
        label1 = str(row.get("label1", "group1"))
        label2 = str(row.get("label2", "group2"))
        dpsi = _to_float(row.get("dPSI", event_row.get("IncLevelDifference", np.nan)))
        fdr = _to_float(row.get("FDR", event_row.get("FDR", np.nan)))
        flipped = str(row.get("splicing_direction_flipped", "")).lower() in {"1", "true", "yes"}
        psi1 = _parse_psi_list(event_row.get("IncLevel1", ""))
        psi2 = _parse_psi_list(event_row.get("IncLevel2", ""))
        ijc1 = _parse_number_list(event_row.get("IJC_SAMPLE_1", ""))
        ijc2 = _parse_number_list(event_row.get("IJC_SAMPLE_2", ""))
        sjc1 = _parse_number_list(event_row.get("SJC_SAMPLE_1", ""))
        sjc2 = _parse_number_list(event_row.get("SJC_SAMPLE_2", ""))
        if flipped:
            psi1, psi2 = psi2, psi1
            ijc1, ijc2 = ijc2, ijc1
            sjc1, sjc2 = sjc2, sjc1

        strict_tracks_1 = [self._sample_track(bam, geometry, strict_mapper, samtools) for bam in b1_paths]
        strict_tracks_2 = [self._sample_track(bam, geometry, strict_mapper, samtools) for bam in b2_paths]
        tracks_1 = [self._sample_track(bam, geometry, readable_mapper, samtools) for bam in b1_paths]
        tracks_2 = [self._sample_track(bam, geometry, readable_mapper, samtools) for bam in b2_paths]
        outdir = self._output_dir(row, output_root) / "native_sashimi"
        outdir.mkdir(parents=True, exist_ok=True)

        event_id = str(row.get("event_uid", row.get("event_id", "")) or "").strip()
        title = f"{gene_symbol or gene_id or 'gene'} {event_type} event"
        metadata = {
            "gene": gene_symbol or gene_id,
            "gene_id": gene_id,
            "event_type": event_type,
            "event_id": event_id,
            "comparison": str(row.get("comparison_name", row.get("comparison_id", ""))),
            "label1": label1,
            "label2": label2,
            "dPSI": dpsi,
            "FDR": fdr,
        }

        rmats_tracks = [
            _rmats_track(label1, ijc1, sjc1, psi1, 0),
            _rmats_track(label2, ijc2, sjc2, psi2, 1),
        ]
        self._write_tables(outdir, metadata, geometry, strict_tracks_1, strict_tracks_2, psi1, psi2)
        self._write_rmats_counts_table(outdir, metadata, rmats_tracks)
        self._write_coverage_notes(outdir)
        self._plot_rmats_evidence(outdir, title, metadata, geometry, strict_mapper, rmats_tracks)
        self._plot_summary(
            outdir,
            title,
            metadata,
            geometry,
            strict_mapper,
            strict_tracks_1,
            strict_tracks_2,
            psi1,
            psi2,
            view_label="genomic coordinate view",
        )
        self._plot_replicates(
            outdir,
            title,
            metadata,
            geometry,
            strict_mapper,
            strict_tracks_1,
            strict_tracks_2,
            psi1,
            psi2,
            view_label="genomic coordinate view",
        )
        self._plot_summary(
            outdir,
            title,
            metadata,
            geometry,
            readable_mapper,
            tracks_1,
            tracks_2,
            psi1,
            psi2,
            output_name="native_summary_readable",
            view_label="compressed event view",
        )
        self._plot_replicates(
            outdir,
            title,
            metadata,
            geometry,
            readable_mapper,
            tracks_1,
            tracks_2,
            psi1,
            psi2,
            output_name="native_replicates_readable",
            view_label="compressed event view",
        )

        return f"Katana native sashimi finished: {metadata['gene']} {event_type} {event_id} -> {outdir}"

    def _samtools_executable(self) -> str | None:
        if Path(self.samtools_path).exists():
            return str(Path(self.samtools_path))
        return shutil.which(self.samtools_path)

    @staticmethod
    def _pysam_available() -> bool:
        try:
            import pysam  # noqa: F401
        except Exception:
            return False
        return True

    @staticmethod
    def _bamnostic_available() -> bool:
        try:
            import bamnostic  # noqa: F401
        except Exception:
            return False
        return True

    @staticmethod
    def _bam_index_path(bam: Path) -> Path | None:
        candidates = [Path(str(bam) + ".bai"), bam.with_suffix(".bai")]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    def _sample_track(self, bam: Path, geometry: EventGeometry, mapper: "_CompressedMapper", samtools: str | None) -> dict[str, object]:
        if samtools is not None:
            return self._sample_track_samtools(bam, geometry, mapper, samtools)
        if self._pysam_available():
            return self._sample_track_pysam(bam, geometry, mapper)
        return self._sample_track_bamnostic(bam, geometry, mapper)

    def _sample_track_samtools(self, bam: Path, geometry: EventGeometry, mapper: "_CompressedMapper", samtools: str) -> dict[str, object]:
        region = f"{geometry.chrom}:{geometry.region_start}-{geometry.region_end}"
        view = subprocess.run(
            [samtools, "view", "-F", "260", str(bam), region],
            capture_output=True,
            text=True,
            check=True,
        )
        mapped_reads = self._mapped_reads(bam, samtools)
        bins = np.zeros(240, dtype=float)
        junction_counts = {junction.label: 0 for junction in geometry.junctions}
        target_lookup = {(junction.start, junction.end): junction.label for junction in geometry.junctions}
        aligned_bases = 0
        aligned_read_count = 0

        for line in view.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 6 or parts[5] == "*":
                continue
            ref_pos = int(parts[3])
            ops = CIGAR_RE.findall(parts[5])
            read_aligned_bases = 0
            for idx, (length_text, op) in enumerate(ops):
                length = int(length_text)
                if op in {"M", "=", "X"}:
                    seg_start = ref_pos
                    seg_end = ref_pos + length - 1
                    self._add_coverage_segment(bins, mapper.map(seg_start), mapper.map(seg_end))
                    read_aligned_bases += length
                    ref_pos += length
                elif op == "D":
                    ref_pos += length
                elif op == "N":
                    intron_start = ref_pos
                    intron_end = ref_pos + length - 1
                    label = target_lookup.get((intron_start, intron_end))
                    if label is not None:
                        left_anchor, right_anchor = _junction_anchors(ops, idx)
                        if min(left_anchor, right_anchor) >= 8:
                            junction_counts[label] = junction_counts.get(label, 0) + 1
                    ref_pos += length
                elif op in {"I", "S", "H", "P"}:
                    continue
            if read_aligned_bases > 0:
                aligned_bases += read_aligned_bases
                aligned_read_count += 1

        coverage, mean_aligned_read_length = self._normalize_coverage_bins(bins, mapped_reads, aligned_bases, aligned_read_count)
        return {
            "sample": bam.stem,
            "bam": str(bam),
            "coverage": coverage,
            "junction_counts": junction_counts,
            "mapped_reads": mapped_reads,
            "mean_aligned_read_length": mean_aligned_read_length,
            "coverage_normalization": "RPKM-style: average read-depth per bin / mapped-reads-per-million * 1000 / mean aligned read length",
        }

    def _sample_track_pysam(self, bam: Path, geometry: EventGeometry, mapper: "_CompressedMapper") -> dict[str, object]:
        import pysam

        bins = np.zeros(240, dtype=float)
        junction_counts = {junction.label: 0 for junction in geometry.junctions}
        target_lookup = {(junction.start, junction.end): junction.label for junction in geometry.junctions}
        aligned_bases = 0
        aligned_read_count = 0
        with pysam.AlignmentFile(str(bam), "rb") as alignment:
            mapped_reads = int(getattr(alignment, "mapped", 0) or 0)
            if mapped_reads <= 0:
                try:
                    mapped_reads = int(sum(item.mapped for item in alignment.get_index_statistics()))
                except Exception:
                    mapped_reads = 0
            for record in alignment.fetch(geometry.chrom, max(0, geometry.region_start - 1), geometry.region_end):
                if record.is_unmapped or record.is_secondary:
                    continue
                ref_pos = int(record.reference_start) + 1
                cigar = list(record.cigartuples or [])
                read_aligned_bases = 0
                for idx, item in enumerate(cigar):
                    op, length = int(item[0]), int(item[1])
                    if op in {0, 7, 8}:  # M, =, X
                        seg_start = ref_pos
                        seg_end = ref_pos + length - 1
                        self._add_coverage_segment(bins, mapper.map(seg_start), mapper.map(seg_end))
                        read_aligned_bases += length
                        ref_pos += length
                    elif op == 2:  # D
                        ref_pos += length
                    elif op == 3:  # N
                        intron_start = ref_pos
                        intron_end = ref_pos + length - 1
                        label = target_lookup.get((intron_start, intron_end))
                        if label is not None:
                            left_anchor, right_anchor = _junction_anchors_bamnostic(cigar, idx)
                            if min(left_anchor, right_anchor) >= 8:
                                junction_counts[label] = junction_counts.get(label, 0) + 1
                        ref_pos += length
                    elif op in {1, 4, 5, 6}:  # I, S, H, P
                        continue
                if read_aligned_bases > 0:
                    aligned_bases += read_aligned_bases
                    aligned_read_count += 1

        coverage, mean_aligned_read_length = self._normalize_coverage_bins(bins, mapped_reads, aligned_bases, aligned_read_count)
        return {
            "sample": bam.stem,
            "bam": str(bam),
            "coverage": coverage,
            "junction_counts": junction_counts,
            "mapped_reads": mapped_reads,
            "mean_aligned_read_length": mean_aligned_read_length,
            "coverage_normalization": "RPKM-style: average read-depth per bin / mapped-reads-per-million * 1000 / mean aligned read length",
        }

    def _sample_track_bamnostic(self, bam: Path, geometry: EventGeometry, mapper: "_CompressedMapper") -> dict[str, object]:
        import bamnostic

        index_path = self._bam_index_path(bam)
        if index_path is None:
            raise RuntimeError(f"Missing BAM index for: {bam}")
        alignment = bamnostic.AlignmentFile(str(bam), "rb", filepath_index=str(index_path))
        bins = np.zeros(240, dtype=float)
        junction_counts = {junction.label: 0 for junction in geometry.junctions}
        target_lookup = {(junction.start, junction.end): junction.label for junction in geometry.junctions}
        aligned_bases = 0
        aligned_read_count = 0

        for record in alignment.fetch(contig=geometry.chrom, start=max(0, geometry.region_start - 1), stop=geometry.region_end):
            if getattr(record, "is_unmapped", False) or getattr(record, "is_secondary", False):
                continue
            ref_pos = int(record.reference_start) + 1
            cigar = list(getattr(record, "cigartuples", None) or getattr(record, "cigar", []) or [])
            read_aligned_bases = 0
            for idx, item in enumerate(cigar):
                op, length = _normalize_bamnostic_cigar_item(item)
                if op in {0, 7, 8}:  # M, =, X
                    seg_start = ref_pos
                    seg_end = ref_pos + length - 1
                    self._add_coverage_segment(bins, mapper.map(seg_start), mapper.map(seg_end))
                    read_aligned_bases += length
                    ref_pos += length
                elif op == 2:  # D
                    ref_pos += length
                elif op == 3:  # N
                    intron_start = ref_pos
                    intron_end = ref_pos + length - 1
                    label = target_lookup.get((intron_start, intron_end))
                    if label is not None:
                        left_anchor, right_anchor = _junction_anchors_bamnostic(cigar, idx)
                        if min(left_anchor, right_anchor) >= 8:
                            junction_counts[label] = junction_counts.get(label, 0) + 1
                    ref_pos += length
                elif op in {1, 4, 5, 6}:  # I, S, H, P
                    continue
            if read_aligned_bases > 0:
                aligned_bases += read_aligned_bases
                aligned_read_count += 1

        mapped_reads = int(getattr(alignment, "mapped", 0) or 0)
        coverage, mean_aligned_read_length = self._normalize_coverage_bins(bins, mapped_reads, aligned_bases, aligned_read_count)
        return {
            "sample": bam.stem,
            "bam": str(bam),
            "coverage": coverage,
            "junction_counts": junction_counts,
            "mapped_reads": mapped_reads,
            "mean_aligned_read_length": mean_aligned_read_length,
            "coverage_normalization": "RPKM-style: average read-depth per bin / mapped-reads-per-million * 1000 / mean aligned read length",
        }

    @staticmethod
    def _normalize_coverage_bins(
        bins: np.ndarray,
        mapped_reads: int,
        aligned_bases: int,
        aligned_read_count: int,
    ) -> tuple[np.ndarray, float]:
        denominator = mapped_reads / 1_000_000.0 if mapped_reads > 0 else 1.0
        mean_aligned_read_length = aligned_bases / aligned_read_count if aligned_read_count > 0 else np.nan
        rpkm_scale = 1000.0 / mean_aligned_read_length if mean_aligned_read_length and not np.isnan(mean_aligned_read_length) else 1.0
        return bins / denominator * rpkm_scale, mean_aligned_read_length

    def _mapped_reads(self, bam: Path, samtools: str) -> int:
        try:
            result = subprocess.run([samtools, "idxstats", str(bam)], capture_output=True, text=True, check=True)
        except Exception:
            return 0
        total = 0
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                try:
                    total += int(parts[2])
                except ValueError:
                    continue
        return total

    @staticmethod
    def _add_coverage_segment(bins: np.ndarray, x1: float, x2: float) -> None:
        lo, hi = sorted((max(0.0, x1), min(100.0, x2)))
        if hi <= lo:
            index = int(round(lo / 100.0 * (len(bins) - 1)))
            bins[max(0, min(len(bins) - 1, index))] += 1.0
            return
        edges = np.linspace(0.0, 100.0, len(bins) + 1)
        bin_width = edges[1] - edges[0]
        start = max(0, min(len(bins) - 1, int(np.floor(lo / bin_width))))
        end = max(0, min(len(bins) - 1, int(np.floor((hi - 1e-9) / bin_width))))
        for index in range(start, end + 1):
            overlap = max(0.0, min(hi, edges[index + 1]) - max(lo, edges[index]))
            if overlap > 0:
                bins[index] += overlap / bin_width

    def _write_tables(
        self,
        outdir: Path,
        metadata: dict[str, object],
        geometry: EventGeometry,
        tracks_1: list[dict[str, object]],
        tracks_2: list[dict[str, object]],
        psi1: list[float],
        psi2: list[float],
    ) -> None:
        rows: list[dict[str, object]] = []
        for group, tracks, psis in ((metadata["label1"], tracks_1, psi1), (metadata["label2"], tracks_2, psi2)):
            for idx, track in enumerate(tracks):
                counts = track["junction_counts"]
                for junction in geometry.junctions:
                    rows.append(
                        {
                            **metadata,
                            "group": group,
                            "sample": track["sample"],
                            "psi": psis[idx] if idx < len(psis) else np.nan,
                            "junction_label": junction.label,
                            "junction_role": junction.role,
                            "junction_start": junction.start,
                            "junction_end": junction.end,
                            "junction_count_anchor8": counts.get(junction.label, 0),
                            "mapped_reads": track["mapped_reads"],
                            "mean_aligned_read_length": track.get("mean_aligned_read_length", np.nan),
                            "coverage_normalization": track.get("coverage_normalization", ""),
                        }
                    )
        pd.DataFrame(rows).to_csv(outdir / "native_sashimi_junction_counts.tsv", sep="\t", index=False)
        pd.DataFrame(
            [
                {"group": metadata["label1"], "mean_PSI": np.nanmean(psi1) if psi1 else np.nan, "replicate_PSI": ",".join(f"{x:.4g}" for x in psi1)},
                {"group": metadata["label2"], "mean_PSI": np.nanmean(psi2) if psi2 else np.nan, "replicate_PSI": ",".join(f"{x:.4g}" for x in psi2)},
            ]
        ).to_csv(outdir / "native_sashimi_psi_summary.tsv", sep="\t", index=False)

    @staticmethod
    def _write_rmats_counts_table(outdir: Path, metadata: dict[str, object], rmats_tracks: list[dict[str, object]]) -> None:
        rows: list[dict[str, object]] = []
        for track in rmats_tracks:
            rows.append(
                {
                    **metadata,
                    "group": track["label"],
                    "mean_IJC": track["mean_ijc"],
                    "mean_SJC": track["mean_sjc"],
                    "replicate_IJC": ",".join(_fmt_float(value) for value in track["ijc"]),
                    "replicate_SJC": ",".join(_fmt_float(value) for value in track["sjc"]),
                    "mean_PSI": track["psi"],
                    "replicate_PSI": ",".join(_fmt_float(value) for value in track["replicate_psis"]),
                    "source": "rMATS IJC/SJC/IncLevel columns",
                }
            )
        pd.DataFrame(rows).to_csv(outdir / "native_sashimi_rmats_event_counts.tsv", sep="\t", index=False)

    @staticmethod
    def _write_coverage_notes(outdir: Path) -> None:
        (outdir / "native_sashimi_coverage_notes.txt").write_text(
            "\n".join(
                [
                    "Katana native sashimi coverage notes",
                    "====================================",
                    "",
                    "Main sashimi signal: coverage is counted from aligned read M/=/X CIGAR segments and normalized as an RPKM-style sashimi coverage track.",
                    "Normalization: average read-depth per plotted bin / mapped-reads-per-million * 1000 / mean aligned read length.",
                    "This is designed to replace legacy sashimiPlot/MISO-style read density panels while keeping rMATS event metadata.",
                    "Y-axis label: RPKM.",
                    "native_summary/native_replicates are the primary sashimiPlot-replacement outputs and use linear genomic coordinates.",
                    "native_summary_readable/native_replicates_readable are optional compressed event-view companions for very long loci.",
                    "Junction counts are split-read counts with >=8 bp anchors on both sides.",
                    "rMATS PSI/dPSI/FDR and IJC/SJC values are inherited from the selected rMATS event and Katana comparison direction.",
                    "native_rmats_evidence summarizes rMATS event-level IJC/SJC/PSI values as an annotation/QC companion, not as a replacement for BAM coverage.",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    def _plot_rmats_evidence(
        self,
        outdir: Path,
        title: str,
        metadata: dict[str, object],
        geometry: EventGeometry,
        mapper: "_CompressedMapper",
        rmats_tracks: list[dict[str, object]],
    ) -> None:
        fig = plt.figure(figsize=(14.5, 8.05))
        gs = fig.add_gridspec(nrows=6, ncols=1, height_ratios=[0.85, 1.15, 1.15, 0.18, 0.72, 0.86], hspace=0.18)
        ax_title = fig.add_subplot(gs[0])
        ax_title.axis("off")
        ax_title.text(0.02, 0.72, title, ha="left", va="center", fontsize=21, fontweight="bold", color="#111827")
        ax_title.text(
            0.02,
            0.38,
            "rMATS event evidence: IJC/SJC counts, replicate PSI, dPSI and FDR",
            ha="left",
            va="center",
            fontsize=12.5,
            color="#4B5563",
        )
        ax_title.text(0.78, 0.72, f"dPSI: {_fmt_signed(metadata.get('dPSI'))}", ha="left", va="center", fontsize=18, fontweight="bold", color="#006B4F")
        ax_title.text(0.78, 0.42, f"FDR: {_fmt_float(metadata.get('FDR'))}", ha="left", va="center", fontsize=16.5, fontweight="bold", color="#111827")

        max_count = max(
            [float(track.get("mean_ijc", 0) or 0) for track in rmats_tracks]
            + [float(track.get("mean_sjc", 0) or 0) for track in rmats_tracks]
            + [1.0]
        )
        for row_index, track in enumerate(rmats_tracks):
            ax = fig.add_subplot(gs[row_index + 1])
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 1)
            ax.axis("off")
            group_index = int(track.get("group_index", row_index))
            line_color = "#4B5563" if group_index == 0 else "#2563EB"
            ax.text(-0.045, 0.52, str(track.get("label", "")), transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=9.4, fontweight="bold", clip_on=False)
            ax.hlines(0.18, 3, 97, color="#111827", linewidth=0.9)
            for exon in geometry.exons:
                color = "#F8C9CF" if exon.role in {"event", "alternative", "retained", "skipped"} else "#DDEBFA"
                self._draw_exon_box(ax, mapper, exon, 0.18, color, show_label=row_index == 0)
            self._draw_rmats_arcs(ax, geometry, mapper, track, max_count)
            psi = _to_float(track.get("psi", np.nan))
            ax.text(
                1.055,
                0.72,
                f"PSI: {_fmt_float(psi)}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=12.5,
                fontweight="bold",
                color=line_color,
                clip_on=False,
            )
            ax.text(
                1.055,
                0.47,
                f"IJC: {_fmt_float(track.get('mean_ijc'))} | SJC: {_fmt_float(track.get('mean_sjc'))}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=10.5,
                color="#374151",
                clip_on=False,
            )
            _draw_psi_strip(ax, psi, track.get("replicate_psis", []), line_color)

        ax_gene = fig.add_subplot(gs[-2])
        structure_mapper = _CompressedMapper.from_exons(geometry.exons, geometry.region_start, geometry.region_end)
        self._draw_miso_transcript_track(ax_gene, geometry, mapper)
        ax_model = fig.add_subplot(gs[-1])
        self._draw_event_schematic(ax_model, geometry, structure_mapper)
        fig.legend(
            handles=[
                Patch(facecolor="#2457D6", label="rMATS inclusion junction mean count"),
                Patch(facecolor="#E41F26", label="rMATS skipping/alternative junction mean count"),
                Patch(facecolor="#DDEBFA", edgecolor="#111827", label="constitutive exon"),
                Patch(facecolor="#F8C9CF", edgecolor="#111827", label="event exon/region"),
            ],
            loc="lower center",
            bbox_to_anchor=(0.5, 0.010),
            ncol=4,
            frameon=False,
            fontsize=9.8,
        )
        fig.subplots_adjust(left=0.105, right=0.860, top=0.95, bottom=0.145)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(outdir / f"native_rmats_evidence.{suffix}", bbox_inches="tight", facecolor="white", dpi=240 if suffix == "png" else None)
        plt.close(fig)

    def _draw_rmats_arcs(
        self,
        ax,
        geometry: EventGeometry,
        mapper: "_CompressedMapper",
        track: dict[str, object],
        max_count: float,
    ) -> None:
        mean_ijc = max(0.0, float(track.get("mean_ijc", 0) or 0))
        mean_sjc = max(0.0, float(track.get("mean_sjc", 0) or 0))
        label_positions: list[tuple[float, float]] = []
        for index, junction in enumerate(geometry.junctions):
            count = mean_ijc if junction.role == "inclusion" else mean_sjc
            if count <= 0:
                continue
            x1 = mapper.map(junction.start - 1)
            x2 = mapper.map(junction.end + 1)
            span = max(2.0, abs(x2 - x1))
            role = "inclusion" if junction.role == "inclusion" else "alternative"
            color = "#2457D6" if role == "inclusion" else "#E41F26"
            lane = index if role == "inclusion" else index + 1
            y0 = 0.30
            apex = min(0.92, 0.52 + 0.10 * lane + 0.08 * min(1.0, span / 60.0))
            frac = max(0.0, min(1.0, count / max(max_count, 1.0)))
            # Keep count differences visible without letting high-count arcs dominate the plot.
            linewidth = 0.55 + 3.9 * (frac ** 1.12)
            pull = min(span * 0.12, 7.0)
            path = MplPath(
                [
                    (x1, y0),
                    (x1 + (pull if x2 >= x1 else -pull), apex),
                    (x2 - (pull if x2 >= x1 else -pull), apex),
                    (x2, y0),
                ],
                [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
            )
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=linewidth, capstyle="round", alpha=0.92))
            label_x = (x1 + x2) / 2.0
            label_y = apex + 0.045
            while any(abs(label_x - old_x) < 5.5 and abs(label_y - old_y) < 0.06 for old_x, old_y in label_positions):
                label_y += 0.065
            label_positions.append((label_x, label_y))
            ax.text(
                label_x,
                label_y,
                _fmt_float(count),
                color=color,
                ha="center",
                va="bottom",
                fontsize=9.5,
                fontweight="bold",
                clip_on=False,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.4},
            )

    def _plot_summary(
        self,
        outdir: Path,
        title: str,
        metadata: dict[str, object],
        geometry: EventGeometry,
        mapper: "_CompressedMapper",
        tracks_1: list[dict[str, object]],
        tracks_2: list[dict[str, object]],
        psi1: list[float],
        psi2: list[float],
        *,
        output_name: str = "native_summary",
        view_label: str = "compressed event view",
    ) -> None:
        group_tracks = [
            _mean_track(metadata["label1"], tracks_1, psi1),
            _mean_track(metadata["label2"], tracks_2, psi2),
        ]
        group_tracks[0]["group_index"] = 0
        group_tracks[1]["group_index"] = 1
        self._plot_tracks(outdir / output_name, title, metadata, geometry, mapper, group_tracks, view_label=view_label)

    def _plot_replicates(
        self,
        outdir: Path,
        title: str,
        metadata: dict[str, object],
        geometry: EventGeometry,
        mapper: "_CompressedMapper",
        tracks_1: list[dict[str, object]],
        tracks_2: list[dict[str, object]],
        psi1: list[float],
        psi2: list[float],
        *,
        output_name: str = "native_replicates",
        view_label: str = "compressed event view",
    ) -> None:
        tracks: list[dict[str, object]] = []
        for group, source, psis in ((metadata["label1"], tracks_1, psi1), (metadata["label2"], tracks_2, psi2)):
            for idx, track in enumerate(source):
                item = dict(track)
                item["label"] = f"{group} {idx + 1}"
                item["psi"] = psis[idx] if idx < len(psis) else np.nan
                item["replicate_psis"] = psis
                item["group_index"] = 0 if group == metadata["label1"] else 1
                tracks.append(item)
        self._plot_tracks(outdir / output_name, title, metadata, geometry, mapper, tracks, view_label=view_label)

    def _plot_tracks(
        self,
        output_stem: Path,
        title: str,
        metadata: dict[str, object],
        geometry: EventGeometry,
        mapper: "_CompressedMapper",
        tracks: list[dict[str, object]],
        *,
        view_label: str,
    ) -> None:
        n_tracks = max(1, len(tracks))
        fig_height = max(8.25, 3.25 + n_tracks * 1.35)
        fig = plt.figure(figsize=(14.5, fig_height))
        gs = fig.add_gridspec(
            nrows=n_tracks + 4,
            ncols=1,
            height_ratios=[0.85] + [1.0] * n_tracks + [0.24, 0.72, 0.86],
            hspace=0.18,
        )
        ax_title = fig.add_subplot(gs[0])
        ax_title.axis("off")
        ax_title.text(0.02, 0.70, title, ha="left", va="center", fontsize=21, fontweight="bold", color="#111827")
        ax_title.text(0.02, 0.38, f"event_id: {metadata.get('event_id', '')}", ha="left", va="center", fontsize=12.5, color="#4B5563")
        dpsi = metadata.get("dPSI")
        fdr = metadata.get("FDR")
        ax_title.text(0.84, 0.72, f"dPSI: {_fmt_signed(dpsi)}", ha="left", va="center", fontsize=18, fontweight="bold", color="#006B4F")
        ax_title.text(0.84, 0.43, f"FDR: {_fmt_float(fdr)}", ha="left", va="center", fontsize=16.5, fontweight="bold", color="#111827")

        coverage_max = max(float(np.nanmax(track["coverage"])) if len(track["coverage"]) else 0.0 for track in tracks) if tracks else 1.0
        coverage_ylim = max(coverage_max * 1.28, 1.0)
        x = np.linspace(0, 100, 240)
        axes = [fig.add_subplot(gs[i + 1]) for i in range(n_tracks)]
        for idx, (ax, track) in enumerate(zip(axes, tracks, strict=False)):
            group_index = int(track.get("group_index", idx))
            fill = "#C9D1DB" if group_index == 0 else "#BFD4F7"
            line = "#4B5563" if group_index == 0 else "#2563EB"
            ax.fill_between(x, 0, track["coverage"], color=fill, linewidth=0, alpha=0.95)
            ax.plot(x, track["coverage"], color=line, linewidth=1.4)
            ax.set_xlim(0, 100)
            ax.set_ylim(0, coverage_ylim)
            ax.set_ylabel("")
            ax.text(-0.050, 0.5, "RPKM", transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=9.5, clip_on=False)
            label = str(track.get("label", ""))
            label_size = 8.2 if len(label) > 20 else 9.2
            ax.text(-0.090, 0.5, label, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=label_size, fontweight="bold", clip_on=False)
            _draw_psi_strip(ax, _to_float(track.get("psi", np.nan)), track.get("replicate_psis", []), line)
            self._draw_junction_arcs(ax, geometry, mapper, track["junction_counts"], np.asarray(track["coverage"], dtype=float), coverage_ylim)
            if idx == len(axes) - 1:
                ticks = np.linspace(0, 100, 5)
                ax.set_xticks(ticks)
                ax.set_xticklabels([f"{int(mapper.inverse(value)):,}" for value in ticks], fontsize=9)
                ax.set_xlabel(f"Genomic coordinate ({geometry.chrom}, {geometry.strand or '?'} strand)", fontsize=9.8)
            else:
                ax.tick_params(axis="x", labelbottom=False)
            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)

        ax_gene = fig.add_subplot(gs[-2])
        structure_mapper = _CompressedMapper.from_exons(geometry.exons, geometry.region_start, geometry.region_end)
        self._draw_miso_transcript_track(ax_gene, geometry, mapper)
        ax_model = fig.add_subplot(gs[-1])
        self._draw_event_schematic(ax_model, geometry, structure_mapper)
        fig.legend(
            handles=[
                Patch(facecolor="#2457D6", label="inclusion junction count"),
                Patch(facecolor="#E41F26", label="alternative/skipping junction count"),
                Patch(facecolor="#DDEBFA", edgecolor="#111827", label="constitutive exon"),
                Patch(facecolor="#F8C9CF", edgecolor="#111827", label="event exon/region"),
            ],
            loc="lower center",
            bbox_to_anchor=(0.5, 0.012),
            ncol=4,
            frameon=False,
            fontsize=9.8,
        )
        fig.subplots_adjust(left=0.145, right=0.860, top=0.96, bottom=0.145)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output_stem.with_suffix(f".{suffix}"), bbox_inches="tight", facecolor="white", dpi=240 if suffix == "png" else None)
        plt.close(fig)

    def _draw_junction_arcs(self, ax, geometry: EventGeometry, mapper: "_CompressedMapper", counts: dict[str, int], coverage: np.ndarray, coverage_ylim: float) -> None:
        positive_counts = [int(value) for value in counts.values() if int(value or 0) > 0]
        max_count = max(positive_counts) if positive_counts else 1
        arc_transform = blended_transform_factory(ax.transData, ax.transAxes)
        ordered = sorted(
            list(geometry.junctions),
            key=lambda item: (0 if item.role == "inclusion" else 1, abs(mapper.map(item.end + 1) - mapper.map(item.start - 1)), item.start),
        )
        occupied_by_role: dict[str, list[tuple[float, float, int]]] = {"inclusion": [], "alternative": []}
        label_positions: list[tuple[float, float]] = []
        for junction in ordered:
            count = int(counts.get(junction.label, 0) or 0)
            if count <= 0:
                continue
            x1 = mapper.map(junction.start - 1)
            x2 = mapper.map(junction.end + 1)
            left = min(x1, x2)
            right = max(x1, x2)
            span = max(2.0, right - left)
            role = "inclusion" if junction.role == "inclusion" else "alternative"
            lane = 0
            occupied = occupied_by_role[role]
            while any(prev_lane == lane and not (right < prev_left - 2 or left > prev_right + 2) for prev_left, prev_right, prev_lane in occupied):
                lane += 1
            occupied.append((left, right, lane))
            color = "#2457D6" if role == "inclusion" else "#E41F26"
            anchor_left = self._coverage_anchor_fraction(coverage, x1, coverage_ylim)
            anchor_right = self._coverage_anchor_fraction(coverage, x2, coverage_ylim)
            local_anchor = max(anchor_left, anchor_right)
            role_floor = 0.52 + 0.10 * lane if role == "inclusion" else 0.72 + 0.08 * lane
            apex = min(0.94, max(local_anchor + 0.24, role_floor))
            count_fraction = max(0.0, min(1.0, count / max_count))
            # Keep count differences visible without letting high-count arcs dominate the plot.
            linewidth = 0.55 + 3.9 * (count_fraction ** 1.12)
            control_pull = min(span * 0.10, 6.5)
            path = MplPath(
                [
                    (x1, anchor_left),
                    (x1 + (control_pull if x2 >= x1 else -control_pull), apex),
                    (x2 - (control_pull if x2 >= x1 else -control_pull), apex),
                    (x2, anchor_right),
                ],
                [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
            )
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=linewidth, capstyle="round", alpha=0.92, transform=arc_transform))
            label_x = (x1 + x2) / 2.0
            label_y = min(0.98, apex + 0.045)
            while any(abs(label_x - old_x) < 5.5 and abs(label_y - old_y) < 0.055 for old_x, old_y in label_positions):
                label_y = min(0.98, label_y + 0.060)
            label_positions.append((label_x, label_y))
            ax.text(
                label_x,
                label_y,
                str(count),
                color=color,
                ha="center",
                va="bottom",
                fontsize=9.5,
                fontweight="bold",
                clip_on=False,
                transform=arc_transform,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.4},
            )

    @staticmethod
    def _coverage_anchor_fraction(coverage: np.ndarray, x_position: float, coverage_ylim: float) -> float:
        if coverage.size == 0:
            return 0.035
        index = int(round(np.clip(x_position, 0, 100) / 100 * (coverage.size - 1)))
        left = max(0, index - 2)
        right = min(coverage.size, index + 3)
        local = float(np.nanmax(coverage[left:right])) if right > left else 0.0
        denominator = max(float(coverage_ylim), 1e-9)
        return min(0.72, max(local / denominator + 0.025, 0.035))

    def _draw_isoforms(self, ax, geometry: EventGeometry, mapper: "_CompressedMapper", *, view_label: str, annotate: bool = True) -> None:
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 100)
        ax.set_yticks([])
        ax.set_xlabel(f"Genomic position, {view_label}", fontsize=12)
        ticks = np.linspace(0, 100, 5)
        labels = [f"{int(mapper.inverse(value)):,}" for value in ticks]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=9)
        y_top, y_bottom = 0.72, 0.35
        for y, label in ((y_top, self._top_isoform_label(geometry.event_type)), (y_bottom, self._bottom_isoform_label(geometry.event_type))):
            ax.hlines(y, 5, 95, color="#111827", linewidth=1.0)
            ax.text(3.3, y, label, ha="right", va="center", fontsize=9.4, color="#4B5563")

        exons = sorted(geometry.exons, key=lambda exon: (exon.start, exon.end))
        for exon in exons:
            color = "#F8C9CF" if exon.role in {"event", "alternative", "retained", "skipped"} else "#DDEBFA"
            self._draw_exon_box(ax, mapper, exon, y_top, color, show_label=True)

        for exon in self._bottom_isoform_exons(geometry):
            color = "#F8C9CF" if exon.role in {"event", "alternative", "retained"} and geometry.event_type != "SE" else "#DDEBFA"
            self._draw_exon_box(ax, mapper, exon, y_bottom, color, show_label=False)

        if annotate:
            self._draw_event_model_annotations(ax, geometry, mapper, y_top, y_bottom)

        strand_text = f"strand {geometry.strand or '?'}"
        ax.annotate(strand_text, xy=(97, y_bottom), xytext=(90, 0.10), arrowprops=dict(arrowstyle="->", linewidth=1.0), ha="center", fontsize=10)
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#9CA3AF")

    def _draw_miso_transcript_track(self, ax, geometry: EventGeometry, mapper: "_CompressedMapper") -> None:
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1)
        ax.axis("off")
        y_top, y_bottom = 0.70, 0.32
        top_exons = sorted(geometry.exons, key=lambda item: (item.start, item.end))
        bottom_exons = sorted(self._bottom_isoform_exons(geometry), key=lambda item: (item.start, item.end))
        ax.hlines(y_top, 0, 100, color="#111827", linewidth=0.95)
        ax.hlines(y_bottom, 0, 100, color="#111827", linewidth=0.95)
        self._draw_fishbone_direction(ax, geometry, mapper, y_top, top_exons)
        self._draw_fishbone_direction(ax, geometry, mapper, y_bottom, bottom_exons)
        for exon in top_exons:
            self._draw_black_exon_box(ax, mapper, exon, y_top)
        for exon in bottom_exons:
            self._draw_black_exon_box(ax, mapper, exon, y_bottom)

    def _draw_fishbone_direction(self, ax, geometry: EventGeometry, mapper: "_CompressedMapper", y: float, exons: list[ExonBlock]) -> None:
        exon_intervals: list[tuple[float, float]] = []
        for exon in exons:
            x1 = mapper.map(exon.start)
            x2 = mapper.map(exon.end)
            exon_intervals.append((min(x1, x2), max(x1, x2)))
        if not exon_intervals:
            return
        exon_intervals = sorted(exon_intervals)
        segments = [(0.0, exon_intervals[0][0]), *[(exon_intervals[i][1], exon_intervals[i + 1][0]) for i in range(len(exon_intervals) - 1)], (exon_intervals[-1][1], 100.0)]
        direction = 1 if str(geometry.strand or "+") != "-" else -1
        for start, end in segments:
            if end - start < 3.0:
                continue
            spacing = 3.8 if end - start < 20 else 5.2
            positions = np.arange(start + spacing, end - spacing * 0.35, spacing)
            for x in positions:
                self._draw_direction_tick(ax, x, y, direction)

    @staticmethod
    def _draw_direction_tick(ax, x: float, y: float, direction: int) -> None:
        dx = 0.55 * direction
        ax.plot([x - dx, x, x - dx], [y + 0.035, y, y - 0.035], color="#111827", linewidth=0.55, solid_capstyle="round")

    def _draw_event_schematic(self, ax, geometry: EventGeometry, mapper: "_CompressedMapper") -> None:
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1)
        ax.axis("off")
        y_top, y_bottom = 0.74, 0.36
        ax.hlines(y_top, 9, 91, color="#111827", linewidth=0.95)
        ax.hlines(y_bottom, 9, 91, color="#111827", linewidth=0.95)
        ax.text(7.2, y_top, self._top_isoform_label(geometry.event_type), ha="right", va="center", fontsize=8.4, color="#4B5563")
        ax.text(7.2, y_bottom, self._bottom_isoform_label(geometry.event_type), ha="right", va="center", fontsize=8.4, color="#4B5563")

        for exon in sorted(geometry.exons, key=lambda item: (item.start, item.end)):
            color = "#F8C9CF" if exon.role in {"event", "alternative", "retained", "skipped"} else "#DDEBFA"
            self._draw_exon_box(ax, mapper, exon, y_top, color, show_label=True)

        for exon in self._bottom_isoform_exons(geometry):
            color = "#F8C9CF" if exon.role in {"event", "alternative", "retained"} and geometry.event_type != "SE" else "#DDEBFA"
            self._draw_exon_box(ax, mapper, exon, y_bottom, color, show_label=False)

        strand_text = f"strand {geometry.strand or '?'}"
        ax.annotate(strand_text, xy=(94, y_bottom), xytext=(86, 0.07), arrowprops=dict(arrowstyle="->", linewidth=0.85), ha="center", fontsize=8.8)

    def _draw_event_model_annotations(self, ax, geometry: EventGeometry, mapper: "_CompressedMapper", y_top: float, y_bottom: float) -> None:
        event_exons = [exon for exon in geometry.exons if exon.role in {"event", "alternative", "retained", "skipped"}]
        if not event_exons:
            return
        left = min(mapper.map(exon.start) for exon in event_exons)
        right = max(mapper.map(exon.end) for exon in event_exons)
        center = (left + right) / 2.0
        width = max(2.0, right - left)
        if geometry.event_type == "SE":
            self._mini_tag(ax, center, y_top + 0.25, "included", "#059669")
            ax.add_patch(
                Rectangle(
                    (left, y_bottom - 0.055),
                    max(1.2, width),
                    0.11,
                    facecolor="none",
                    edgecolor="#E41F26",
                    linewidth=1.3,
                    linestyle=(0, (3, 2)),
                )
            )
            self._mini_tag(ax, center, y_bottom - 0.18, "skipped", "#E41F26")
        elif geometry.event_type == "RI":
            self._mini_tag(ax, center, y_top + 0.25, "retained", "#059669")
            ax.plot([left, right], [y_bottom + 0.11, y_bottom + 0.11], color="#E41F26", linewidth=2.0, linestyle=(0, (3, 2)))
            self._mini_tag(ax, center, y_bottom + 0.19, "spliced out", "#E41F26")
        elif geometry.event_type in {"A3SS", "A5SS"}:
            self._mini_tag(ax, center, y_top + 0.25, "long site", "#2457D6")
            self._mini_tag(ax, center + max(4.0, width * 0.35), y_bottom + 0.19, "short site", "#E41F26")
        elif geometry.event_type == "MXE":
            first = [exon for exon in event_exons if exon.role == "event"]
            second = [exon for exon in event_exons if exon.role == "alternative"]
            first_center = _exon_group_center(first, mapper) if first else center
            second_center = _exon_group_center(second, mapper) if second else center
            self._mini_tag(ax, first_center, y_top + 0.25, "choice 1", "#2457D6")
            self._mini_tag(ax, second_center, y_bottom + 0.19, "choice 2", "#E41F26")

    @staticmethod
    def _mini_tag(ax, x: float, y: float, text: str, color: str) -> None:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=8.2,
            color=color,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": color, "linewidth": 0.8},
        )

    @staticmethod
    def _event_interpretation_text(geometry: EventGeometry) -> str:
        event_label = NativeSashimiAdapter._event_region_label(geometry)
        if geometry.event_type == "SE":
            return f"Event interpretation: skipped-exon event. The highlighted exon/region ({event_label}) is included in the inclusion isoform and skipped in the skipping isoform."
        if geometry.event_type == "RI":
            return f"Event interpretation: retained-intron event. The highlighted intronic region ({event_label}) is retained in one isoform and removed by splicing in the spliced isoform."
        if geometry.event_type == "A3SS":
            return f"Event interpretation: alternative 3' splice-site event. The highlighted region ({event_label}) marks the long/short 3' splice-site choice."
        if geometry.event_type == "A5SS":
            return f"Event interpretation: alternative 5' splice-site event. The highlighted region ({event_label}) marks the long/short 5' splice-site choice."
        if geometry.event_type == "MXE":
            return f"Event interpretation: mutually exclusive exon event. The highlighted exons/regions ({event_label}) are used as alternative exon choices."
        return f"Event interpretation: highlighted event region ({event_label}) defines the alternative splicing event."

    @staticmethod
    def _event_region_label(geometry: EventGeometry) -> str:
        labels: list[str] = []
        for exon in geometry.exons:
            if exon.role in {"event", "alternative", "retained", "skipped"}:
                label = _display_exon_label(exon)
                if label and label not in labels:
                    labels.append(label)
        return ", ".join(labels) if labels else "event region"

    @staticmethod
    def _draw_exon_box(ax, mapper: "_CompressedMapper", exon: ExonBlock, y: float, color: str, show_label: bool) -> None:
        x1 = mapper.map(exon.start)
        x2 = mapper.map(exon.end)
        width = max(1.2, abs(x2 - x1))
        left = min(x1, x2)
        ax.add_patch(Rectangle((left, y - 0.05), width, 0.10, facecolor=color, edgecolor="#111827", linewidth=0.95))
        if show_label:
            ax.text(left + width / 2, y + 0.13, _display_exon_label(exon), ha="center", va="center", fontsize=9.4, fontweight="bold")
            ax.text(left + width / 2, y - 0.14, _display_exon_role(exon), ha="center", va="center", fontsize=8.1, color="#4B5563")

    @staticmethod
    def _draw_black_exon_box(ax, mapper: "_CompressedMapper", exon: ExonBlock, y: float) -> None:
        x1 = mapper.map(exon.start)
        x2 = mapper.map(exon.end)
        width = max(1.2, abs(x2 - x1))
        left = min(x1, x2)
        ax.add_patch(Rectangle((left, y - 0.055), width, 0.11, facecolor="black", edgecolor="black", linewidth=0.8))

    @staticmethod
    def _top_isoform_label(event_type: str) -> str:
        return {
            "SE": "inclusion isoform",
            "RI": "retention isoform",
            "A3SS": "long 3SS isoform",
            "A5SS": "long 5SS isoform",
            "MXE": "MXE1 isoform",
        }.get(event_type, "isoform A")

    @staticmethod
    def _bottom_isoform_label(event_type: str) -> str:
        return {
            "SE": "skipping isoform",
            "RI": "spliced isoform",
            "A3SS": "short 3SS isoform",
            "A5SS": "short 5SS isoform",
            "MXE": "MXE2 isoform",
        }.get(event_type, "isoform B")

    @staticmethod
    def _bottom_isoform_exons(geometry: EventGeometry) -> list[ExonBlock]:
        if geometry.event_type == "SE":
            return [exon for exon in geometry.exons if exon.role != "skipped"]
        if geometry.event_type == "RI":
            return [exon for exon in geometry.exons if exon.role != "retained"]
        if geometry.event_type == "MXE":
            return [exon for exon in geometry.exons if exon.role != "event"] + [
                ExonBlock(exon.start, exon.end, exon.label, exon.role)
                for exon in geometry.exons
                if exon.role == "alternative"
            ]
        return [exon for exon in geometry.exons if exon.role != "event"] + [
            exon for exon in geometry.exons if exon.role == "alternative"
        ]

    def _event_geometry(self, row: pd.Series, event_type: str) -> EventGeometry:
        chrom = str(row.get("chr", "")).strip()
        strand = str(row.get("strand", "")).strip()
        if not chrom:
            raise RuntimeError("rMATS event row is missing chromosome.")
        exons: list[ExonBlock] = []
        junctions: list[JunctionDef] = []

        if event_type == "SE":
            upstream = _block(row, "upstreamES", "upstreamEE", "upstream", "constitutive")
            skipped = _block(row, "exonStart_0base", "exonEnd", "skipped", "skipped")
            downstream = _block(row, "downstreamES", "downstreamEE", "downstream", "constitutive")
            exons = [upstream, skipped, downstream]
            junctions = [
                _junction(upstream, skipped, "upstream-skipped", "inclusion"),
                _junction(skipped, downstream, "skipped-downstream", "inclusion"),
                _junction(upstream, downstream, "upstream-downstream", "alternative"),
            ]
        elif event_type == "RI":
            retained = _block(row, "riExonStart_0base", "riExonEnd", "retained", "retained")
            upstream = _block(row, "upstreamES", "upstreamEE", "upstream", "constitutive", required=False)
            downstream = _block(row, "downstreamES", "downstreamEE", "downstream", "constitutive", required=False)
            exons = [item for item in [upstream, retained, downstream] if item is not None]
            if upstream and downstream:
                junctions = [_junction(upstream, downstream, "spliced", "alternative")]
        elif event_type in {"A3SS", "A5SS"}:
            long_exon = _block(row, "longExonStart_0base", "longExonEnd", "long", "event")
            short_exon = _block(row, "shortES", "shortEE", "short", "alternative")
            flank = _block(row, "flankingES", "flankingEE", "flanking", "constitutive")
            exons = sorted([long_exon, short_exon, flank], key=lambda item: item.start)
            junctions = [
                _junction(flank, long_exon, "long junction", "alternative"),
                _junction(flank, short_exon, "short junction", "inclusion"),
            ]
        elif event_type == "MXE":
            upstream = _block(row, "upstreamES", "upstreamEE", "upstream", "constitutive")
            first = _block(row, "1stExonStart_0base", "1stExonEnd", "MXE1", "event")
            second = _block(row, "2ndExonStart_0base", "2ndExonEnd", "MXE2", "alternative")
            downstream = _block(row, "downstreamES", "downstreamEE", "downstream", "constitutive")
            exons = [upstream, first, second, downstream]
            junctions = [
                _junction(upstream, first, "upstream-MXE1", "inclusion"),
                _junction(first, downstream, "MXE1-downstream", "inclusion"),
                _junction(upstream, second, "upstream-MXE2", "alternative"),
                _junction(second, downstream, "MXE2-downstream", "alternative"),
            ]

        exons = [exon for exon in exons if exon is not None]
        if not exons:
            raise RuntimeError(f"Could not parse {event_type} event coordinates.")
        region_start = max(1, min(exon.start for exon in exons) - 250)
        region_end = max(exon.end for exon in exons) + 250
        return EventGeometry(chrom, strand, event_type, exons, junctions, region_start, region_end)

    def _augment_geometry_with_gtf(self, geometry: EventGeometry, gtf_path: Path, gene_symbol: str, gene_id: str) -> EventGeometry:
        gtf_exons = self._gtf_exons_for_region(gtf_path, geometry.chrom, geometry.region_start, geometry.region_end, gene_symbol, gene_id)
        if not gtf_exons:
            return geometry
        exons: list[ExonBlock] = []
        for exon in geometry.exons:
            matched_label = _matching_gtf_exon_label(exon.start, exon.end, gtf_exons)
            label = matched_label if matched_label and exon.role in {"event", "alternative", "retained", "skipped"} else exon.label
            exons.append(ExonBlock(exon.start, exon.end, label, exon.role))
        return EventGeometry(geometry.chrom, geometry.strand, geometry.event_type, sorted(exons, key=lambda exon: (exon.start, exon.end)), geometry.junctions, geometry.region_start, geometry.region_end)

    @staticmethod
    def _gtf_exons_for_region(gtf_path: Path, chrom: str, start: int, end: int, gene_symbol: str, gene_id: str) -> dict[str, list[tuple[int, int, str]]]:
        opener = gzip.open if gtf_path.suffix == ".gz" else open
        gene_symbol = gene_symbol.strip()
        gene_id_core = gene_id.split(".", 1)[0].strip()
        rows: dict[str, dict[tuple[int, int], set[str]]] = {}
        with opener(gtf_path, "rt", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[0] != chrom or parts[2] != "exon":
                    continue
                exon_start, exon_end = int(parts[3]), int(parts[4])
                if exon_end < start or exon_start > end:
                    continue
                attrs = _parse_gtf_attrs(parts[8])
                row_gene = str(attrs.get("gene_name", ""))
                row_gene_id = str(attrs.get("gene_id", "")).split(".", 1)[0]
                if (gene_symbol and row_gene == gene_symbol) or (gene_id_core and row_gene_id == gene_id_core):
                    label = _gtf_exon_label(attrs)
                    transcript_key = str(attrs.get("transcript_id", "") or attrs.get("transcript_name", "") or "unknown_transcript")
                    rows.setdefault(transcript_key, {}).setdefault((exon_start, exon_end), set()).add(label)
        return {
            transcript: [(start, end, _collapse_gtf_labels(labels)) for (start, end), labels in sorted(exons.items())]
            for transcript, exons in rows.items()
        }

    @staticmethod
    def _read_bam_list(path: Path) -> list[Path]:
        if not path.exists():
            return []
        values: list[Path] = []
        for raw_line in path.read_text(encoding="utf-8").splitlines():
            for token in raw_line.split(","):
                token = token.strip()
                if token:
                    values.append(Path(token))
        return values

    @staticmethod
    def _output_dir(row: pd.Series, output_root: Path) -> Path:
        outdir = Path(str(row.get("outdir", "native_sashimi")))
        if not outdir.is_absolute():
            outdir = output_root / "06_sashimi" / "plots" / outdir
        return outdir


class _CompressedMapper:
    def __init__(self, segments: list[tuple[int, int, float, float]]) -> None:
        self.segments = segments

    @classmethod
    def from_exons(cls, exons: list[ExonBlock], region_start: int, region_end: int) -> "_CompressedMapper":
        merged = _merge_intervals([(exon.start, exon.end) for exon in exons])
        raw_segments: list[tuple[int, int, float]] = []
        cursor = region_start
        for start, end in merged:
            if start > cursor:
                raw_segments.append((cursor, start - 1, min(8.0, max(2.5, math.log10(max(2, start - cursor)) * 2.5))))
            raw_segments.append((start, end, max(3.4, min(11.5, (end - start + 1) / 115.0))))
            cursor = end + 1
        if cursor < region_end:
            raw_segments.append((cursor, region_end, min(8.0, max(2.5, math.log10(max(2, region_end - cursor + 1)) * 2.5))))
        total = sum(item[2] for item in raw_segments) or 1.0
        segments: list[tuple[int, int, float, float]] = []
        pos = 0.0
        for start, end, width in raw_segments:
            next_pos = pos + width / total * 100.0
            segments.append((start, end, pos, next_pos))
            pos = next_pos
        return cls(segments)

    def map(self, position: int | float) -> float:
        pos = float(position)
        for start, end, x1, x2 in self.segments:
            if start <= pos <= end:
                if end <= start:
                    return x1
                return x1 + (pos - start) / (end - start) * (x2 - x1)
        if pos < self.segments[0][0]:
            return 0.0
        return 100.0

    def inverse(self, x: float) -> float:
        for start, end, x1, x2 in self.segments:
            if x1 <= x <= x2:
                if x2 <= x1:
                    return float(start)
                return start + (x - x1) / (x2 - x1) * (end - start)
        return float(self.segments[-1][1])


class _LinearMapper:
    def __init__(self, region_start: int, region_end: int) -> None:
        self.region_start = int(region_start)
        self.region_end = max(int(region_end), self.region_start + 1)

    def map(self, position: int | float) -> float:
        span = max(1.0, float(self.region_end - self.region_start))
        return max(0.0, min(100.0, (float(position) - self.region_start) / span * 100.0))

    def inverse(self, x: float) -> float:
        span = max(1.0, float(self.region_end - self.region_start))
        return self.region_start + max(0.0, min(100.0, float(x))) / 100.0 * span


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    sorted_intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged: list[tuple[int, int]] = []
    for start, end in sorted_intervals:
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def _block(row: pd.Series, start_col: str, end_col: str, label: str, role: str, required: bool = True) -> ExonBlock | None:
    if start_col not in row.index or end_col not in row.index or pd.isna(row[start_col]) or pd.isna(row[end_col]):
        if required:
            raise RuntimeError(f"Missing rMATS coordinate columns: {start_col}, {end_col}")
        return None
    start = int(float(row[start_col]))
    end = int(float(row[end_col]))
    if start_col.endswith("_0base") or start_col.endswith("ES") or start_col in {"shortES", "flankingES"}:
        start += 1
    return ExonBlock(min(start, end), max(start, end), label, role)


def _junction(left: ExonBlock, right: ExonBlock, label: str, role: str) -> JunctionDef:
    ordered = sorted([left, right], key=lambda exon: exon.start)
    return JunctionDef(ordered[0].end + 1, ordered[1].start - 1, label, role)


def _junction_anchors(ops: list[tuple[str, str]], n_index: int) -> tuple[int, int]:
    left_anchor = 0
    right_anchor = 0
    for i in range(n_index - 1, -1, -1):
        length, op = ops[i]
        if op in {"M", "=", "X"}:
            left_anchor = int(length)
            break
        if op in {"N", "S", "H"}:
            break
    for i in range(n_index + 1, len(ops)):
        length, op = ops[i]
        if op in {"M", "=", "X"}:
            right_anchor = int(length)
            break
        if op in {"N", "S", "H"}:
            break
    return left_anchor, right_anchor


def _normalize_bamnostic_cigar_item(item: tuple[int, int] | tuple[str, str]) -> tuple[int, int]:
    first, second = item
    try:
        first_int = int(first)
        second_int = int(second)
    except (TypeError, ValueError):
        return 0, 0
    # bamnostic exposes BAM-style cigartuples as (operation_code, length).
    if first_int in {0, 1, 2, 3, 4, 5, 6, 7, 8}:
        return first_int, second_int
    return second_int, first_int


def _junction_anchors_bamnostic(ops: list[tuple[int, int]], n_index: int) -> tuple[int, int]:
    left_anchor = 0
    right_anchor = 0
    for i in range(n_index - 1, -1, -1):
        op, length = _normalize_bamnostic_cigar_item(ops[i])
        if op in {0, 7, 8}:
            left_anchor = int(length)
            break
        if op in {3, 4, 5}:
            break
    for i in range(n_index + 1, len(ops)):
        op, length = _normalize_bamnostic_cigar_item(ops[i])
        if op in {0, 7, 8}:
            right_anchor = int(length)
            break
        if op in {3, 4, 5}:
            break
    return left_anchor, right_anchor


def _parse_psi_list(value: object) -> list[float]:
    if value is None or pd.isna(value):
        return []
    values: list[float] = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if not token or token.upper() == "NA":
            continue
        try:
            values.append(float(token))
        except ValueError:
            continue
    return values


def _parse_number_list(value: object) -> list[float]:
    if value is None or pd.isna(value):
        return []
    values: list[float] = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if not token or token.upper() in {"NA", "NAN"}:
            continue
        try:
            values.append(float(token))
        except ValueError:
            continue
    return values


def _to_float(value: object) -> float:
    try:
        if value is None or pd.isna(value):
            return float("nan")
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _fmt_float(value: object) -> str:
    number = _to_float(value)
    if np.isnan(number):
        return "NA"
    return f"{number:.3g}"


def _fmt_signed(value: object) -> str:
    number = _to_float(value)
    if np.isnan(number):
        return "NA"
    return f"{number:+.3g}"


def _mean_or_nan(values: list[float]) -> float:
    clean = [float(value) for value in values if pd.notna(value)]
    return float(np.nanmean(clean)) if clean else float("nan")


def _rmats_track(label: object, ijc: list[float], sjc: list[float], psi: list[float], group_index: int) -> dict[str, object]:
    return {
        "label": str(label),
        "ijc": ijc,
        "sjc": sjc,
        "mean_ijc": _mean_or_nan(ijc),
        "mean_sjc": _mean_or_nan(sjc),
        "psi": _mean_or_nan(psi),
        "replicate_psis": psi,
        "group_index": group_index,
    }


def _mean_track(label: object, tracks: list[dict[str, object]], psis: list[float]) -> dict[str, object]:
    if tracks:
        coverage = np.vstack([track["coverage"] for track in tracks]).mean(axis=0)
    else:
        coverage = np.zeros(240, dtype=float)
    junction_labels = set()
    for track in tracks:
        junction_labels.update(track["junction_counts"].keys())
    counts = {
        label_: int(round(float(np.mean([track["junction_counts"].get(label_, 0) for track in tracks]))))
        for label_ in junction_labels
    }
    return {
        "label": str(label),
        "coverage": coverage,
        "junction_counts": counts,
        "psi": float(np.nanmean(psis)) if psis else np.nan,
        "replicate_psis": psis,
        "group_index": 0,
    }


def _draw_psi_strip(ax, mean: float, reps: object, color: str) -> None:
    center_x = 1.055
    left, right = 1.000, 1.110
    if np.isnan(mean):
        ax.text(center_x, 0.86, "PSI: NA", transform=ax.transAxes, ha="center", va="center", fontsize=12, fontweight="bold")
        return
    ax.text(center_x, 0.86, f"PSI: {mean:.2f}", transform=ax.transAxes, ha="center", va="center", fontsize=12.5, fontweight="bold")
    values = [float(value) for value in reps if pd.notna(value)] if isinstance(reps, list) else []
    if not values:
        return
    base = 0.755
    ax.plot([left, right], [base, base], transform=ax.transAxes, color="#CBD5E1", linewidth=1.1, clip_on=False)
    for value in values:
        px = left + (right - left) * max(0.0, min(1.0, value))
        ax.scatter([px], [base], transform=ax.transAxes, s=24, color=color, edgecolor="white", linewidth=0.6, clip_on=False, zorder=5)
    mx = left + (right - left) * max(0.0, min(1.0, mean))
    ax.plot([mx, mx], [base - 0.025, base + 0.035], transform=ax.transAxes, color="#111827", linewidth=1.3, clip_on=False)
    ax.text(center_x, base - 0.07, "replicate PSI", transform=ax.transAxes, ha="center", va="top", fontsize=8, color="#4B5563")


def _display_exon_label(exon: ExonBlock) -> str:
    label = str(exon.label or "").strip()
    if exon.role == "skipped" and label.lower() == "skipped":
        return "skipped exon"
    return label or "exon"


def _exon_group_center(exons: list[ExonBlock], mapper: "_CompressedMapper") -> float:
    if not exons:
        return 50.0
    left = min(mapper.map(exon.start) for exon in exons)
    right = max(mapper.map(exon.end) for exon in exons)
    return (left + right) / 2.0


def _display_exon_role(exon: ExonBlock) -> str:
    role_names = {
        "constitutive": "constitutive exon",
        "event": "event exon",
        "alternative": "alternative exon",
        "retained": "retained intron",
        "skipped": "skipped exon",
    }
    role = role_names.get(exon.role, exon.role or "exon")
    label = str(exon.label or "")
    if exon.role == "skipped":
        if _label_has_multiple_exons(label):
            return "skipped exons"
        if not re.search(r"\bE\d+\b", label):
            return f"{role} | {exon.start:,}-{exon.end:,}"
    return role


def _parse_gtf_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def _gtf_exon_label(attrs: dict[str, str]) -> str:
    exon_number = str(attrs.get("exon_number", "")).strip()
    if exon_number:
        return f"E{exon_number}"
    exon_id = str(attrs.get("exon_id", "")).strip()
    if exon_id:
        return exon_id.split(".", 1)[0]
    return "annotated exon"


def _collapse_gtf_labels(labels: set[str]) -> str:
    clean = sorted({label for label in labels if label}, key=_label_sort_key)
    if not clean:
        return "annotated exon"
    if len(clean) <= 3:
        return "/".join(clean)
    return "/".join(clean[:3]) + "+"


def _label_sort_key(label: str) -> tuple[int, str]:
    match = re.search(r"\d+", label)
    if match:
        return int(match.group(0)), label
    return 10**9, label


def _matching_gtf_exon_label(start: int, end: int, gtf_exons_by_transcript: dict[str, list[tuple[int, int, str]]]) -> str:
    best_multi_labels: list[str] = []
    best_multi_score: tuple[int, int, int] = (0, 0, -10**9)
    best_single_label = ""
    best_single_score = 0.0
    query_length = max(1, end - start + 1)
    for transcript_exons in gtf_exons_by_transcript.values():
        labels: list[str] = []
        covered_exon_bp = 0
        boundary_distance = 10**9
        for exon_start, exon_end, label in transcript_exons:
            overlap = max(0, min(end, exon_end) - max(start, exon_start) + 1)
            if overlap <= 0:
                continue
            exon_length = max(1, exon_end - exon_start + 1)
            exon_fraction = overlap / exon_length
            query_fraction = overlap / query_length
            exact_match = exon_start == start and exon_end == end
            if exon_fraction >= 0.50 or exact_match:
                labels.extend(_split_exon_labels(label))
                covered_exon_bp += overlap
                boundary_distance = min(boundary_distance, abs(exon_start - start) + abs(exon_end - end))
            if exact_match or query_fraction >= 0.55:
                single_score = query_fraction + (1.0 if exact_match else 0.0)
                if single_score > best_single_score:
                    best_single_score = single_score
                    best_single_label = _format_exon_label_set(_split_exon_labels(label))
        unique_labels = sorted({label for label in labels if label}, key=_label_sort_key)
        if len(unique_labels) >= 2:
            # Multi-exon skipped blocks include introns, so genomic-span coverage is
            # intentionally not used here; the robust signal is multiple annotated
            # exons from the same transcript lying inside the event block.
            score = (len(unique_labels), covered_exon_bp, -boundary_distance)
            if score > best_multi_score:
                best_multi_score = score
                best_multi_labels = unique_labels
    if best_multi_labels:
        return _format_exon_label_set(best_multi_labels)
    return best_single_label


def _split_exon_labels(label: str) -> list[str]:
    return [token for token in re.split(r"[/,\s]+", str(label or "")) if token]


def _format_exon_label_set(labels: list[str]) -> str:
    clean = sorted({label for label in labels if label}, key=_label_sort_key)
    if not clean:
        return ""
    numbers: list[int] = []
    non_numeric: list[str] = []
    for label in clean:
        match = re.fullmatch(r"E(\d+)", label)
        if match:
            numbers.append(int(match.group(1)))
        else:
            non_numeric.append(label)
    if numbers and not non_numeric:
        numbers = sorted(set(numbers))
        if len(numbers) == 1:
            return f"E{numbers[0]}"
        if numbers == list(range(numbers[0], numbers[-1] + 1)):
            return f"E{numbers[0]}-E{numbers[-1]}"
        return "/".join(f"E{number}" for number in numbers)
    return "/".join(clean[:4]) + ("+" if len(clean) > 4 else "")


def _label_has_multiple_exons(label: str) -> bool:
    return bool(re.search(r"E\d+\s*-\s*E\d+", label) or len(re.findall(r"\bE\d+\b", label)) > 1)
