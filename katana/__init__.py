"""Katana Unified Python Suite — Precision-Aware RNA Splicing, Network, Cis-Motif & Proteome Platform.

Supports two independent yet seamlessly converging workflows:
1. Standalone Python Package / API (`katana.run_bulk`, `katana.run_single_cell`, `katana.open_workspace`):
   Run directly from completed upstream splicing/expression tables (e.g. rMATS, DESeq2, Smart-seq3 junction matrices)
   with zero Nextflow, Docker, or R requirement, producing full Parquet/TSV outputs, publication PNG/SVG figures,
   and a self-contained portable `REPORT.html` + `REPORT.zip`.
2. End-to-End Automated Nextflow Pipeline (`katana.run_pipeline` / `nextflow run pipeline/main.nf`):
   Start directly from raw `FASTQ` or aligned `BAM` files, execute upstream alignment/QC/MultiQC/quantification,
   and converge directly into `katana.run_bulk` for identical downstream statistical inference and reporting.
"""
from __future__ import annotations

from .suite import open_workspace, run_bulk, run_pipeline, run_single_cell

__version__ = "2.7.0"

__all__ = [
    "__version__",
    "open_workspace",
    "run_bulk",
    "run_pipeline",
    "run_single_cell",
]
