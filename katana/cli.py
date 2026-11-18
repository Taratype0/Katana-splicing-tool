"""Unified `katana` CLI: Standalone Python Suite (`run-bulk`, `run-sc`, `view`) + Nextflow FASTQ/BAM Pipeline (`pipeline`)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from . import __version__
from .suite import open_workspace, run_bulk, run_pipeline, run_single_cell


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="katana",
        description=(
            "Katana v2.7 Unified RNA Splicing Suite: run directly as a pure-Python package on rMATS / Single-Cell "
            "tables (`run-bulk`, `run-sc`), inspect interactively (`view`), or execute end-to-end from FASTQ / BAM (`pipeline`)."
        ),
    )
    parser.add_argument("--version", action="version", version=f"katana {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # 1. katana run-bulk
    p_bulk = sub.add_parser(
        "run-bulk",
        help="Run the standalone Python Bulk Splicing Suite (PS-DAS, Beyond-Jutils PCA/Clustermap, SpliceNet, SpliceTrace 170-PWM, ORA & REPORT.zip).",
    )
    p_bulk.add_argument("--bundle", type=Path, default=None, help="Completed Nextflow upstream bundle directory (converges FASTQ/BAM pipeline into Python suite)")
    p_bulk.add_argument("--rmats", type=str, default=None, help="rMATS directory (single contrast dir, parent dir of contrasts, or c1=dir1,c2=dir2)")
    p_bulk.add_argument("--samples", type=Path, default=None, help="Sample metadata TSV/CSV (sample_id, condition[, biological_sample_id])")
    p_bulk.add_argument("--contrasts", type=Path, default=None, help="Optional contrasts TSV/CSV (contrast_id, treatment, control)")
    p_bulk.add_argument("--deseq2", type=Path, default=None, help="Optional DESeq2 directory or TSV")
    p_bulk.add_argument("--multiqc", type=Path, default=None, help="Optional MultiQC directory for instant 5'-3' coverage extraction")
    p_bulk.add_argument("--fasta", type=Path, default=None, help="Optional genome FASTA for 11-channel SpliceTrace & mature ORF/NMD translation")
    p_bulk.add_argument("--gtf", type=Path, default=None, help="Optional GTF annotation for SpliceTrace & domain mapping")
    p_bulk.add_argument("--bams", type=Path, default=None, help="Optional BAM directory or sample->BAM TSV for on-demand Sashimi")
    p_bulk.add_argument("--study-name", default="Katana Bulk Splicing Study", help="Study display title")
    p_bulk.add_argument("--species", choices=["auto", "human", "mouse"], default="auto", help="Species namespace (default: auto)")
    p_bulk.add_argument("--fdr", type=float, default=0.05, help="Significance FDR threshold (default: 0.05)")
    p_bulk.add_argument("--dpsi", type=float, default=0.10, help="Minimum |ΔPSI| effect size (default: 0.10)")
    p_bulk.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory for tables, figures, REPORT.html & REPORT.zip")

    # 2. katana run-sc
    p_sc = sub.add_parser(
        "run-sc",
        help="Run the standalone Python Single-Cell Splicing Suite (Precision-Kernel Event-PSI UMAP, SC-SpliceNet, ORA & REPORT.zip).",
    )
    p_sc.add_argument("--events", type=Path, required=True, help="Single-cell splicing contrast table (.parquet/.tsv/.csv)")
    p_sc.add_argument("--cells", type=Path, required=True, help="Cell metadata table (cell_id, condition, donor_id[, cell_type])")
    p_sc.add_argument("--psi-matrix", type=Path, default=None, help="Optional event x cell PSI matrix (.parquet/.tsv)")
    p_sc.add_argument("--support-matrix", type=Path, default=None, help="Optional event x cell junction read-support matrix (.parquet/.tsv)")
    p_sc.add_argument("--expression", type=Path, default=None, help="Optional gene x cell logTPM matrix for cross-modality divergence Δ_kNN")
    p_sc.add_argument("--study-name", default="Katana Single-Cell Splicing Study", help="Study display title")
    p_sc.add_argument("--species", choices=["human", "mouse"], default="human", help="Species namespace")
    p_sc.add_argument("--fdr", type=float, default=0.05, help="Significance FDR threshold")
    p_sc.add_argument("--dpsi", type=float, default=0.10, help="Minimum |ΔPSI| effect size")
    p_sc.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory for tables, figures, REPORT.html & REPORT.zip")

    # 3. katana view
    p_view = sub.add_parser(
        "view",
        help="Launch the local interactive workspace & live compute labs on a REPORT.zip or report directory.",
    )
    p_view.add_argument("report", nargs="?", default=None, help="Optional path to REPORT.zip or report directory to pre-open")
    p_view.add_argument("--host", default="127.0.0.1", help="Bind host (default: 127.0.0.1)")
    p_view.add_argument("--port", type=int, default=18775, help="Bind port (default: 18775)")
    p_view.add_argument("--no-browser", action="store_true", help="Start server without opening web browser")

    # 4. katana pipeline
    p_pipe = sub.add_parser(
        "pipeline",
        help="Run the end-to-end automated Nextflow pipeline starting from FASTQ or BAM (`pipeline/main.nf`), with optional convergence into `run-bulk`.",
    )
    p_pipe.add_argument("--nextflow", default="nextflow", help="Nextflow executable (default: nextflow)")
    p_pipe.add_argument("--suite-outdir", type=Path, default=None, help="Optional directory to automatically run `katana run-bulk` on the completed upstream bundle")
    p_pipe.add_argument("nextflow_args", nargs=argparse.REMAINDER, help="Arguments passed directly to `nextflow run pipeline/main.nf` (e.g. --mode complete --input samples.csv ...)")

    args = parser.parse_args(argv)
    try:
        if args.command == "run-bulk":
            receipt = run_bulk(
                outdir=args.outdir,
                bundle=args.bundle,
                rmats=args.rmats,
                samples=args.samples,
                contrasts=args.contrasts,
                deseq2=args.deseq2,
                multiqc=args.multiqc,
                fasta=args.fasta,
                gtf=args.gtf,
                bams=args.bams,
                study_name=args.study_name,
                species=args.species,
                fdr_cutoff=args.fdr,
                dpsi_cutoff=args.dpsi,
            )
            print(json.dumps(receipt, indent=2))
            return 0
        if args.command == "run-sc":
            receipt = run_single_cell(
                outdir=args.outdir,
                events=args.events,
                cells=args.cells,
                psi_matrix=args.psi_matrix,
                support_matrix=args.support_matrix,
                expression=args.expression,
                study_name=args.study_name,
                species=args.species,
                fdr_cutoff=args.fdr,
                dpsi_cutoff=args.dpsi,
            )
            print(json.dumps(receipt, indent=2))
            return 0
        if args.command == "view":
            open_workspace(args.report, host=args.host, port=args.port, open_browser=not args.no_browser)
            return 0
        if args.command == "pipeline":
            nf_args = [a for a in args.nextflow_args if a != "--"]
            return run_pipeline(nf_args, nextflow_executable=args.nextflow, converge_suite_outdir=args.suite_outdir)
    except Exception as exc:
        print(json.dumps({"status": "failed", "error": str(exc)}), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
