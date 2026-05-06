from __future__ import annotations

from pathlib import Path

import subprocess


class IsoformSwitchRunner:
    def __init__(self, rscript_path: Path | None = None) -> None:
        self.rscript_path = rscript_path

    def write_script(self, output_dir: Path, manifest_file: Path, design_file: Path, gtf_file: Path, fasta_file: Path | None) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        script = output_dir / "run_isoformswitch.R"
        script.write_text(
            "\n".join(
                [
                    "args <- commandArgs(trailingOnly = TRUE)",
                    "manifest <- args[1]",
                    "design_file <- args[2]",
                    "gtf_file <- args[3]",
                    "fasta_file <- args[4]",
                    "outdir <- args[5]",
                    "dir.create(outdir, recursive = TRUE, showWarnings = FALSE)",
                    "suppressPackageStartupMessages({",
                    "  library(IsoformSwitchAnalyzeR)",
                    "  library(readr)",
                    "})",
                    "manifest_df <- read.delim(manifest, sep='\\t', stringsAsFactors = FALSE)",
                    "design_df <- read.delim(design_file, sep='\\t', stringsAsFactors = FALSE)",
                    "sample_vector <- manifest_df$quant_sf",
                    "names(sample_vector) <- manifest_df$sample_id",
                    "iso_quant <- importIsoformExpression(sampleVector = sample_vector)",
                    "fasta_arg <- if (nchar(fasta_file) > 0 && file.exists(fasta_file)) fasta_file else NULL",
                    "switch_list <- importRdata(",
                    "  isoformCountMatrix = iso_quant$counts,",
                    "  isoformRepExpression = iso_quant$abundance,",
                    "  designMatrix = design_df,",
                    "  isoformExonAnnoation = gtf_file,",
                    "  isoformNtFasta = fasta_arg,",
                    "  quiet = TRUE",
                    ")",
                    "switch_list <- preFilter(",
                    "  switch_list,",
                    "  geneExpressionCutoff = 1,",
                    "  isoformExpressionCutoff = 0,",
                    "  removeSingleIsoformGenes = TRUE",
                    ")",
                    "switch_list <- isoformSwitchTestDEXSeq(",
                    "  switchAnalyzeRlist = switch_list,",
                    "  quiet = TRUE",
                    ")",
                    "switch_list <- extractTopSwitches(",
                    "  switchAnalyzeRlist = switch_list,",
                    "  n = 100,",
                    "  filterForConsequences = FALSE",
                    ")",
                    "write.table(",
                    "  switch_list,",
                    "  file = file.path(outdir, 'isoform_switch_top.tsv'),",
                    "  sep = '\\t',",
                    "  quote = FALSE,",
                    "  row.names = FALSE",
                    ")",
                    "writeLines(c(",
                    "  'IsoformSwitchAnalyzeR completed',",
                    "  paste('Manifest:', manifest),",
                    "  paste('Design:', design_file),",
                    "  paste('GTF:', gtf_file)",
                    "), con = file.path(outdir, 'isoformswitch_run.log'))",
                ]
            ),
            encoding="utf-8",
        )
        return script

    def run(
        self,
        script_path: Path,
        manifest_file: Path,
        design_file: Path,
        gtf_file: Path,
        fasta_file: Path | None,
        output_dir: Path,
    ) -> subprocess.CompletedProcess:
        if self.rscript_path is None:
            raise RuntimeError("Rscript is not configured. Set the Rscript path in Settings.")
        return subprocess.run(
            [
                str(self.rscript_path),
                str(script_path),
                str(manifest_file),
                str(design_file),
                str(gtf_file),
                str(fasta_file) if fasta_file else "",
                str(output_dir),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
