#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

# Parse command arguments

argv <- commandArgs(trailingOnly = TRUE)

argc <- length(argv)

samples <- argv[1]

contrasts <- argv[2]

counts <- argv[3]

ntop <- as.integer(argv[4])


# Load required packages

library(BiocParallel)

library(DEXSeq)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_path)), "katana_dexseq_design.R"))

workers <- Sys.getenv("R_BIOC_WORKERS", "")
if (nzchar(workers)) {
    workers <- suppressWarnings(as.integer(workers))
    if (is.na(workers) || workers < 1L) stop("R_BIOC_WORKERS must be positive")
    register(if (workers == 1L) SerialParam() else MulticoreParam(workers = workers))
}


# Define helper functions

splitByContrast <- function(object, contrast) {

    keep <- colData(object)$condition %in% contrast

    object <- object[, keep, drop = FALSE]

    colData(object)$condition <- factor(
        x = colData(object)$condition,
        levels = contrast
    )

    object

}

DEXSeq <- function(object) {

    object <- estimateSizeFactors(object)

    object <- estimateDispersions(object)

    object <- testForDEU(object)

    object <- estimateExonFoldChanges(object)

    object

}

keepValidColumns <- function(x) {

    keep <- colnames(x)[sapply(x, class) %in% c("numeric", "character")]

    x <- x[ , colnames(x) %in% keep, drop = FALSE]

    x

}

vectorToDataFrame <- function(x) {

    x <- data.frame(groupID = names(x), padj = unname(x))

    x

}


# Read samples table

samples <- read.delim(samples, stringsAsFactors = TRUE, check.names = FALSE)

colnames(samples)[colnames(samples) == "sample_id"] <- "sample"


# Read contrasts table

contrasts <- read.csv(contrasts, check.names = FALSE)

contrasts <- contrasts[, intersect(c("contrast", "treatment", "control", "paired", "design"), names(contrasts)), drop = FALSE]
if (Sys.getenv("KANATA_SINGLE_CONTRAST", "") == "1" && nrow(contrasts) != 1L)
    stop("One contrast is required per Nextflow task")


# Read counts table

counts <- read.table(counts, sep = "\t", header = TRUE, check.names = FALSE)

annotation <- data.frame(
    featureID = counts$feature_id,
    groupID = counts$gene_id
)

keep <- seq(3, ncol(counts), by = 1)

counts <- counts[, keep, drop = FALSE]

counts <- round(counts)


# Reorder count matrix by sample table

stopifnot(all(samples$sample %in% colnames(counts)))

stopifnot(all(colnames(counts) %in% samples$sample))

counts <- counts[, samples$sample, drop = FALSE]


# Create DEXSeqDataSet object

object <- DEXSeqDataSet(
    countData = counts,
    sampleData = samples,
    featureID = annotation$featureID,
    groupID = annotation$groupID
)


# Define contrasts

groups <- levels(colData(object)$condition)

contrast_records <- contrasts
contrasts <- data.frame(A = contrasts$treatment, B = contrasts$control)

contrasts <- contrasts[contrasts$A != contrasts$B, , drop = FALSE]

contrasts <- asplit(contrasts, MARGIN = 1)

names <- sapply(contrasts, paste, collapse = "-")


# Create objects list

objects <- mapply(
    FUN = splitByContrast,
    contrast = contrasts,
    MoreArgs = list(object = object),
    SIMPLIFY = FALSE
)


# Run DEXSeq workflow

objects <- lapply(seq_along(objects), function(i) katana_dexseq_fit(objects[[i]], contrast_records[i, , drop=FALSE]))

names(objects) <- names


# Run DEXSeqResults

results.exon <- lapply(objects, DEXSeqResults)

names(results.exon) <- names


# Run perGeneQValue

results.gene <- lapply(results.exon, perGeneQValue)

names(results.gene) <- names


# Save objects to disk

mapply(
    saveRDS,
    object = objects,
    file = paste0("DEXSeqDataSet.", names, ".rds")
)

mapply(
    saveRDS,
    object = results.exon,
    file = paste0("DEXSeqResults.", names, ".rds")
)

mapply(
    saveRDS,
    object = results.gene,
    file = paste0("perGeneQValue.", names, ".rds")
)


# Save results to disk

mapply(
    write.table,
    x = lapply(results.exon, keepValidColumns),
    file = paste0("DEXSeqResults.", names, ".tsv"),
    MoreArgs = list(quote = FALSE, sep = "\t", row.names = FALSE)
)

mapply(
    write.table,
    x = lapply(results.gene, vectorToDataFrame),
    file = paste0("perGeneQValue.", names, ".tsv"),
    MoreArgs = list(quote = FALSE, sep = "\t", row.names = FALSE)
)


# Print session information

citation("DEXSeq")

sessionInfo()
