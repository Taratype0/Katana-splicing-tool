# Paired DEXSeq adaptation: subject-specific exon usage in BOTH LRT models.
# Native estimateExonFoldChanges remains a descriptive condition-only estimate.
katana_dexseq_fit <- function(object, record) {
    paired <- "paired" %in% names(record) && tolower(as.character(record$paired)) %in% c("true", "1", "yes")
    if (!paired) return(DEXSeq(object))  # Preserve the existing independent route.
    full <- ~ sample + exon + condition:exon
    reduced <- ~ sample + exon
    if (paired) {
        cd <- as.data.frame(colData(object))
        if (!"subject" %in% names(cd) || any(is.na(cd$subject) | !nzchar(as.character(cd$subject))))
            stop("Paired DEXSeq requires subject identity for every sample")
        units <- unique(cd[, c("sample", "condition", "subject")])
        if (any(table(units$subject, units$condition) != 1L))
            stop("Paired DEXSeq requires exactly one sample per subject per condition")
        colData(object)$subject <- droplevels(factor(colData(object)$subject))
        colData(object)$sample <- droplevels(factor(colData(object)$sample))
        full <- ~ sample + exon + subject:exon + condition:exon
        reduced <- ~ sample + exon + subject:exon
    }
    design(object) <- full
    cd <- as.data.frame(colData(object))
    full_rank <- qr(model.matrix(full, cd))$rank
    reduced_rank <- qr(model.matrix(reduced, cd))$rank
    if (full_rank >= nrow(cd) || full_rank - reduced_rank != 1L)
        stop("DEXSeq condition effect is not estimable with residual degrees of freedom")
    object <- estimateSizeFactors(object)
    object <- estimateDispersions(object, formula=full)
    object <- testForDEU(object, fullModel=full, reducedModel=reduced)
    object <- estimateExonFoldChanges(object)
    metadata(object)$katana_design <- list(
        paired=paired, full=paste(deparse(full), collapse=" "),
        reduced=paste(deparse(reduced), collapse=" "),
        native_fold_change="condition-only descriptive estimate; paired adjustment applies to the likelihood-ratio P/FDR",
        treatment=as.character(record$treatment), control=as.character(record$control))
    object
}
