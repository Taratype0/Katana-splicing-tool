args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7L) stop('Expected design, output, refA, refB, library, min_reads, min_replicates')
if (nzchar(args[5])) .libPaths(c(args[5], .libPaths()))
suppressPackageStartupMessages(library(DEXSeq))
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
meta <- read.delim(args[1], stringsAsFactors=FALSE, check.names=FALSE)
refA <- args[3]; refB <- args[4]
stopifnot(!anyDuplicated(meta$sample_id), length(unique(meta$factor_a))==2L,
          length(unique(meta$factor_b))==2L, refA %in% meta$factor_a, refB %in% meta$factor_b)
counts_by_sample <- lapply(meta$count_file, function(path) {
  x <- read.delim(path, header=FALSE, col.names=c('bin_id','count'), stringsAsFactors=FALSE)
  x <- x[!grepl('^_', x$bin_id), ]
  if (anyDuplicated(x$bin_id) || any(!grepl(':',x$bin_id))) stop('Invalid or duplicate exon bin IDs')
  if (anyNA(x$count) || any(!is.finite(x$count)) || any(x$count < 0) || any(x$count != floor(x$count))) stop('Counts must be nonnegative integers')
  x
})
ids <- counts_by_sample[[1]]$bin_id
mat <- vapply(counts_by_sample, function(x) {
  if (!setequal(ids, x$bin_id)) stop('Exon bin universe differs between samples')
  x$count[match(ids, x$bin_id)]
}, numeric(length(ids)))
rownames(mat) <- ids; colnames(mat) <- meta$sample_id
gene <- sub(':[^:]+$', '', ids)
feature <- sub('^.*:', 'E', ids)
meta$factor_a <- factor(meta$factor_a, levels=c(refA, setdiff(unique(meta$factor_a),refA)))
meta$factor_b <- factor(meta$factor_b, levels=c(refB, setdiff(unique(meta$factor_b),refB)))
rownames(meta) <- meta$sample_id
if (any(table(meta$factor_a,meta$factor_b)<2L)) stop('Incomplete or unreplicated 2x2 design')
cat('Constructing DEXSeqDataSet:', nrow(mat), 'bins;', ncol(mat), 'samples\n')
dxd <- DEXSeqDataSet(mat, sampleData=meta, design=~sample+exon, featureID=feature, groupID=gene)
# Numeric usage terms avoid sample-level aliases in the expanded this/others design.
cd <- colData(dxd)
is_this <- as.numeric(cd$exon=='this')
colData(dxd)$a_usage <- as.numeric(cd$factor_a != refA)*is_this
colData(dxd)$b_usage <- as.numeric(cd$factor_b != refB)*is_this
colData(dxd)$ab_usage <- colData(dxd)$a_usage*colData(dxd)$b_usage
full <- ~sample+exon+a_usage+b_usage+ab_usage
reduced <- ~sample+exon+a_usage+b_usage
if ('batch' %in% names(meta)) {
  if (anyNA(meta$batch) || any(meta$batch=='')) stop('Incomplete batch metadata')
  if (length(unique(meta$batch))>1L) {
    batches <- model.matrix(~factor(batch), as.data.frame(colData(dxd)))[,-1,drop=FALSE]
    names_added <- paste0('batch_usage_',seq_len(ncol(batches)))
    for (i in seq_along(names_added)) colData(dxd)[[names_added[i]]] <- batches[,i]*is_this
    full <- reformulate(c('sample','exon','a_usage','b_usage','ab_usage',names_added))
    reduced <- reformulate(c('sample','exon','a_usage','b_usage',names_added))
  }
}
for (formula in list(full,reduced)) {
  mm <- model.matrix(formula, as.data.frame(colData(dxd)))
  if (qr(mm)$rank != ncol(mm)) stop('Experimental design is confounded or rank deficient')
}
design(dxd) <- full
supported <- rowSums(mat>=as.numeric(args[6])) >= as.integer(args[7])
multiple_bins <- as.integer(table(gene)[gene]) >= 2L
keep <- supported & multiple_bins
audit <- data.frame(bin_id=ids,gene_id=gene,read_supported=supported,multiple_bins=multiple_bins,test_eligible=keep)
write.table(audit,file.path(out,'filter_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if (!any(keep)) stop('No eligible exon bins; thresholds were not relaxed')
dxd <- dxd[keep,]
cat('Eligible bins:',sum(keep),'Estimating size factors\n')
dxd <- estimateSizeFactors(dxd)
cat('Estimating dispersions with full design\n')
checkpoint <- file.path(out,'dispersions_checkpoint.rds')
signature <- list(inputs=tools::md5sum(c(args[1],meta$count_file)),
                  references=c(refA,refB), eligibility=args[6:7],
                  model=deparse(full), DEXSeq=as.character(packageVersion('DEXSeq')),
                  DESeq2=as.character(packageVersion('DESeq2')), R=R.version.string)
cached <- if (file.exists(checkpoint)) tryCatch(readRDS(checkpoint),error=function(e) NULL) else NULL
if (!is.null(cached) && identical(cached$signature,signature)) {
  dxd <- cached$dxd
  cat('Dispersion checkpoint: cache hit\n')
} else {
  dxd <- estimateDispersions(dxd,formula=full,BPPARAM=BiocParallel::SerialParam())
  saveRDS(list(signature=signature,dxd=dxd),paste0(checkpoint,'.tmp'))
  if (file.exists(checkpoint)) unlink(checkpoint)
  if (!file.rename(paste0(checkpoint,'.tmp'),checkpoint)) stop('Could not commit dispersion checkpoint')
  cat('Dispersion checkpoint saved\n')
}
cat('Testing one genotype-by-treatment exon-usage interaction\n')
dxd <- testForDEU(dxd,fullModel=full,reducedModel=reduced,BPPARAM=BiocParallel::SerialParam())
results <- as.data.frame(DEXSeqResults(dxd,independentFiltering=FALSE))
results$bin_id <- ids[keep]
results$interaction_FDR <- p.adjust(results$pvalue,method='BH')
# Fractions are descriptive exon-bin / gene-bin counts, NOT rMATS PSI.
totals <- rowsum(mat,gene,reorder=FALSE)
fraction <- mat / totals[match(gene,rownames(totals)),,drop=FALSE]
means <- sapply(c('00','10','01','11'), function(cell) {
  select <- as.integer(meta$factor_a)-1L == as.integer(substr(cell,1,1)) &
            as.integer(meta$factor_b)-1L == as.integer(substr(cell,2,2))
  rowMeans(fraction[,select,drop=FALSE],na.rm=TRUE)
})
results$descriptive_fraction_interaction <- (means[,'11']-means[,'01']-means[,'10']+means[,'00'])[keep]
results$row_id <- results$bin_id
flat <- vapply(results,function(x) is.atomic(x) && is.null(dim(x)),logical(1))
write.table(results[,flat,drop=FALSE],file.path(out,'dexseq_interaction.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(bin_id=ids[keep],fraction[keep,,drop=FALSE],check.names=FALSE),file.path(out,'exon_bin_fractions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(dxd,file.path(out,'DEXSeq_interaction.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
writeLines(c(paste('Full:',deparse(full)),paste('Reduced:',deparse(reduced)),
 'One 2x2 interaction; BH across all valid tested exon-bin p-values.',
 'Exon-bin fraction is descriptive and is not junction-derived PSI.'),file.path(out,'model.txt'))
cat('DEXSeq interaction complete\n')
