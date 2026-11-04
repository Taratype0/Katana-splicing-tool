#!/usr/bin/env Rscript
# Generated counts only. Exercise the installed DEXSeq adapter and DESeq2 GLM.
args <- commandArgs(trailingOnly=TRUE)
bin <- normalizePath(args[1]); out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE); setwd(out)
suppressPackageStartupMessages(library(DEXSeq))
set.seed(20261028)
samples <- data.frame(sample=paste0("S",1:6), condition=rep(c("control","treated"), each=3),
                      subject=rep(paste0("mouse",1:3),2))
comparison <- data.frame(contrast="treated_vs_control", treatment="treated",control="control",
                         paired=TRUE, design="~ subject + condition")
write.csv(samples,"samples.csv",row.names=FALSE,quote=FALSE)
write.table(setNames(samples,c("sample_id","condition","subject")),"samples.tsv",sep="\t",row.names=FALSE,quote=FALSE)
write.csv(comparison,"contrasts.csv",row.names=FALSE,quote=FALSE)
genes <- rep(paste0("G",1:200), each=3)
feature <- paste0("T",seq_along(genes))
mu <- matrix(rep(runif(600,50,800),6),600,6)
mu[,c(1,4)] <- mu[,c(1,4)]*1.7
mu[,c(2,5)] <- mu[,c(2,5)]*.7
mu[seq(1,60,3),4:6] <- mu[seq(1,60,3),4:6]*4
counts <- matrix(rnbinom(length(mu),mu=as.vector(mu),size=40),600,6,dimnames=list(feature,samples$sample))
write.table(data.frame(gene_id=genes,feature_id=feature,counts),"counts.tsv",sep="\t",row.names=FALSE,quote=FALSE)
run <- function(script, arguments) {
    status <- system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(bin,script)),arguments))
    if(status != 0) stop(script," failed: ", status)
}
run("run_dexseq_dtu.R",c("samples.tsv","contrasts.csv","counts.tsv","0"))
fit <- readRDS("DEXSeqDataSet.treated-control.rds")
proof <- metadata(fit)$katana_design
stopifnot(isTRUE(proof$paired), grepl("subject:exon",proof$full,fixed=TRUE),grepl("subject:exon",proof$reduced,fixed=TRUE))
results <- read.delim("DEXSeqResults.treated-control.tsv")
stopifnot(nrow(results)==600, sum(is.finite(results$pvalue))>500)
# Independent native reference fit, with sample order deliberately permuted.
order <- c(6,2,4,1,5,3)
ref <- DEXSeqDataSet(counts[,order], samples[order,], design=~sample+exon+subject:exon+condition:exon,
                     featureID=feature,groupID=genes)
ref <- estimateSizeFactors(ref)
ref <- estimateDispersions(ref,formula=~sample+exon+subject:exon+condition:exon)
ref <- testForDEU(ref,fullModel=~sample+exon+subject:exon+condition:exon,reducedModel=~sample+exon+subject:exon)
expected <- as.data.frame(DEXSeqResults(ref))
stopifnot(isTRUE(all.equal(results$pvalue,expected$pvalue,tolerance=1e-5)))
# Exercise the exon-count entry with the same matched design.
dir.create("exon-counts")
gff <- unlist(lapply(1:200,function(i) {
    start <- (i-1)*1000+1
    c(paste("chr1","fixture","aggregate_gene",start,start+600,".","+",".",paste0('gene_id "G',i,'";'),sep="\t"),
      vapply(1:3,function(j)paste("chr1","fixture","exonic_part",start+(j-1)*200,start+(j-1)*200+99,".","+",".",
          paste0('gene_id "G',i,'"; transcripts "tx',i,'"; exonic_part_number "',sprintf('%03d',j),'";'),sep="\t"),""))
}))
writeLines(gff,"annotation.gff")
for(i in 1:6) write.table(data.frame(id=paste0(genes,":",sprintf("%03d",rep(1:3,200))),count=counts[,i]),
    paste0("exon-counts/",samples$sample[i],".clean.count.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
run("run_dexseq_exon.R",c("exon-counts","annotation.gff","samples.csv","contrasts.csv","0"))
exon <- readRDS("DEXSeqDataSet.treated_vs_control.rds")
stopifnot(isTRUE(metadata(exon)$katana_design$paired))
# Check the production DEG entry on actual synthetic tximport counts.
gene_counts <- rowsum(counts,genes)
txi <- list(counts=gene_counts, abundance=gene_counts+0., length=matrix(1000,nrow(gene_counts),6,dimnames=dimnames(gene_counts)),countsFromAbundance="no")
saveRDS(txi,"gene_txi.rds")
run("katana_differential.R",c("--stage","differential","--samples","samples.csv","--contrasts","contrasts.csv",
                            "--gene-rds","gene_txi.rds","--outdir","deg"))
# One missing treated observation; retain its control elsewhere in the study.
# Per-contrast exclusion must match fitting only the four complete observations.
dir.create("subset"); setwd("subset")
subset_samples <- samples[1:5, , drop=FALSE]
subset_comparison <- comparison
subset_comparison$exclude_samples <- "S3"
write.csv(subset_samples, "samples.csv", row.names=FALSE, quote=FALSE)
write.table(setNames(subset_samples,c("sample_id","condition","subject")),"samples.tsv",sep="\t",row.names=FALSE,quote=FALSE)
write.csv(subset_comparison,"contrasts.csv",row.names=FALSE,quote=FALSE)
subset_txi <- txi
for (key in c('counts','abundance','length')) subset_txi[[key]] <- txi[[key]][,1:5,drop=FALSE]
saveRDS(subset_txi,'gene_txi.rds')
write.table(data.frame(gene_id=genes,feature_id=feature,counts[,1:5]),"counts.tsv",sep="\t",row.names=FALSE,quote=FALSE)
run("run_dexseq_dtu.R",c("samples.tsv","contrasts.csv","counts.tsv","0"))
subset_fit <- readRDS("DEXSeqDataSet.treated-control.rds")
stopifnot(setequal(as.character(unique(colData(subset_fit)$katana_sample_id)), c("S1","S2","S4","S5")))
run("run_dexseq_exon.R",c("../exon-counts","../annotation.gff","samples.csv","contrasts.csv","0"))
subset_exon <- readRDS("DEXSeqDataSet.treated_vs_control.rds")
stopifnot(setequal(as.character(unique(colData(subset_exon)$katana_sample_id)), c("S1","S2","S4","S5")))
run("katana_differential.R",c("--stage","differential","--samples","samples.csv","--contrasts","contrasts.csv",
                            "--gene-rds","gene_txi.rds","--outdir","deg"))
selected_dds <- readRDS("deg/treated_vs_control.dds.rds")
stopifnot(setequal(colnames(selected_dds), c("S1","S2","S4","S5")))
write.csv(subset_samples[subset_samples$sample != 'S3',],"manual-pairs.csv",row.names=FALSE,quote=FALSE)
write.csv(comparison,"manual-contrast.csv",row.names=FALSE,quote=FALSE)
manual_txi <- txi
for (key in c('counts','abundance','length')) manual_txi[[key]] <- txi[[key]][,c(1,2,4,5),drop=FALSE]
saveRDS(manual_txi,'manual_txi.rds')
run("katana_differential.R",c("--stage","differential","--samples","manual-pairs.csv","--contrasts","manual-contrast.csv",
                            "--gene-rds","manual_txi.rds","--outdir","manual-deg"))
got <- read.delim("deg/treated_vs_control.deseq2.results.tsv")
expected <- read.delim("manual-deg/treated_vs_control.deseq2.results.tsv")
stopifnot(isTRUE(all.equal(got,expected,tolerance=1e-12)))
writeLines(c('status=passed','selected=S1,S2,S4,S5','S3_retained_in_study=true','DEG_matches_explicit_four_sample_fit=true'),
           'sample-subset.receipt.txt')
setwd('..')
writeLines(c("status=passed",paste0("DEXSeq=",packageVersion("DEXSeq")),paste0("DESeq2=",packageVersion("DESeq2")),
             "subjects=3","observations=6","real_study_data_used=false","DEXSeq_order_invariance=passed",
             "native_DTU_DEU_DEG_fits=passed"),"paired-models.receipt.txt")
