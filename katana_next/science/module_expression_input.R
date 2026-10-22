args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[5],.libPaths()))
suppressPackageStartupMessages(library(DESeq2))
meta<-read.delim(args[1]);rownames(meta)<-meta$sample_id
txi<-readRDS(args[2]);stopifnot(!anyDuplicated(meta$sample_id),all(meta$sample_id %in% colnames(txi$counts)))
for(k in c('counts','abundance','length')) txi[[k]]<-txi[[k]][,meta$sample_id,drop=FALSE]
dds<-DESeqDataSetFromTximport(txi,colData=meta,design=~condition)
dds<-dds[rowSums(counts(dds)>=10)>=2,]
dds<-estimateSizeFactors(dds)
ann<-read.delim(args[3]);symbol<-setNames(ann$gene_name,ann$gene_id)
x<-log2(counts(dds,normalized=TRUE)+1)
write.table(data.frame(gene_id=rownames(x),gene_symbol=symbol[rownames(x)],x,check.names=FALSE),args[4],sep='\t',row.names=FALSE,quote=FALSE)
