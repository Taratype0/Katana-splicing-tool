args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(7L,9L)) stop('design, tximport RDS, gene names TSV, GMT, output, library, seed; optional downstream rules and policy')
.libPaths(c(args[6],.libPaths()))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(fgsea))
set.seed(as.integer(args[7]))
meta <- read.delim(args[1],check.names=FALSE)
if(anyDuplicated(meta$sample_id)||anyNA(meta[,c('sample_id','factor_a','factor_b')])) stop('Invalid design')
meta$factor_a <- factor(meta$factor_a,levels=c('reference','target'))
meta$factor_b <- factor(meta$factor_b,levels=c('reference','target'))
if(anyNA(meta$factor_a)||anyNA(meta$factor_b)||any(table(meta$factor_a,meta$factor_b)<2)) stop('Need replicated 2x2 design')
rownames(meta)<-meta$sample_id
txi<-readRDS(args[2]); stopifnot(all(c('counts','abundance','length','countsFromAbundance') %in% names(txi)))
for(n in c('counts','abundance','length')) {
 if(anyDuplicated(rownames(txi[[n]]))||anyDuplicated(colnames(txi[[n]]))||!all(meta$sample_id %in% colnames(txi[[n]]))) stop('Invalid matrix identities')
 txi[[n]]<-txi[[n]][,meta$sample_id,drop=FALSE]
}
if(!identical(dimnames(txi$counts),dimnames(txi$length))) stop('Counts and effective lengths differ')
full<-~factor_a*factor_b
if('batch' %in% names(meta) && length(unique(meta$batch))>1) {meta$batch<-factor(meta$batch);full<-~batch+factor_a*factor_b}
mm<-model.matrix(full,meta)
if(qr(mm)$rank<ncol(mm)) stop('Confounded design')
dds<-DESeqDataSetFromTximport(txi,colData=meta,design=full)
keep<-rowSums(counts(dds)>=10)>=2
out<-args[5];dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.table(data.frame(gene_id=rownames(dds),eligible=keep),file.path(out,'gene_filter_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
dds<-DESeq(dds[keep,],quiet=FALSE)
cn<-resultsNames(dds); ia<-grep('factor_a.*factor_b',cn,value=TRUE);stopifnot(length(ia)==1)
a<-grep('^factor_a_',cn,value=TRUE);b<-grep('^factor_b_',cn,value=TRUE)
stopifnot(length(a)==1,length(b)==1)
contrasts<-list(A_at_Breference=a,A_at_Btarget=c(a,ia),B_at_Areference=b,B_at_Atarget=c(b,ia),interaction=ia)
ann<-read.delim(args[3],check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c('gene_id','gene_name') %in% names(ann)),!anyDuplicated(ann$gene_id))
symbols<-setNames(ann$gene_name,ann$gene_id)
scope_keep <- function(gene_ids) rep(TRUE,length(gene_ids))
if(length(args)==9L) {
 rules<-read.delim(args[8],stringsAsFactors=FALSE)
 canonical<-function(x) tolower(sub('\\.[0-9]+$','',as.character(x)))
 allowed<-canonical(rules$identifier[rules$action=='keep'])
 blocked<-canonical(rules$identifier[rules$action=='exclude'])
 scope_keep<-function(gene_ids) {
  key<-canonical(gene_ids);symbol<-canonical(symbols[gene_ids])
  keep<-!key %in% blocked & !symbol %in% blocked
  if(args[9]=='strict') keep<-keep & key %in% allowed
  keep
 }
}
paths<-gmtPathways(args[4]); allgsea<-list()
for(n in names(contrasts)) {
 cc<-contrasts[[n]]
 res<-if(length(cc)==1) results(dds,name=cc,independentFiltering=FALSE) else results(dds,contrast=list(cc),independentFiltering=FALSE)
 tab<-as.data.frame(res);tab$gene_id<-rownames(tab);tab$gene_symbol<-symbols[tab$gene_id]
 if(length(args)==9L) {
  write.table(tab,file.path(out,paste0('unfiltered_model_',n,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  tab<-tab[scope_keep(tab$gene_id),,drop=FALSE]
 }
 write.table(tab,file.path(out,paste0(n,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 valid<-is.finite(tab$stat)&!is.na(tab$gene_symbol)&nzchar(tab$gene_symbol)
 stats<-tapply(tab$stat[valid],tab$gene_symbol[valid],mean)
 stats<-sort(stats,decreasing=TRUE)
 if(length(stats)<15) stop('Too few mapped ranked genes')
 gs<-fgseaMultilevel(paths,stats,minSize=15,maxSize=500,eps=0,nproc=1)
 gs$contrast<-n; gs$leadingEdge<-vapply(gs$leadingEdge,paste,character(1),collapse=';')
 allgsea[[n]]<-as.data.frame(gs)
}
gs<-do.call(rbind,allgsea); gs$study_family_FDR<-p.adjust(gs$pval,'BH')
write.table(gs,file.path(out,'ranked_pathways.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
norm<-log2(counts(dds,normalized=TRUE)+1)
if(length(args)==9L) norm<-norm[scope_keep(rownames(norm)),,drop=FALSE]
write.table(data.frame(gene_id=rownames(norm),gene_symbol=symbols[rownames(norm)],norm,check.names=FALSE),file.path(out,'log_normalized_expression.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(dds,file.path(out,'expression_model.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
writeLines(c(paste('design:',deparse(full)),paste('countsFromAbundance:',txi$countsFromAbundance),
 'Official DESeqDataSetFromTximport: effective-length handling retained.',
 'Eligibility >=10 estimated counts in >=2 samples; BH gene FDR per contrast.',
 'fgsea: signed Wald statistic; duplicate symbols averaged; 15-500 mapped genes; BH study_family_FDR across five contrasts and all pathways.',
 'Enrichment is competitive association, not independent validation or proof of pathway activity.'),file.path(out,'METHODS.txt'))
cat('Expression interaction and ranked pathways complete\n')
