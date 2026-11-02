# Actual native EBMF / WGCNA estimators on the identical audited input matrix.
args=commandArgs(TRUE);root=args[1];method=args[2];input=args[3];out=args[4]
.libPaths(c(file.path(root,'r_library32'),file.path(root,'r_library'),.libPaths()))
set.seed(32041);options(stringsAsFactors=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
keys=readLines(file.path(input,'event_keys.txt'));samples=read.delim(file.path(input,'samples.tsv'),check.names=FALSE)
con=file(file.path(input,'input.bin'),'rb');y=matrix(readBin(con,'double',n=length(keys)*nrow(samples),size=8,endian='little'),nrow=nrow(samples),ncol=length(keys));close(con)
stopifnot(all(is.finite(y)));colnames(y)=keys;rownames(y)=samples$sample_id
if(method %in% c('ebmf','ebmf_adaptive','ebmf_counts')){
 suppressPackageStartupMessages(library(flashier));suppressPackageStartupMessages(library(ebnm))
 rank=50
 known_sd=NULL;Y=t(y)
 if(method=='ebmf_counts'){
  con=file(file.path(input,'sampling_sd.bin'),'rb');known_sd=t(matrix(readBin(con,'double',n=length(y),size=8,endian='little'),nrow=nrow(y)));close(con)
  con=file(file.path(input,'observed_mask.bin'),'rb');obs=t(matrix(readBin(con,'integer',n=length(y),size=1,signed=FALSE),nrow=nrow(y)));close(con)
  Y[obs==0]=NA
 }
 if(method %in% c('ebmf_adaptive','ebmf_counts')){
  source(file.path(out,'rank32.R'));rank=select_rank32(y,input,out,known_sd=known_sd,paired=method=='ebmf_counts')
  if(rank==0)stop('Donor-block rank validation did not support a nonzero shared factor')
 }
 fit=flash_init(Y,S=known_sd,var_type=if(is.null(known_sd))1L else 0L)
 fit=flash_greedy(fit,Kmax=rank,ebnm_fn=list(ebnm_point_normal,ebnm_normal),maxiter=500,verbose=1)
 saveRDS(fit,file.path(out,'greedy_model.rds'))
 fit=flash_backfit(fit,maxiter=200,verbose=1)
 fit=flash_nullcheck(fit,verbose=1)
 saveRDS(fit,file.path(out,'native_model.rds'))
 write.table(data.frame(event_key=keys,fit$L_pm),file.path(out,'loadings.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(event_key=keys,fit$L_psd),file.path(out,'loading_sd.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(event_key=keys,fit$L_lfsr),file.path(out,'loading_lfsr.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(sample_id=samples$sample_id,fit$F_pm),file.path(out,'factor_scores.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(factor=seq_len(fit$n_factors),pve=fit$pve),file.path(out,'factors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 writeLines(c(paste('selected_factors',fit$n_factors),paste('selected_rank_budget',rank),paste('elbo',fit$elbo),paste('cap_only_development',method=='ebmf')),file.path(out,'fit_summary.txt'))
}else if(method=='wgcna'){
 suppressPackageStartupMessages(library(WGCNA))
 # Strong-resolution native configuration, fixed before stability/biology.
 fit=blockwiseModules(y,power=6,corType='pearson',networkType='signed',TOMType='signed',
  minModuleSize=20,mergeCutHeight=0,numericLabels=TRUE,pamRespectsDendro=FALSE,
  deepSplit=4,maxBlockSize=5000,nPreclusteringCenters=50,randomSeed=32041,verbose=2,nThreads=4,
  saveTOMs=TRUE,saveTOMFileBase=file.path(out,'TOM'))
 saveRDS(fit,file.path(out,'native_model.rds'))
 write.table(data.frame(event_key=keys,module=fit$colors,block=fit$blocks),file.path(out,'assignments.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(fit$MEs,file.path(out,'scores.tsv'),sep='\t',quote=FALSE)
 source(file.path(out,'export_tom32.R'));export_tom32(out)
}else stop('Unknown native method')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines('completed',file.path(out,'COMPLETED'))
