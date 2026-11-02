# Native count-aware flashier estimator with a separately audited rank selector.
args=commandArgs(TRUE);root=args[1];method=args[2];input=args[3];out=args[4]
.libPaths(c(file.path(root,'r_library32'),file.path(root,'r_library'),.libPaths()))
stopifnot(method=='ebmf_counts');set.seed(32041);options(stringsAsFactors=FALSE)
keys=readLines(file.path(input,'event_keys.txt'));samples=read.delim(file.path(input,'samples.tsv'),check.names=FALSE)
con=file(file.path(input,'input.bin'),'rb');y=matrix(readBin(con,'double',n=length(keys)*nrow(samples),size=8,endian='little'),nrow=nrow(samples),ncol=length(keys));close(con)
stopifnot(all(is.finite(y)));colnames(y)=keys;rownames(y)=samples$sample_id
suppressPackageStartupMessages(library(flashier));suppressPackageStartupMessages(library(ebnm))
con=file(file.path(input,'sampling_sd.bin'),'rb');known_sd=t(matrix(readBin(con,'double',n=length(y),size=8,endian='little'),nrow=nrow(y)));close(con)
con=file(file.path(input,'observed_mask.bin'),'rb');obs=t(matrix(readBin(con,'integer',n=length(y),size=1,signed=FALSE),nrow=nrow(y)));close(con)
Y=t(y);Y[obs==0]=NA
path=read.delim(file.path(out,'parallel_rank.tsv'));rank=sum(as.character(path$retained)%in%c('TRUE','True','true'))
path$selected=path$rank==rank;write.table(path,file.path(out,'rank_selection.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('selected_rank',rank),paste('searched_through_rank',max(path$rank)),paste('safety_cap_reached',rank==max(path$rank)),
 'Rule: consecutive eigenvalues above 95th-percentile gene-by-donor block-sign null; 999 randomizations; no fixed requested factor count',
 'Borrowed from FlipPA; the structured blocks are a CoSplice adaptation without claiming the original independent-entry theorem'),file.path(out,'rank_selection.txt'))
if(rank==0)stop('Block-sign rank selection did not support a shared cross-gene program')
fit=flash_init(Y,S=known_sd,var_type=0L)
fit=flash_greedy(fit,Kmax=rank,ebnm_fn=list(ebnm_point_normal,ebnm_normal),maxiter=500,verbose=1)
saveRDS(fit,file.path(out,'greedy_model.rds'))
fit=flash_backfit(fit,maxiter=200,verbose=1)
fit=flash_nullcheck(fit,verbose=1)
if(fit$n_factors==0)stop('Native null-check retained no supported shared factor')
saveRDS(fit,file.path(out,'native_model.rds'))
write.table(data.frame(event_key=keys,fit$L_pm),file.path(out,'loadings.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(event_key=keys,fit$L_psd),file.path(out,'loading_sd.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(event_key=keys,fit$L_lfsr),file.path(out,'loading_lfsr.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(sample_id=samples$sample_id,fit$F_pm),file.path(out,'factor_scores.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(factor=seq_len(fit$n_factors),pve=fit$pve),file.path(out,'factors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('selected_factors',fit$n_factors),paste('selected_rank_budget',rank),paste('elbo',fit$elbo),'cap_only_development FALSE'),file.path(out,'fit_summary.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines('completed',file.path(out,'COMPLETED'))
