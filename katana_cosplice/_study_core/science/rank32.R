# Donor-by-gene bi-cross-validation selects rank before the final full-data fit.
# This is a CoSplice selection wrapper around unmodified flashier estimators.
select_rank32=function(y,input,out,known_sd=NULL,paired=FALSE){
 genes=read.delim(file.path(input,'events.tsv'))$ensembl_gene_id
 samples=read.delim(file.path(input,'samples.tsv'),colClasses=c(biological_sample_id='character'))
 nr=ncol(y);nc=nrow(y);Y=t(y)
 readmask=function(name){con=file(file.path(input,name),'rb');m=t(matrix(readBin(con,'integer',n=nr*nc,size=1,signed=FALSE),nrow=nc,ncol=nr));close(con);m>0}
 mask=readmask('inner_mask.bin');observed=readmask('observed_mask.bin');train=Y;train[mask]=NA
 if(!is.null(known_sd))train[!observed]=NA
 unit=split(seq_len(nc),samples$biological_sample_id);gene_groups=match(genes,unique(genes))
 evaluate=function(pred){
  err=(Y-pred)^2;err[!mask|!observed]=NA
  v=sapply(unit,function(ix){
   a=rowMeans(err[,ix,drop=FALSE],na.rm=TRUE);valid=is.finite(a)
   sums=rowsum(replace(a,!valid,0),gene_groups,reorder=FALSE)
   n=rowsum(as.numeric(valid),gene_groups,reorder=FALSE)
   mean(sums[n>0]/n[n>0])
  })
  c(mean=mean(v),se=sd(v)/sqrt(length(v)),v)
 }
 cv=flash_init(train,S=known_sd,var_type=if(is.null(known_sd))1L else 0L);record=list(c(rank=0,evaluate(matrix(0,nr,nc))))
 best=record[[1]]['mean'];stale=0;cap=min(64,nc-length(unit)-1)
 for(k in seq_len(cap)){
  before=cv$n_factors
  cv=flash_greedy(cv,Kmax=1,ebnm_fn=list(ebnm_point_normal,ebnm_normal),maxiter=500,verbose=0)
  if(cv$n_factors<=before)break
  risk=evaluate(tcrossprod(cv$L_pm,cv$F_pm));record[[length(record)+1]]=c(rank=cv$n_factors,risk)
  message('Inner donor-block rank ',cv$n_factors,': mean=',round(risk['mean'],6),' SE=',round(risk['se'],6))
  if(risk['mean']<best*(1-1e-4)){best=risk['mean'];stale=0}else stale=stale+1
  if(stale>=5)break
 }
 path=as.data.frame(do.call(rbind,record));imin=which.min(path$mean);threshold=path$mean[imin]+path$se[imin]
 path$within_unpaired_one_se=path$mean<=threshold
 losses=as.matrix(path[,4:(3+length(unit))]);difference=sweep(losses,2,losses[imin,],'-')
 path$paired_difference=rowMeans(difference);path$paired_se=apply(difference,1,sd)/sqrt(length(unit))
 path$within_paired_one_se=path$paired_difference<=path$paired_se+1e-12
 path$within_one_se=if(paired)path$within_paired_one_se else path$within_unpaired_one_se
 selected=min(path$rank[path$within_one_se]);path$selected=path$rank==selected
 write.table(path,file.path(out,'rank_selection.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 saveRDS(cv,file.path(out,'rank_selection_model.rds'))
 writeLines(c(paste('selected_rank',selected),paste('searched_through_rank',max(path$rank)),paste('safety_cap_reached',max(path$rank)>=cap),
  paste('paired_donor_standard_error',paired),
  'Rule: smallest rank within one donor-level standard error of the best internal gene-by-donor heldout loss; paired differences when paired=TRUE'),file.path(out,'rank_selection.txt'))
 selected
}
