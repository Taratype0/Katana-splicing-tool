args <- commandArgs(TRUE)
root <- normalizePath(args[1],winslash='/'); method <- args[2]; input <- args[3]; out <- args[4]
threads <- as.integer(Sys.getenv('KATANA_CPUS', '2'))
if (is.na(threads) || threads < 1L) stop('KATANA_CPUS must be positive')
set.seed(31)
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dat <- t(as.matrix(read.delim(input,row.names=1,check.names=FALSE)))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
if(method %in% c('wgcna','wgcna_single','wgcna_enhanced','wgcna_robust')) {
  block <- if(method=='wgcna_single') ncol(dat) else if(method %in% c('wgcna_enhanced','wgcna_robust')) 5000L else 2000L
  if(method=='wgcna_enhanced') dat <- apply(dat,2,rank,ties.method='average')
  # No global correlation matrix is calculated before blockwiseModules.
  net <- blockwiseModules(dat,power=6,corType=if(method=='wgcna_enhanced') 'pearson' else 'bicor',networkType='signed',TOMType='signed',
    minModuleSize=as.integer(Sys.getenv('KATANA_MIN_MODULE_SIZE','20')),mergeCutHeight=as.numeric(Sys.getenv('KATANA_MERGE_CUT_HEIGHT','.25')),numericLabels=TRUE,pamRespectsDendro=FALSE,
    deepSplit=2,maxBlockSize=block,randomSeed=31,verbose=0,nThreads=threads,
    maxPOutliers=if(method=='wgcna_robust') .1 else 1)
  km <- signedKME(dat,net$MEs,corFnc='cor',corOptions="use = 'p'")
  own <- rep(NA_real_,ncol(dat))
  for(i in seq_along(own)) if(net$colors[i]!=0) own[i] <- km[i,paste0('kME',net$colors[i])]
  write.table(data.frame(event_key=colnames(dat),module=net$colors,block=net$blocks,kME=own),
    file.path(out,'assignments.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(net$MEs,file.path(out,'scores.tsv'),sep='\t',quote=FALSE)
  saveRDS(net,file.path(out,'native_model.rds'))
} else if(method=='megena') {
  suppressPackageStartupMessages(library(MEGENA))
  suppressPackageStartupMessages(library(igraph))
  el <- calculate.correlation(t(dat),doPerm=100,doPar=FALSE,method='pearson',
    FDR.cutoff=.05,is.signed=TRUE)
  write.table(el,file.path(out,'candidate_edges.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  if(nrow(el)<3) stop('Fewer than three eligible MEGENA edges')
  pfn <- calculate.PFN(el[,1:3],doPar=FALSE)
  write.table(pfn,file.path(out,'pfn.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  g <- graph_from_data_frame(pfn,directed=FALSE)
  fit <- do.MEGENA(g,doPar=FALSE,n.perm=100,remove.unsig=FALSE,min.size=20)
  saveRDS(fit,file.path(out,'native_model.rds'))
  # Retain all hierarchical memberships, including parents, without choosing a
  # level using the test data. Modules are descriptive when remove.unsig=FALSE.
  memb <- fit$module.output$modules
  tab <- do.call(rbind,lapply(names(memb),function(m) data.frame(module=m,event_key=as.character(memb[[m]]))))
  write.table(tab,file.path(out,'memberships.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
} else if(method %in% c('joint','independent')) {
  suppressPackageStartupMessages(library(JGL))
  meta <- read.delim(args[5],check.names=FALSE)
  meta <- meta[match(rownames(dat),meta$sample_id),]
  eligible <- names(which(table(meta$condition)>=15))
  if(length(eligible)<2) stop('Fewer than two conditions have 15 training samples')
  centers <- colMeans(dat); scales <- apply(dat,2,sd);scales[scales<1e-8]<-1
  z <- sweep(sweep(dat,2,centers,'-'),2,scales,'/')
  Y <- lapply(eligible,function(k)z[meta$condition==k,,drop=FALSE]); names(Y)<-eligible
  # Trace diagnostics only; original author solver body and iterates unchanged.
  solver_audit <- list()
  trace('admm.iters',where=asNamespace('JGL'),print=FALSE,exit=quote({
    .GlobalEnv$solver_audit[[length(.GlobalEnv$solver_audit)+1]] <-
      data.frame(iterations=iter,relative_change=diff_value,primal_l1=diff)
  }))
  fit <- JGL(Y,penalty='fused',lambda1=.15,lambda2=if(method=='joint').1 else 0,
    return.whole.theta=TRUE,maxiter=500,tol=1e-5)
  if(length(solver_audit)) write.table(do.call(rbind,solver_audit),file.path(out,'solver_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  saveRDS(fit,file.path(out,'native_model.rds'))
  write.table(data.frame(event_key=colnames(dat),center=centers,scale=scales),
    file.path(out,'transform.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  for(k in seq_along(eligible)) {
    write.table(fit$theta[[k]],file.path(out,paste0('precision_',k,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  }
  write.table(data.frame(index=seq_along(eligible),condition=eligible,n=sapply(Y,nrow)),
    file.path(out,'conditions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
} else stop('Unknown method')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
