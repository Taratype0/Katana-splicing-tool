# Read native optimizer flags without updating or refitting any saved model.
args=commandArgs(TRUE);root=args[1];out=args[2]
.libPaths(c(file.path(root,'r_library32'),file.path(root,'r_library'),.libPaths()))
base=file.path(root,'evidence/revision32/formal/ebmf_counts')
splits=c('full',sprintf('refit_%02d',0:19),sprintf('fold_%02d',0:4))
rows=list()
for(s in splits){
 p=file.path(base,s);stopifnot(file.exists(file.path(p,'COMPLETED')))
 f=readRDS(file.path(p,'native_model.rds'));log=readLines(file.path(p,'engine.log'),warn=FALSE)
 rows[[s]]=data.frame(split=s,retained_factors=f$n_factors,elbo=f$elbo,
  backfit_iteration_limit_reached=isTRUE(f$flash_fit[['maxiter.reached']]),
  native_timeout_reached=isTRUE(f$flash_fit[['timeout.reached']]),
  iteration_limit_messages=sum(grepl('Maximum number of iterations reached',log,fixed=TRUE)),
  convergence_tolerance=f$flash_fit[['conv.tol']])
 rm(f);invisible(gc())
}
write.table(do.call(rbind,rows),out,sep='\t',quote=FALSE,row.names=FALSE)
