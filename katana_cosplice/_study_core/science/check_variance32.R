args=commandArgs(TRUE);root=args[1]
.libPaths(c(file.path(root,'r_library32'),file.path(root,'r_library'),.libPaths()))
library(flashier);source(file.path(root,'science/variance32.R'));functions=install_variance32(root);set.seed(32317)
rows=list()
for(truth in c(0,.2,5)){
 n=1000000;s=exp(runif(n,-8,2));r=(s+truth)*rchisq(n,1)
 a=system.time(reference<-functions$original(r,s))[3]
 b=system.time(candidate<-functions$replacement(r,s))[3]
 original_objective=function(x)sum(log(x+s))+sum(r/(x+s))
 objective_error=abs(original_objective(reference)-.Call(functions$objective,r,s,1,reference))
 stopifnot(abs(reference-candidate)<1e-5,objective_error<1e-7)
 rows[[length(rows)+1]]=data.frame(true_extra_variance=truth,reference=reference,candidate=candidate,objective_error=objective_error,reference_seconds=a,candidate_seconds=b)
}
tab=do.call(rbind,rows);print(tab);write.table(tab,file.path(root,'evidence/revision32/acceptance/variance_arithmetic.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
