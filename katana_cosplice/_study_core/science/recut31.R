# Resolution diagnostic on independently refitted native WGCNA TOM trees.
# This is a declared fixed-k tree-cut adapter, not default WGCNA module output.
args=commandArgs(TRUE); root=args[1]; output=args[2]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
folders=list.dirs(root,recursive=FALSE,full.names=TRUE)
for (f in folders) {
  p=file.path(f,'model','native_model.rds')
  if (!file.exists(p)) next
  fit=readRDS(p); keys=readLines(file.path(f,'model','event_keys.txt'))
  stopifnot(length(keys)==length(fit$colors))
  sizes=sapply(fit$blockGenes,length); allocation=pmax(1,floor(128*sizes/sum(sizes)))
  remainder=128-sum(allocation)
  if(remainder>0)allocation[order(-(128*sizes/sum(sizes)-allocation))[seq_len(remainder)]]=allocation[order(-(128*sizes/sum(sizes)-allocation))[seq_len(remainder)]]+1
  labels=rep(0,length(keys));offset=0
  for(i in seq_along(sizes)){
    lab=cutree(fit$dendrograms[[i]],k=allocation[i]);labels[fit$blockGenes[[i]]]=lab+offset;offset=offset+max(lab)
  }
  write.table(data.frame(event_key=keys,module=labels),file.path(output,paste0(basename(f),'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
}
