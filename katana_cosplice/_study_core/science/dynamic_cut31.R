# Native dynamicTreeCut tree-mode resolution diagnostic, not a new network fit.
# Examine shape only; no stability/holdout/pathway outcome chooses this adapter.
args=commandArgs(TRUE);root=args[1];output=args[2]
suppressPackageStartupMessages(library(dynamicTreeCut))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
folders=list.dirs(root,recursive=FALSE,full.names=TRUE)
if(length(args)>2)folders=folders[basename(folders)==args[3]]
for(f in folders){
 fit=readRDS(file.path(f,'model','native_model.rds'))
 keys=readLines(file.path(f,'model','event_keys.txt'));labels=rep(0,length(keys));offset=0
 for(i in seq_along(fit$blockGenes)){
  lab=cutreeDynamic(dendro=fit$dendrograms[[i]],method='tree',deepSplit=TRUE,minClusterSize=80,verbose=0)
  lab[lab>0]=lab[lab>0]+offset;labels[fit$blockGenes[[i]]]=lab;offset=max(labels)
 }
 write.table(data.frame(event_key=keys,module=labels),file.path(output,paste0(basename(f),'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 print(c(split=basename(f),modules=offset,assigned=sum(labels>0),median_size=median(table(labels[labels>0]))))
}
