args=commandArgs(trailingOnly=TRUE);out=args[1]
suppressPackageStartupMessages(library(WGCNA))
x=read.csv(file.path(out,'input.csv'),row.names=1,check.names=FALSE)
a=read.delim(file.path(out,'original_assignments.tsv'),check.names=FALSE)
colors=a$module[match(colnames(x),a$event_key)]
stopifnot(!anyNA(colors))
MEs=moduleEigengenes(x,colors=colors,excludeGrey=TRUE,align='along average')$eigengenes
k=signedKME(x,MEs,corFnc='cor',corOptions="use = 'p'")
values=rep(0,ncol(x))
for (i in seq_along(values)) if(colors[i]!='grey') values[i]=k[i,paste0('kME',colors[i])]
write.table(data.frame(event_key=colnames(x),kME=values),file.path(out,'native_kme.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
