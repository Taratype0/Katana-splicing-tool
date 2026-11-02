args <- commandArgs(TRUE)
root <- args[1]; input <- args[2]; assignments <- args[3]; out <- args[4]
suppressPackageStartupMessages(library(WGCNA))
set.seed(31)
dat <- t(as.matrix(read.delim(input,row.names=1,check.names=FALSE)))
tab <- read.delim(assignments,check.names=FALSE)
colors <- tab$module[match(colnames(dat),tab$event_key)]
stopifnot(!anyNA(colors))
fit <- mergeCloseModules(dat,colors,cutHeight=.25,iterate=TRUE,unassdColor=0,verbose=0)
write.table(data.frame(event_key=colnames(dat),module=fit$colors,native_module=colors),file.path(out,'assignments.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
saveRDS(fit,file.path(out,'native_merge.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
