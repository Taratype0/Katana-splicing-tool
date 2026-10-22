# Full fixed-module metrics. Export thresholds never truncate kWithin.
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(TRUE)
plan <- fromJSON(args[1], simplifyVector=FALSE)
x <- t(as.matrix(read.delim(plan$matrix, row.names=1, check.names=FALSE)))
events <- read.delim(plan$events, check.names=FALSE, stringsAsFactors=FALSE, na.strings="", colClasses='character')
rownames(events) <- events$event_key
capture.output(sessionInfo(), file=plan$session)
corfun <- function(a,b=NULL) {
  if (plan$method == 'bicor') bicor(a,b,use='pairwise.complete.obs',quick=0,pearsonFallback='individual')
  else stats::cor(a,b,use='pairwise.complete.obs',method=plan$method)
}
metrics <- function(z,units,genes) {
  p <- ncol(z); n <- nrow(z)
  kw <- cross <- rep(0,p); eligible <- crosseligible <- degree <- crossdegree <- rep(0,p)
  known <- !is.na(genes) & nzchar(genes)
  # Exact biological-unit counts, including repeated measures within units.
  groups <- split(seq_len(n),units)
  for (start in seq.int(1,p,by=as.integer(plan$block_size))) {
    ix <- start:min(p,start+as.integer(plan$block_size)-1L)
    r <- suppressWarnings(corfun(z[,ix,drop=FALSE],z))
    r <- matrix(r,nrow=length(ix),ncol=p)
    valid <- is.finite(z); counts <- matrix(0,length(ix),p)
    if (length(groups)==n) counts <- t(valid[,ix,drop=FALSE]) %*% valid
    else for (g in groups) counts <- counts + (t(valid[g,ix,drop=FALSE]) %*% valid[g,,drop=FALSE] > 0)
    r[counts < as.integer(plan$minimum_samples)] <- NA_real_
    r[] <- pmax(-1,pmin(1,r))
    a <- if (plan$signed) ((1+r)/2)^plan$power else abs(r)^plan$power
    a[cbind(seq_along(ix),ix)] <- NA_real_
    cg <- outer(genes[ix],genes,'!=') & outer(known[ix],known,'&')
    cg[is.na(cg)] <- FALSE
    kw[ix] <- rowSums(a,na.rm=TRUE)
    eligible[ix] <- rowSums(is.finite(a))
    degree[ix] <- rowSums(a >= plan$edge_threshold,na.rm=TRUE)
    if (b==1L) for (j in seq_along(ix)) {
      targets <- which(is.finite(a[j,]) & a[j,]>=plan$edge_threshold)
      if (length(targets)) {
        targets <- targets[order(-a[j,targets])][seq_len(min(8L,length(targets)))]
        display <- data.frame(module=module$name,condition=scope$name,event_a=colnames(z)[ix[j]],
          event_b=colnames(z)[targets],adjacency=a[j,targets],correlation=r[j,targets],pair_n=counts[j,targets])
        write.table(display,edgeout,sep='\t',quote=TRUE,row.names=FALSE,col.names=FALSE,na='NA')
      }
    }
    a[!cg] <- NA_real_
    cross[ix] <- rowSums(a,na.rm=TRUE)
    crosseligible[ix] <- rowSums(is.finite(a))
    crossdegree[ix] <- rowSums(a >= plan$edge_threshold,na.rm=TRUE)
  }
  kw[eligible==0] <- NA_real_; cross[crosseligible==0] <- NA_real_
  # Mean imputation is used ONLY for PC1, refitted independently in each bootstrap.
  # Adjacency always uses observed, pairwise-complete PSI.
  means <- colMeans(z,na.rm=TRUE); sdv <- apply(z,2,sd,na.rm=TRUE)
  good <- is.finite(means) & is.finite(sdv) & sdv > 1e-12
  kme <- rep(NA_real_,p)
  if (sum(good)>=2) {
    complete <- z[,good,drop=FALSE]
    for (i in seq_len(ncol(complete))) complete[!is.finite(complete[,i]),i] <- means[good][i]
    me <- tryCatch(moduleEigengenes(complete,rep('blue',ncol(complete)),impute=FALSE,subHubs=FALSE,
                                   trapErrors=TRUE,align='along average')$eigengenes[,1],error=function(e) NULL)
    if (!is.null(me) && all(is.finite(me)) && sd(me)>1e-12) {
      kme <- as.numeric(suppressWarnings(corfun(z,matrix(me,ncol=1))))
      kme[colSums(is.finite(z)) < plan$minimum_samples] <- NA_real_
    }
  }
  cross_possible <- vapply(seq_len(p),function(i) if(known[i]) sum(known & genes!=genes[i],na.rm=TRUE) else 0L,integer(1))
  data.frame(event_key=colnames(z),kWithin=kw,kWithin_no_same_gene=cross,kME=kme,
    valid_pairs=eligible,pair_coverage=eligible/max(1,p-1),supported_degree=degree,
    cross_gene_valid_pairs=crosseligible,cross_gene_possible_pairs=cross_possible,
    cross_gene_pair_coverage=ifelse(cross_possible>0,crosseligible/cross_possible,NA_real_),
    cross_gene_supported_degree=crossdegree,observed_samples=colSums(is.finite(z)),
    gene_known=known,gene_annotation_fraction=mean(known),stringsAsFactors=FALSE)
}
out <- file(plan$output,'wt'); first <- TRUE
edgeout <- file(plan$edges,'wt')
writeLines('module\tcondition\tevent_a\tevent_b\tadjacency\tcorrelation\tpair_n',edgeout)
for (scope in plan$scopes) for (module in plan$modules) {
  keys <- unlist(module$events); genes <- events[keys,'host_gene']
  reps <- c(list(list(indices=scope$indices,units=scope$units)),scope$bootstraps)
  for (b in seq_along(reps)) {
    item <- reps[[b]]; indices <- as.integer(unlist(item$indices))+1L; units <- unlist(item$units)
    z <- x[indices,keys,drop=FALSE]
    if (length(unique(units)) < plan$minimum_samples || length(keys)<3) {
      result <- data.frame(event_key=keys,kWithin=NA_real_,kWithin_no_same_gene=NA_real_,kME=NA_real_,
        valid_pairs=0,pair_coverage=0,supported_degree=0,cross_gene_valid_pairs=0,cross_gene_possible_pairs=0,
        cross_gene_pair_coverage=NA_real_,cross_gene_supported_degree=0,observed_samples=colSums(is.finite(z)),
        gene_known=!is.na(genes)&nzchar(genes),gene_annotation_fraction=mean(!is.na(genes)&nzchar(genes)))
    } else result <- metrics(z,units,genes)
    result$module <- module$name; result$condition <- scope$name; result$bootstrap <- b-1L
    result$biological_units <- length(unique(units)); result$module_events <- length(keys)
    result$unique_source_biological_units <- length(unique(unlist(plan$source_units)[indices]))
    write.table(result,out,sep='\t',quote=TRUE,row.names=FALSE,col.names=first,na='NA')
    first <- FALSE
  }
}
close(out)
close(edgeout)
