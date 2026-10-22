args<-commandArgs(TRUE);sourceRoot<-normalizePath(args[1]);input<-normalizePath(args[2]);out<-normalizePath(args[3]);setwd(out)
suppressPackageStartupMessages({library(WGCNA);library(limma);library(dplyr);library(magrittr);library(tidyr);library(jsonlite)})
for(name in c('goodSamplesJunct.R','JCNA1pass.R','JCNAgenePrepare.R','JCNA2pass.R'))source(file.path(sourceRoot,'R',name),keep.source=TRUE)
options(stringsAsFactors=FALSE,device=function(...)pdf(file=tempfile('DJExpress_',tmpdir=out,fileext='.pdf'),...))
readmat<-function(name)as.matrix(read.delim(file.path(input,name),row.names=1,check.names=FALSE))
counts<-readmat('counts.tsv');traits<-as.data.frame(readmat('traits.tsv'))
norm<-limma::voom(counts,lib.size=colSums(counts),normalize.method='none',plot=FALSE)
saveRDS(norm,'voom_preparation.rds');x<-as.data.frame(t(norm$E))
traits<-traits[rownames(x),,drop=FALSE]
sft<-pickSoftThreshold(x,powerVector=c(1:10,seq(12,20,2)),verbose=0)
saveRDS(sft,'soft_threshold.rds')
write.table(t(x),'normalized_junction_expression.tsv',sep='\t',quote=FALSE,col.names=NA)
ledger<-list();status<-function(){write_json(list(engine='DJExpress',scope='unchanged author core; Katana prepared-object adapter',states=ledger),'METHOD_STATUS.json',auto_unbox=TRUE,pretty=TRUE)}
run<-function(name,fun){
 ledger[[name]]<<-list(status='running');status()
 result<-tryCatch(fun(),error=function(e){ledger[[name]]<<-list(status='failed',error=conditionMessage(e));status();stop(e)})
 saveRDS(result,paste0(name,'.rds'));ledger[[name]]<<-list(status='complete');status();result
}
export<-function(result,prefix){
 write.table(data.frame(event_key=colnames(result$datExpr),module=as.character(result$net$colors)),paste0(prefix,'_assignments.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(result$net$MEs,paste0(prefix,'_scores.tsv'),sep='\t',quote=FALSE,col.names=NA)
 write.table(result$moduleTraitCor,paste0(prefix,'_module_trait_cor.tsv'),sep='\t',quote=FALSE,col.names=NA)
 write.table(result$moduleTraitPvalue,paste0(prefix,'_module_trait_p.tsv'),sep='\t',quote=FALSE,col.names=NA)
}
pdf('JCNA1pass.pdf',width=12,height=8)
first<-run('JCNA1pass',function()JCNA1pass(list(datExpr=x,datTraits=traits,sft=sft),cor.method='bicor'))
# Upstream genePrepare reads J1pass$sig.cors.trait from its enclosing session.
# Bind the actual original fit under the author's documented example name.
J1pass<-first
export(first,'JCNA1pass')
if(file.exists(file.path(input,'gene_expression.tsv'))){
 genes<-readmat('gene_expression.tsv')
 pdf('JCNAgenePrepare.pdf',width=12,height=8)
 gp<-run('JCNAgenePrepare',function()JCNAgenePrepare(first,genExpr=as.data.frame(genes),gtf=file.path(input,'reference.gtf.gz'),cor.method='bicor'))
 traitIndex<-data.frame(trait=names(traits),stage=sprintf('JCNA2pass_%03d',seq_along(traits)))
 write.table(traitIndex,'SECOND_PASS_INDEX.tsv',sep='\t',quote=FALSE,row.names=FALSE)
 for(i in seq_along(traits)){
  name<-traitIndex$stage[i];pdf(paste0(name,'.pdf'),width=12,height=8)
  result<-run(name,function()JCNA2pass(first,gp,trait=names(traits)[i],workDir=out,cor.method='bicor'))
  if(!is.null(result$net)){
   keep<-rownames(gp$JunctGeneTrait[[names(traits)[i]]])[gp$JunctGeneTrait[[names(traits)[i]]]$gene.association=='gene_independent'];keep<-keep[!is.na(keep)]
   write.table(data.frame(event_key=keep,module=as.character(result$net$colors)),paste0(name,'_assignments.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
   write.table(result$net$MEs,paste0(name,'_scores.tsv'),sep='\t',quote=FALSE,col.names=NA)
  }else{ledger[[name]]<-list(status='complete_no_eligible_network');status()}
 }
}
while(dev.cur()>1)dev.off()
capture.output(sessionInfo(),file='R_SESSION.txt')
assigned<-sum(first$net$colors!=0)
write_json(list(status='complete',network_status=if(assigned>0)'modules_available' else 'no_estimable_modules',assigned_features=assigned,retained_features=ncol(first$datExpr),engine='DJExpress',primary='JCNA1pass',second_pass='separate per trait, no automatic winner or union',native_raw_JCNAprepare='not called; known undefined min.expressed in fixed author revision'),'ENGINE_COMPLETE.json',auto_unbox=TRUE,pretty=TRUE)
