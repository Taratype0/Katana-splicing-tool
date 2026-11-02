export_tom32=function(out){
 fit=readRDS(file.path(out,'native_model.rds'));blocks=list()
 for(i in seq_along(fit$blockGenes)){
  env=new.env();load(fit$TOMFiles[i],envir=env);a=as.matrix(env$TOM);diag(a)=0
  name=paste0('tom-block-',i,'.bin');con=file(file.path(out,name),'wb');writeBin(as.numeric(a),con,size=4,endian='little');close(con)
  blocks[[i]]=list(file=name,n=nrow(a),indices=fit$blockGenes[[i]]-1)
 }
 jsonlite::write_json(blocks,file.path(out,'tom_blocks.json'),auto_unbox=TRUE,pretty=TRUE)
}
if(sys.nframe()==0)export_tom32(commandArgs(TRUE)[1])
