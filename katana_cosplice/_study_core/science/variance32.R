# Optional arithmetic implementation; original native optimizer remains Brent.
install_variance32=function(root,adapter_dir=NULL){
 lib=dyn.load('D:/anaconda3/Library/bin/mkl_rt.2.dll')
 dll=dyn.load(if(is.null(adapter_dir))file.path(root,'runtime/variance32/variance32.dll') else file.path(adapter_dir,'variance32.dll'))
 objective=getNativeSymbolInfo('variance_objective32',PACKAGE=dll)
 ns=asNamespace('flashier');original=get('optimize.noisy',ns)
 replacement=function(R2,S2,wts=1){
  if(length(wts)!=1||length(S2)!=length(R2))return(original(R2,S2,wts))
  interval.max=max((R2-S2)/wts)
  if(interval.max<=0)return(0)
  optimize(function(x).Call(objective,R2,S2,wts,x),interval=c(0,interval.max),tol=sqrt(.Machine$double.eps))$minimum
 }
 unlockBinding('optimize.noisy',ns);assign('optimize.noisy',replacement,envir=ns);lockBinding('optimize.noisy',ns)
 invisible(list(original=original,replacement=replacement,objective=objective))
}
