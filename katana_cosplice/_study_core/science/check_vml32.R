# Timing-only candidate: retain native scalar Brent optimization and tolerance.
root=commandArgs(TRUE)[1]
lib=dyn.load('D:/anaconda3/Library/bin/mkl_rt.2.dll')
ln=getNativeSymbolInfo('VDLN',PACKAGE=lib)
fastlog=function(x).C(ln,n=as.integer(length(x)),a=as.double(x),y=double(length(x)))$y
x=exp(seq(-27,20,length.out=2000000));difference=max(abs(log(x)-fastlog(x)));stopifnot(difference<1e-12)
s=system.time(for(i in 1:5) a=sum(log(x)))[3];f=system.time(for(i in 1:5) b=sum(fastlog(x)))[3]
writeLines(c(paste('maximum_log_error',difference),paste('reference_seconds',s),paste('vml_seconds',f),paste('sum_error',abs(a-b))),file.path(root,'evidence/revision32/development/vml_microbenchmark.txt'))
print(c(error=difference,reference=s,vml=f))
