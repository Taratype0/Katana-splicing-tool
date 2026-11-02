# Engineering-only timing profile. Never used as a scientific completed fit.
args=commandArgs(TRUE);root=args[1];input=args[2];out=args[3]
.libPaths(c(file.path(root,'r_library32'),file.path(root,'r_library'),.libPaths()))
library(flashier);library(ebnm);set.seed(32041)
n=nrow(read.delim(file.path(input,'samples.tsv')));p=length(readLines(file.path(input,'event_keys.txt')))
readmat=function(name,type='double',size=8){con=file(file.path(input,name),'rb');v=readBin(con,type,n=n*p,size=size,endian='little');close(con);t(matrix(v,nrow=n))}
y=readmat('input.bin');s=readmat('sampling_sd.bin');mask=readmat('observed_mask.bin','integer',1);y[mask==0]=NA
Rprof(file.path(out,'native_profile.out'),interval=.01)
fit=flash_init(y,S=s,var_type=0L);fit=flash_greedy(fit,Kmax=1,ebnm_fn=list(ebnm_point_normal,ebnm_normal),maxiter=20,verbose=1)
Rprof(NULL);pr=summaryRprof(file.path(out,'native_profile.out'))
write.table(pr$by.total,file.path(out,'profile_by_total.tsv'),sep='\t',quote=FALSE)
write.table(pr$by.self,file.path(out,'profile_by_self.tsv'),sep='\t',quote=FALSE)
writeLines('Engineering timing only; incomplete fit is not validation evidence',file.path(out,'PROFILE_ONLY'))
