# The frozen native engine and rank helper run unchanged. Only its scalar
# likelihood arithmetic is supplied by a separately audited implementation.
args=commandArgs(TRUE);root=args[1];out=args[4]
.libPaths(c(file.path(root,'r_library32'),file.path(root,'r_library'),.libPaths()))
suppressPackageStartupMessages(library(flashier))
source(file.path(out,'variance32.R'));install_variance32(root,out)
source(file.path(out,'native_engine.R'))
