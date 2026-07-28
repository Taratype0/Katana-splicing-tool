options(repos = c(CRAN = "https://cloud.r-project.org"))

binary_helpers <- c("cpp11", "withr", "rstudioapi", "bslib", "tinytex", "gridExtra")
missing_helpers <- binary_helpers[!vapply(binary_helpers, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_helpers)) {
  install.packages(missing_helpers, type = "binary")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", type = "binary")
}
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  BiocManager::install("WGCNA", ask = FALSE, update = FALSE)
}
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("WGCNA installation did not complete successfully.")
}
cat("WGCNA_OK ", as.character(packageVersion("WGCNA")), "\n", sep = "")
