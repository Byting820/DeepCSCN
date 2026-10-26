#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DeepCSCN))
cat("R:", R.version.string, "\n")
print(check_deepcscn_python())
cat("\nR session information:\n")
print(sessionInfo())
