#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DeepCSCN))
set.seed(123)
n_genes <- 18L
n_cells <- 24L
true_module <- rep(seq_len(3L), each = 6L)
signal <- matrix(rnorm(3L * n_cells), nrow = 3L)
expression <- t(vapply(
  seq_len(n_genes),
  function(gene) signal[true_module[gene], ] + rnorm(n_cells, sd = 0.2),
  numeric(n_cells)
))
expression <- expression - min(expression)
rownames(expression) <- paste0("gene_", seq_len(n_genes))
colnames(expression) <- paste0("cell_", seq_len(n_cells))
output_dir <- file.path(tempdir(), "deepcscn_r_smoke_test")
result <- run_genecluster(
  expression,
  output_dir = output_dir,
  n_clusters = 3,
  epochs = 1,
  batch_size = 6,
  seed = 123,
  device = "cpu",
  patience = -1,
  overwrite = TRUE,
  verbose = FALSE
)
stopifnot(
  identical(dim(result$features), c(n_genes, 512L)),
  nrow(result$clusters) == n_genes,
  file.exists(result$files$metadata),
  file.exists(result$files$model)
)
print(result)
cat("R-to-Python smoke test passed.\n")
