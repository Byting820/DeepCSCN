#' Select the Python interpreter used by DeepCSCN
#'
#' Configure the Python interpreter or Conda environment before the first
#' Python call in an R session. If neither argument is supplied, reticulate's
#' normal interpreter discovery is used. The project Conda environments are
#' named `deepcscn` by default.
#'
#' @param python Optional path to a Python executable.
#' @param condaenv Optional Conda environment name or path.
#' @param required Whether failure to locate the requested interpreter should
#'   raise an error.
#'
#' @return Invisibly returns the selected Python configuration.
#' @export
configure_deepcscn_python <- function(python = NULL, condaenv = NULL, required = TRUE) {
  if (!is.null(python) && !is.null(condaenv)) {
    stop("Supply only one of `python` and `condaenv`.", call. = FALSE)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The reticulate package is required. Install it with install.packages('reticulate').", call. = FALSE)
  }
  if (!is.null(python)) {
    python <- normalizePath(path.expand(python), mustWork = TRUE)
    reticulate::use_python(python, required = required)
  } else if (!is.null(condaenv)) {
    reticulate::use_condaenv(condaenv, required = required)
  }
  invisible(reticulate::py_config())
}

.deepcscn_python_dir <- function() {
  directory <- system.file("python", package = "DeepCSCN")
  if (!nzchar(directory)) {
    development_path <- file.path(getwd(), "inst", "python")
    if (dir.exists(development_path)) {
      directory <- normalizePath(development_path, mustWork = TRUE)
    }
  }
  if (!nzchar(directory) || !dir.exists(directory)) {
    stop("Could not locate the Python backend bundled with DeepCSCN.", call. = FALSE)
  }
  directory
}

.deepcscn_backend <- function() {
  reticulate::import_from_path(
    module = "deepcscn_backend",
    path = .deepcscn_python_dir(),
    convert = TRUE,
    delay_load = FALSE
  )
}

#' Check the DeepCSCN Python environment
#'
#' Imports the bundled GeneCluster backend and reports the exact Python,
#' PyTorch, NumPy, pandas, and scikit-learn versions. This is a quick setup
#' check and does not train a model.
#'
#' @param python Optional path to a Python executable.
#' @param condaenv Optional Conda environment name or path.
#'
#' @return A named list of dependency versions and CUDA availability.
#' @export
check_deepcscn_python <- function(python = NULL, condaenv = NULL) {
  configure_deepcscn_python(python = python, condaenv = condaenv)
  versions <- .deepcscn_backend()$dependency_versions()
  class(versions) <- c("deepcscn_python_check", class(versions))
  versions
}

#' Train GeneCluster from R through the bundled Python backend
#'
#' `run_genecluster()` is the R entry point for DeepCSCN's CNN-SPP and
#' DeepCluster feature-learning stage. Genes are observations and columns are
#' ordered cells. The supplied values are used directly; this function does
#' not perform CP10K or log1p normalization.
#'
#' The Python backend writes the feature matrix, final pseudo-clusters,
#' training history, PyTorch checkpoint, and a JSON audit record containing
#' parameters, software versions, device information, and an SHA-256 checksum
#' of the exact matrix used for training.
#'
#' @param expression A numeric gene-by-cell matrix/data frame, or a path to a
#'   delimited file whose first column contains unique gene names.
#' @param output_dir Directory in which the canonical input and model outputs
#'   are written.
#' @param n_clusters Number of DeepCluster pseudo-clusters.
#' @param epochs Maximum training epochs.
#' @param batch_size Mini-batch size.
#' @param learning_rate SGD learning rate.
#' @param momentum SGD momentum.
#' @param weight_decay L2 weight decay.
#' @param focal_gamma Focal-loss gamma; use zero for cross-entropy behavior.
#' @param seed Random seed used by R, NumPy, scikit-learn, and PyTorch.
#' @param device One of `"auto"`, `"cpu"`, `"cuda"`, or `"cuda:<index>"`.
#' @param num_workers Number of PyTorch data-loader workers. Zero is the most
#'   portable and reproducible setting.
#' @param patience Stop after this many consecutive non-improving epochs plus
#'   one. A negative value disables early stopping.
#' @param deterministic Request deterministic PyTorch algorithms when
#'   available.
#' @param python Optional path to a Python executable, used before reticulate
#'   initializes Python.
#' @param condaenv Optional Conda environment name or path, used before
#'   reticulate initializes Python.
#' @param input_sep Field separator used when `expression` is a file path.
#' @param overwrite Whether existing GeneCluster outputs may be replaced.
#' @param verbose Whether the Python backend prints epoch-level progress.
#' @param return_features Whether to read the final feature matrix into R.
#'
#' @return An object of class `deepcscn_genecluster_result` containing paths,
#'   final clusters, training history, run metadata, and optionally features.
#' @export
run_genecluster <- function(expression,
                            output_dir,
                            n_clusters = 15L,
                            epochs = 100L,
                            batch_size = 32L,
                            learning_rate = 0.01,
                            momentum = 0.9,
                            weight_decay = 1e-5,
                            focal_gamma = 2,
                            seed = 0L,
                            device = "auto",
                            num_workers = 0L,
                            patience = 10L,
                            deterministic = TRUE,
                            python = NULL,
                            condaenv = NULL,
                            input_sep = "\t",
                            overwrite = FALSE,
                            verbose = TRUE,
                            return_features = TRUE) {
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create `output_dir`: ", output_dir, call. = FALSE)
  }

  if (is.character(expression) && length(expression) == 1L) {
    input_source <- normalizePath(path.expand(expression), mustWork = TRUE)
    matrix_data <- utils::read.table(
      input_source,
      header = TRUE,
      row.names = 1L,
      sep = input_sep,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
  } else if (is.matrix(expression) || is.data.frame(expression)) {
    matrix_data <- expression
  } else {
    stop("`expression` must be a numeric matrix/data frame or a file path.", call. = FALSE)
  }

  matrix_data <- as.matrix(matrix_data)
  storage.mode(matrix_data) <- "double"
  if (!is.numeric(matrix_data) || any(!is.finite(matrix_data))) {
    stop("`expression` must contain only finite numeric values.", call. = FALSE)
  }
  if (nrow(matrix_data) < 3L || ncol(matrix_data) < 2L) {
    stop("`expression` must contain at least three genes and two cells.", call. = FALSE)
  }
  if (is.null(rownames(matrix_data)) || anyNA(rownames(matrix_data)) ||
      any(!nzchar(rownames(matrix_data))) || anyDuplicated(rownames(matrix_data))) {
    stop("`expression` must have unique, non-empty gene row names.", call. = FALSE)
  }
  if (is.null(colnames(matrix_data))) {
    colnames(matrix_data) <- paste0("cell_", seq_len(ncol(matrix_data)))
  }
  if (anyNA(colnames(matrix_data)) || any(!nzchar(colnames(matrix_data))) ||
      anyDuplicated(colnames(matrix_data))) {
    stop("Cell column names must be unique and non-empty.", call. = FALSE)
  }
  n_clusters <- as.integer(n_clusters)
  if (length(n_clusters) != 1L || is.na(n_clusters) || n_clusters < 2L || n_clusters >= nrow(matrix_data)) {
    stop("`n_clusters` must be between 2 and nrow(expression) - 1.", call. = FALSE)
  }

  generated_files <- c(
    "input_expression.tsv", "features.tsv", "features.npy", "clusters.tsv",
    "training_history.tsv", "model.pt", "run_metadata.json"
  )
  existing <- file.path(output_dir, generated_files)
  existing <- existing[file.exists(existing)]
  if (length(existing) && !isTRUE(overwrite)) {
    stop(
      "GeneCluster output files already exist. Set `overwrite = TRUE` to replace them: ",
      paste(basename(existing), collapse = ", "),
      call. = FALSE
    )
  }

  input_path <- file.path(output_dir, "input_expression.tsv")
  utils::write.table(
    matrix_data,
    file = input_path,
    sep = "\t",
    quote = FALSE,
    col.names = NA,
    row.names = TRUE
  )

  set.seed(as.integer(seed))
  configure_deepcscn_python(python = python, condaenv = condaenv)
  backend <- .deepcscn_backend()
  result <- backend$train_genecluster(
    input_path = input_path,
    output_dir = output_dir,
    n_clusters = n_clusters,
    epochs = as.integer(epochs),
    batch_size = as.integer(batch_size),
    learning_rate = as.numeric(learning_rate),
    momentum = as.numeric(momentum),
    weight_decay = as.numeric(weight_decay),
    focal_gamma = as.numeric(focal_gamma),
    seed = as.integer(seed),
    device = as.character(device),
    num_workers = as.integer(num_workers),
    patience = as.integer(patience),
    deterministic = isTRUE(deterministic),
    overwrite = isTRUE(overwrite),
    verbose = isTRUE(verbose)
  )

  features <- NULL
  if (isTRUE(return_features)) {
    features <- as.matrix(utils::read.delim(
      result$features_path,
      row.names = 1L,
      check.names = FALSE
    ))
  }
  clusters <- utils::read.delim(result$clusters_path, check.names = FALSE)
  history <- utils::read.delim(result$history_path, check.names = FALSE)

  output <- list(
    features = features,
    clusters = clusters,
    history = history,
    files = list(
      input = input_path,
      features = result$features_path,
      features_npy = result$features_npy_path,
      clusters = result$clusters_path,
      history = result$history_path,
      model = result$model_path,
      metadata = result$metadata_path
    ),
    summary = result,
    call = match.call()
  )
  class(output) <- "deepcscn_genecluster_result"
  output
}

#' @export
print.deepcscn_genecluster_result <- function(x, ...) {
  cat("DeepCSCN GeneCluster result\n")
  cat("  genes:", x$summary$n_genes, "\n")
  cat("  cells:", x$summary$n_cells, "\n")
  cat("  epochs run:", x$summary$epochs_run, "\n")
  cat("  device:", x$summary$device, "\n")
  cat("  features:", x$files$features, "\n")
  invisible(x)
}
