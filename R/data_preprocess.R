#' Plot Expression Density
#'
#' This function generates a density plot of gene expression levels,
#' highlighting the peak value with a dashed red line.
#'
#' @param plot_df A data frame containing at least a column `mean` representing
#' the average expression level of each gene.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#'   expr_data <- data.frame(mean = runif(100, 0, 10))
#'   p <- expr_density_plot(expr_data)
#'   print(p)
#' }
#' @export
expr_density_plot <- function(plot_df) {
  library(ggplot2)

  # Calculate density
  density_info <- density(plot_df$mean)
  peak_x <- density_info$x[which.max(density_info$y)]
  peak_y <- max(density_info$y)
  y_max <- peak_y * 1.1

  p <- ggplot(plot_df, aes(x = mean)) +
    geom_density(fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = peak_x, linetype = "dashed", color = "red") +
    annotate("text", x = peak_x, y = peak_y,
             label = round(peak_x, 2),
             vjust = -0.5, color = "red") +
    labs(x = "Gene Expression level", y = "Density") +
    ylim(0, y_max) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    )
  return(list(plot=p,
              peak_x=peak_x))
}

#' Filter Expression Matrix by Metadata and Cutoff
#'
#' This function filters an expression matrix based on metadata grouping
#' and a specified expression cutoff. Genes are retained if their mean
#' expression is above the cutoff in at least one group.
#'
#' @param expr A numeric matrix or data frame of gene expression values
#' with genes in rows and samples in columns.
#' @param meta_path Path to a CSV file containing metadata.
#' The first column should match the column names of \code{expr}.
#' Must contain a column named \code{time_info} for grouping.
#' @param cutoff Numeric threshold for mean expression filtering.
#'
#' @return A filtered expression matrix containing only genes
#' passing the cutoff in at least one group.
#'
#' @examples
#' \dontrun{
#'   filtered_expr <- filter_expr(expr, "meta.csv", 0.15)
#'   write.csv(filtered_expr, "filtered_expr.csv")
#' }
#' @export
filter_expr <- function(expr, meta_path, cutoff) {
  library(dplyr)

  meta <- read.csv(meta_path, sep = ",", row.names = 1)

  merged_data <- merge(t(expr), meta, by = "row.names")
  rownames(merged_data) <- merged_data$Row.names
  merged_data$Row.names <- NULL

  mean_expression <- merged_data %>%
    group_by(time_info) %>%
    summarise_all(mean) %>%
    t()

  colnames(mean_expression) <- mean_expression[1, ]
  mean_expression <- mean_expression[-1, ]

  filtered_matrix <- mean_expression[rowSums(mean_expression >= cutoff) > 0, ]
  filtered_expr <- expr[rownames(filtered_matrix), ]

  return(filtered_expr)
}


#' Find Highly Variable Genes using Python (Scanpy)
#'
#' This function calls a Python implementation of HVG selection via \code{scanpy},
#' using \code{reticulate} to run Python code within R.
#'
#' @param expr A data frame or matrix (genes as columns, cells as rows).
#' @param num_HVG Integer number of highly variable genes to return.
#' @param python_env Path to the Python environment containing scanpy and numpy.
#'
#' @return A data frame containing only the HVGs.
#' @examples
#' \dontrun{
#'   res <- findHVG_py(expr, num_HVG = 1000, python_env = "/path/to/conda_env")
#' }
#' @export
findHVG_py <- function(expr, num_HVG, python_env = NULL) {
  library(reticulate)

  if (!is.null(python_env)) {
    use_python(python_env, required = TRUE)
  }

  # load Python script
  hv_module <- import_from_path(
    module = "findHVG",
    path = system.file("python", package = "DeepCSCN")
  )

  # transfer expr to pandas.DataFrame（cell x gene）
  pd <- import("pandas")
  expr_df <- pd$DataFrame(expr)

  # call Python function
  result <- hv_module$findHVG(expr_df, as.integer(num_HVG))

  # transfer R data.frame
  return(as.data.frame(result))
}
