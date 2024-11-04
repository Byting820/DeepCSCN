#' Data Preprocessing
#'
#' This script is used to read and process expression data.
#' 
#' @title Data Preprocess
#' @author Byting
#' @param input_file A string representing the path to the input file.
#' @param num_HVG An integer specifying the number of highly variable genes to find.
#' @return A data frame containing the normalized expression data.
#' @export

process_expression_data <- function(input_file, num_HVG) {
    library(Seurat)
    library(dplyr)
    library(Matrix)

    # Quality control filter function
    qcFilter <- function(adata) {
        adata <- adata[rowSums(adata) > 0, ]  # Remove non-expressed genes
        adata <- adata[, colSums(adata) > 0]  # Remove non-expressed cells
        return(as.data.frame(adata))
    }

    # Find highly variable genes function
    findHVG <- function(expr, num_HVG) {
        adata <- CreateSeuratObject(counts = expr)
        adata <- NormalizeData(adata)
        adata <- FindVariableFeatures(adata, nfeatures = num_HVG)
        hvg_genes <- VariableFeatures(adata)
        
        return(expr[, hvg_genes])
    }

    # Min-max normalization function
    min_max_normalization <- function(data) {
        normalized_data <- apply(data, 2, function(x) (x - min(x)) / (max(x) - min(x)))
        return(as.data.frame(t(normalized_data)))
    }

    # Main processing logic
    if (grepl("\\.h5ad$", input_file)) {
        # Load h5ad data
        adata <- ReadH5AD(input_file)
        expr <- as.data.frame(adata@assays$RNA@counts)  # Assuming RNA counts
    } else if (grepl("\\.csv$", input_file)) {
        # Load CSV data
        expr <- read.csv(input_file, sep='\t', row.names = 1)
    } else {
        stop("Unsupported file format. Please provide a .h5ad or .csv file.")
    }
    
    # Quality control
    expr <- qcFilter(expr)
    
    # Find highly variable genes
    expr_hvg <- findHVG(expr, num_HVG)
    
    # Normalize data
    normalized_expr <- min_max_normalization(expr_hvg)
    
    return(normalized_expr)
}

# Example usage
# normalized_data <- process_expression_data("data/raw/pbmc1-Drop-1000hvg.csv", num_HVG = 1000)
# write.csv(normalized_data, "data/processed/normal_data.csv", sep = '\t')
