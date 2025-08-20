#' Cell-Type-Specific Gene Co-expression Network
#'
#' This package provides functions to construct and analyze cell-type-specific networks.
#'
#' @docType package
#' @name DeepCSCN


#' Identify Marker Genes Using Seurat
#'
#' This function identifies marker genes for each group in the metadata
#' using Seurat's `FindAllMarkers()` function.
#'
#' @param expr A numeric matrix or data frame of gene expression (genes x cells).
#' @param metadata_path Path to a CSV file containing metadata.
#'   The row names should match the column names of `expr`.
#'   Must contain a column named `time_info` for grouping.
#' @param output_dir Path to save the marker genes CSV file.
#'
#' @return A data frame containing marker genes identified by Seurat.
#' @examples
#' \dontrun{
#' markers <- marker_identify(expr, "meta.csv", "markers.csv")
#' }
#' @export
marker_identify <- function(expr, metadata_path, output_dir) {
  # seurat version: V4
  library(Seurat)

  # Read metadata
  meta_data <- read.csv(metadata_path, row.names = 1, header = TRUE)

  # Keep columns in expr consistent with metadata
  new_expr <- expr[, rownames(meta_data)]

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = new_expr, meta.data = meta_data)

  # Set identity classes
  Idents(seurat_obj) <- meta_data$time_info

  # Find marker genes
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,          # Only find upregulated genes
    min.pct = 0.25,           # Expressed in at least 25% of cells
    logfc.threshold = 0.5     # Minimum log fold change threshold
  )

  # Save results
  write.csv(markers, output_dir, row.names = FALSE)

  return(markers)
}


#' Get the Top 10% Features for Each Cell Type
#' 
#' @name CelltypeFeat
#' @param DEG_path Marker genes for different cell types.
#' @param feat A feat matrix from GeneCluster model.
#' @param expr expression matrix.
#'
#' @export
CelltypeFeat <- function(DEG, feat, expr, meta, CellType_col) {
    if (!dir.exists("data/res/feat/")) {
    dir.create("data/res/feat/", recursive = TRUE)
    }

    new_feat <- feat[DEG$gene, ]
    # feat_filtered <- new_feat[rowSums(new_feat > 0) > 0, colSums(new_feat > 0) > 0] 
    new_count <- expr[DEG$gene, rownames(meta)] 
    # count_filtered <- new_count[rowSums(new_count > 0) > 0, colSums(new_count > 0) > 0]
    FeatExprCor <- cor(new_feat, new_count, use = "p")    

    abs_FeatExprCor <- abs(FeatExprCor)
    for (i in unique(meta[[CellType_col]])) {
        sample <- meta[meta[[CellType_col]] == i, , drop = FALSE]
        xcell_FeatExprCor <- abs_FeatExprCor[, colnames(abs_FeatExprCor) %in% rownames(sample)]
        top_feat <- data.frame(cor = sort(rowMeans(xcell_FeatExprCor), decreasing = TRUE)[1:51])
        xcell_feat <- feat[, rownames(top_feat)]
        write.csv(xcell_feat, paste("data/res/feat/", i, "_feat.csv", sep = ''))
    }
}


#' Compute Cosine Similarity Correlation Matrices
#'
#' This function reads feature matrices (CSV files) from an input directory,
#' computes pairwise cosine similarity for each file, 
#' and saves the correlation matrices as CSV files in the output directory.
#'
#' @param input_dir Path to the input directory containing feature CSV files.
#' @param output_dir Path to the output directory where correlation matrices will be saved.
#'
#' @return A list of correlation matrices (one for each input file).
#' @examples
#' \dontrun{
#' compute_cosine_cormat(
#'   input_dir = "hESC_GCN/time_feat/v2/",
#'   output_dir = "hESC_GCN/time_cormat/v2/"
#' )
#' }
#' @export
compute_cosine_cormat <- function(input_dir, output_dir) {
  library(lsa)

  # List input files
  filelist <- list.files(input_dir, full.names = TRUE)
  results <- list()

  for (file in filelist) {
    # Read features
    celltype_feat <- read.csv(file, row.names = 1)

    # Compute cosine similarity
    celltype_cor <- cosine(as.matrix(t(celltype_feat)))

    # Generate output file name
    name <- sub("_feat.csv", "", basename(file))
    outfile <- file.path(output_dir, paste0(name, "_cormat.csv"))

    # Save correlation matrix
    write.csv(celltype_cor, outfile)

    # Store result in list
    results[[name]] <- celltype_cor
  }

  return(results)
}


#' Perform First Clustering
#' @name FirstCluster
#' @param cor_matrix A correlation matrix of genes.
#' @param celltype_feat Selected cell type features.
#'
#' @return A data frame with gene clusters. Each gene is assigned a module number.
#' @export
FirstCluster <- function(celltype_feat){
    cor_matrix <- cosine(as.matrix(t(celltype_feat)))
    module_num <- 2
    hclust_res <- cutree(hclust(as.dist(1-cor_matrix)),module_num)
    hclust_res <- as.data.frame(hclust_res)
    colnames(hclust_res) <- "cluster"
    sort_df <- arrange(hclust_res,cluster)
    return(list(celltype_cor = cor_matrix, celltype_cluster = sort_df))
}


#' Compute Cell Type Specificity Scores for Gene Modules
#'
#' This function calculates z-like specificity scores for two gene modules 
#' (e.g., cluster 1 and cluster 2) at a given time point, based on relative gene expression.
#'
#' @param expr A gene-by-cell expression matrix (numeric matrix or data frame).
#' @param meta A metadata data frame, rownames must match cell names in \code{expr}, 
#'   and must contain a column \code{time_info}.
#' @param celltype_module A data frame of gene modules, with rownames as gene symbols 
#'   and a column \code{cluster} indicating module assignment (e.g., 1 or 2).
#' @param target_time The time point (string) in \code{meta$time_info} for which 
#'   the scores should be computed.
#'
#' @return A named list with two elements:
#' \item{M1_score}{Specificity score for module 1.}
#' \item{M2_score}{Specificity score for module 2.}
#'
#' @examples
#' \dontrun{
#' scores <- CellTypeSpecificScore(
#'   expr = expr_matrix,
#'   meta = metadata,
#'   celltype_module = module_df,
#'   target_time = "24h"
#' )
#' print(scores$M1_score)
#' print(scores$M2_score)
#' }
#' @export
CellTypeSpecificScore <- function(expr, meta, celltype_module, target_time){
    # relative expression per gene
    relative_expr = expr/rowSums(expr)

    # merge expression + metadata
    merged_data <- merge(t(relative_expr), meta, by = "row.names")
    rownames(merged_data) <- merged_data$Row.names
    merged_data$Row.names <- NULL    

    module1 = subset(celltype_module,cluster==1)
    module1_genes = rownames(module1)
    module2 = subset(celltype_module,cluster==2)
    module2_genes = rownames(module2)

    ## --- Module 1 ---
    module1_expr <- merged_data[merged_data$time_info == target_time, module1_genes]
    all_expr1 <- merged_data[, module1_genes]
    mu_target1 <- mean(as.numeric(as.matrix(module1_expr)))
    mu_all1 <- mean(as.numeric(as.matrix(all_expr1)))
    sigma_all2 <- sd(as.numeric(as.matrix(all_expr1)))
    score1 = abs(mu_target1 - mu_all1) / sigma_all2
    ## --- Module 2 ---
    module2_expr <- merged_data[merged_data$time_info == target_time, module2_genes]
    all_expr2 <- merged_data[, module2_genes]
    mu_target2 <- mean(as.numeric(as.matrix(module2_expr)))
    mu_all2 <- mean(as.numeric(as.matrix(all_expr2)))
    sigma_all2 <- sd(as.numeric(as.matrix(all_expr2)))
    score2 = abs(mu_target2 - mu_all2) / sigma_all2

    return(list(M1_score=score1, M2_score=score2)) 
}


#' Identify Specific Clusters from Correlation Matrix
#'
#' This function applies hierarchical clustering and dynamic tree cutting 
#' (via \pkg{WGCNA} and \pkg{dynamicTreeCut}) to identify clusters from a correlation matrix.  
#'
#' @param cor A correlation matrix (square numeric matrix or data frame, 
#'   rownames and colnames should match).
#' @param minClusterSize Minimum cluster size (default = 10).
#'
#' @return A data frame with one column:
#' \item{cluster}{Cluster assignment for each variable (columns of \code{cor}).}
#' Row names correspond to variable names.
#'
#' @examples
#' \dontrun{
#' cor_mat <- cor(matrix(rnorm(1000), nrow = 50))
#' clusters <- specific_cluster(cor_mat, minClusterSize = 5)
#' head(clusters)
#' }
#' @export
specific_cluster <- function(cor,minClusterSize=10) {
    library(WGCNA)
    dist_cor <- 1 - cor
    hclust_dist <- hclust(as.dist(dist_cor), method = "average")
    memb <- dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                            distM = dist_cor, 
                                            deepSplit = 2, 
                                            pamRespectsDendro = FALSE, 
                                            minClusterSize = minClusterSize)
    names(memb) <- colnames(cor)
    memb_df <- as.data.frame(memb)
    colnames(memb_df) <- "cluster"
    sort_df <- arrange(memb_df, cluster) 
    return(sort_df)
}


#' Heatmap Plotting
#'
#' This function generates a heatmap for the correlation matrix.
#'
#' @name heatmap_plot
#' @param cor_matrix A correlation matrix.
#' @param marker A data frame containing marker genes.
#' @return A heatmap object.
#' @export
heatmap_plot <- function(cor_matrix, marker) {
    genelist <- marker$gene
    index <- which(rownames(cor_matrix) %in% genelist)
    labs <- rownames(cor_matrix)[index]
    lab2 <- rowAnnotation(foo = anno_mark(at = index, labels = labs, labels_gp = gpar(fontsize = 9), lines_gp = gpar()))
    col_fun <- colorRamp2(c(-1, 0, 1), c("#0084ff", "white", "#ff0404"))
    
    p <- Heatmap(cor_matrix, col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE, 
                  show_row_names = FALSE, show_column_names = FALSE, 
                  name = "Cor", width = 0.5, height = 0.5, right_annotation = lab2)
    return(p)
}



#' Run GO Enrichment Analysis
#'
#' This function performs GO enrichment analysis on clustered genes.
#'
#' @name Run_GO
#' @param Cluster_res A data frame containing cluster results.
#' @param universe A character vector of universe genes.
#' @return A data frame of GO enrichment results.
#' @export
Run_GO <- function(Cluster_res) {
    df1 <- data.frame()
    for (g in unique(Cluster_res$cluster)) {
        module_gene <- subset(Cluster_res, cluster == g)
        genelist <- rownames(module_gene)
        GOres <- enrichGO(genelist, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "ALL", 
                          pvalueCutoff = 0.05, qvalueCutoff = 0.05)
        if (dim(GOres)[1] > 0) {
            go_resdf <- data.frame(GOres)
            # top3go <- go_resdf[order(go_resdf$p.adjust), ][1:3, ]
            go_resdf$module <- rep(g, nrow(go_resdf))
            df1 <- rbind(df1, go_resdf)            
        }
    }
    return(df1)
}

# GO Heatmap
GO_heatmap <- function(plot_df) {
    my_palette <- colorRampPalette(colors = c("#7b67ed", "white"))(length(bk <- seq(0, 1e-3, by = 1e-4)))
    p <- pheatmap(t(plot_df), cellwidth = 20, cellheight = 20, cluster_rows = FALSE, 
                   cluster_cols = FALSE, display_numbers = FALSE, number_color = "white", 
                   fontsize_row = 15, fontsize_col = 15, angle_col = 45, 
                   fontsize_number = 15, color = my_palette, breaks = bk, 
                   legend_breaks = seq(0, 1e-3, 1e-4))
    return(p) 
}



