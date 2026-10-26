#' Cell-Type-Specific Gene Co-expression Network
#'
#' This package provides functions to construct and analyze cell-type-specific networks.


#' Identify Marker Genes Using Seurat
#'
#' This function identifies marker genes for each group in the metadata
#' using Seurat's FindAllMarkers function.
#'
#' @param expr A numeric matrix or data frame of gene expression (genes x cells).
#' @param metadata_path Path to a CSV file containing metadata.
#'   The row names should match the column names of `expr`.
#'   Must contain a column named `time_info` for grouping.
#' @param output_dir Path to save the marker genes CSV file.
#'
#' @return A data frame containing marker genes identified by Seurat.
#'
#' @examples
#' \dontrun{
#' markers <- marker_identify(expr, "meta.csv", "markers.csv")
#' }
#'
#' @export
marker_identify <- function(expr, metadata_path, output_dir) {
  # seurat version: V4
  # Read metadata
  meta_data <- utils::read.csv(metadata_path, row.names = 1, header = TRUE)

  # Keep columns in expr consistent with metadata
  new_expr <- expr[, rownames(meta_data)]

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = new_expr, meta.data = meta_data)

  # Set identity classes
  Seurat::Idents(seurat_obj) <- meta_data$time_info

  # Find marker genes
  markers <- Seurat::FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,          # Only find upregulated genes
    min.pct = 0.25,           # Expressed in at least 25% of cells
    logfc.threshold = 0.5     # Minimum log fold change threshold
  )

  # Save results
  utils::write.csv(markers, output_dir, row.names = FALSE)

  return(markers)
}


#' Get the Top Ten Percent of Features for Each Cell Type
#'
#' @name CelltypeFeat
#' @param DEG A data frame of marker genes containing a `gene` column.
#' @param feat A gene-by-feature matrix from the GeneCluster model.
#' @param expr A gene-by-cell expression matrix.
#' @param meta A cell metadata data frame whose row names match expression columns.
#' @param CellType_col Name of the metadata column containing cell types or states.
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
    FeatExprCor <- stats::cor(new_feat, new_count, use = "p")

    abs_FeatExprCor <- abs(FeatExprCor)
    for (i in unique(meta[[CellType_col]])) {
        sample <- meta[meta[[CellType_col]] == i, , drop = FALSE]
        xcell_FeatExprCor <- abs_FeatExprCor[, colnames(abs_FeatExprCor) %in% rownames(sample)]
        top_feat <- data.frame(cor = sort(rowMeans(xcell_FeatExprCor), decreasing = TRUE)[1:51])
        xcell_feat <- feat[, rownames(top_feat)]
        utils::write.csv(xcell_feat, paste("data/res/feat/", i, "_feat.csv", sep = ''))
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
  # List input files
  filelist <- list.files(input_dir, full.names = TRUE)
  results <- list()

  for (file in filelist) {
    # Read features
    celltype_feat <- utils::read.csv(file, row.names = 1)

    # Compute cosine similarity
    celltype_cor <- lsa::cosine(as.matrix(t(celltype_feat)))

    # Generate output file name
    name <- sub("_feat.csv", "", basename(file))
    outfile <- file.path(output_dir, paste0(name, "_cormat.csv"))

    # Save correlation matrix
    utils::write.csv(celltype_cor, outfile)

    # Store result in list
    results[[name]] <- celltype_cor
  }

  return(results)
}


#' Perform First Clustering
#' @name FirstCluster
#' @param celltype_feat Selected cell type features.
#'
#' @return A data frame with gene clusters. Each gene is assigned a module number.
#' @export
FirstCluster <- function(celltype_feat){
    cor_matrix <- lsa::cosine(as.matrix(t(celltype_feat)))
    module_num <- 2
    hclust_res <- stats::cutree(stats::hclust(stats::as.dist(1 - cor_matrix)),module_num)
    hclust_res <- as.data.frame(hclust_res)
    colnames(hclust_res) <- "cluster"
    sort_df <- hclust_res[order(hclust_res$cluster), , drop = FALSE]
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
    relative_expr <- expr / rowSums(expr)

    # merge expression + metadata
    merged_data <- merge(t(relative_expr), meta, by = "row.names")
    rownames(merged_data) <- merged_data$Row.names
    merged_data$Row.names <- NULL

    module1 <- celltype_module[celltype_module$cluster == 1, , drop = FALSE]
    module1_genes <- rownames(module1)
    module2 <- celltype_module[celltype_module$cluster == 2, , drop = FALSE]
    module2_genes <- rownames(module2)

    ## --- Module 1 ---
    module1_expr <- merged_data[merged_data$time_info == target_time, module1_genes]
    all_expr1 <- merged_data[, module1_genes]
    mu_target1 <- mean(as.numeric(as.matrix(module1_expr)))
    mu_all1 <- mean(as.numeric(as.matrix(all_expr1)))
    sigma_all2 <- stats::sd(as.numeric(as.matrix(all_expr1)))
    score1 = abs(mu_target1 - mu_all1) / sigma_all2
    ## --- Module 2 ---
    module2_expr <- merged_data[merged_data$time_info == target_time, module2_genes]
    all_expr2 <- merged_data[, module2_genes]
    mu_target2 <- mean(as.numeric(as.matrix(module2_expr)))
    mu_all2 <- mean(as.numeric(as.matrix(all_expr2)))
    sigma_all2 <- stats::sd(as.numeric(as.matrix(all_expr2)))
    score2 = abs(mu_target2 - mu_all2) / sigma_all2

    ## --- Compare scores ---
    if(score1 > score2){
        message("M1 has the higher specificity score (", signif(score1 / score2, 4), "-fold over M2).")
    }else if(score2 > score1){
        message("M2 has the higher specificity score (", signif(score2 / score1, 4), "-fold over M1).")
    }else{
        message("M1 and M2 have equal specificity scores.")
    }

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
    dist_cor <- 1 - cor
    hclust_dist <- stats::hclust(stats::as.dist(dist_cor), method = "average")
    memb <- dynamicTreeCut::cutreeDynamic(dendro = hclust_dist,
                                            distM = dist_cor,
                                            deepSplit = 2,
                                            pamRespectsDendro = FALSE,
                                            minClusterSize = minClusterSize)
    names(memb) <- colnames(cor)
    memb_df <- as.data.frame(memb)
    colnames(memb_df) <- "cluster"
    sort_df <- memb_df[order(memb_df$cluster), , drop = FALSE]
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
    lab2 <- ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = index, labels = labs, labels_gp = grid::gpar(fontsize = 9), lines_gp = grid::gpar()))
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#0084ff", "white", "#ff0404"))

    p <- ComplexHeatmap::Heatmap(cor_matrix, col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE,
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
#' @return A data frame of GO enrichment results.
#' @export
Run_GO <- function(Cluster_res) {
    df1 <- data.frame()
    for (g in unique(Cluster_res$cluster)) {
        module_gene <- Cluster_res[Cluster_res$cluster == g, , drop = FALSE]
        genelist <- rownames(module_gene)
        GOres <- clusterProfiler::enrichGO(genelist, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = 'SYMBOL', ont = "ALL",
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
    my_palette <- grDevices::colorRampPalette(colors = c("#7b67ed", "white"))(length(bk <- seq(0, 1e-3, by = 1e-4)))
    p <- pheatmap::pheatmap(t(plot_df), cellwidth = 20, cellheight = 20, cluster_rows = FALSE,
                   cluster_cols = FALSE, display_numbers = FALSE, number_color = "white",
                   fontsize_row = 15, fontsize_col = 15, angle_col = 45,
                   fontsize_number = 15, color = my_palette, breaks = bk,
                   legend_breaks = seq(0, 1e-3, 1e-4))
    return(p)
}



