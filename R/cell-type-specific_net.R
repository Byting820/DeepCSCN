#' Cell Type-Specific Gene Co-expression Network
#'
#' This module contains functions for constructing cell type-specific gene co-expression networks.
#'
#' @title Cell Type-Specific Gene Co-expression Network
#' @author Byting
#' @import dplyr
#' @import lsa
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import pheatmap
#' @import ComplexHeatmap
#' @export


#' Get the Top 100 Features for Each Cell Type
#' 
#' @name CelltypeFeat
#' @param DEG_path Marker genes for different cell types.
#' @param feat A feat matrix from GeneCluster model.
#' @param count Raw count matrix.
#'
#' @export
CelltypeFeat <- function(DEG_path, feat, count,sort_meta) {
    if (!dir.exists("res/")) {
    dir.create("res/", recursive = TRUE)
    }

    DEG <- read.csv(DEG_path, row.names = 1)
    new_feat <- feat[rownames(DEG), ]
    new_count <- count[rownames(DEG), rownames(sort_meta)]
    FeatExprCor <- cor(new_feat, new_count, use = "p")    

    abs_FeatExprCor <- abs(FeatExprCor)
    for (i in unique(sort_meta$CellType)) {
        sample <- sort_meta[sort_meta$CellType == i, ]
        xcell_FeatExprCor <- abs_FeatExprCor[, colnames(abs_FeatExprCor) %in% rownames(sample)]
        top100_feat <- data.frame(cor = sort(rowMeans(xcell_FeatExprCor), decreasing = TRUE)[1:100])
        xcell_feat <- feat[, colnames(feat) %in% rownames(top100_feat)]
        write.csv(xcell_feat, paste("res/", i, "_feat.csv", sep = ''))
    }
}

#' Perform First Clustering
#' @name selectCellModules
#' @param cor_matrix A correlation matrix of genes.
#' @param celltype_feat Selected cell type features.
#'
#' @return A data frame with gene clusters. Each gene is assigned a module number.
#' @export
selectCellModules <- function(celltype_feat){
    cor_matrix <- cosine(as.matrix(t(celltype_feat)))
    module_num <- 2
    hclust_res <- cutree(hclust(as.dist(1-cor_matrix)),module_num)
    hclust_res <- as.data.frame(hclust_res)
    colnames(hclust_res) <- "cluster"
    return(list(celltype_cor = cor_matrix, celltype_cluster = hclust_res))
}


#' Calculate Module-Specific Scores for a Given Cell Type
#'
#' This function calculates module-specific scores (expression sums) for two gene modules
#' within a specific cell type, using normalized counts.
#' @name calcModuleScore
#' @param count A matrix of gene expression counts with genes as rows and samples as columns.
#' @param celltype_module The clustering result is obtained by selectCellModules function. 
#' for the specified cell type.
#' @param celltye_name A character string specifying the cell type to calculate scores for.
#' @param meta A metadata data frame with cell type information for each sample.
#'
#' @return A list with `M1_expr` and `M2_expr`, which represent the average expression sums for modules 1 and 2, respectively.
#' @export
#'
calcModuleScore <- function(count,celltype_module,celltye_name,meta){
    relative_count = count
    sum_count = data.frame(sum_count = rowSums(relative_count))
    for (i in 1:ncol(relative_count)) {
        relative_count[,i] <- relative_count[,i]/sum_count
    }

    # celltype_module = read.csv(celltype_first_cluster_path,row.names=1)
    celltype_meta = meta[meta$CellType == celltye_name,]
    module1 = subset(celltype_module,cluster==1)
    module2 = subset(celltype_module,cluster==2)

    M1_count = relative_count[row.names(module1),row.names(celltype_meta)]
    M2_count = relative_count[row.names(module2),row.names(celltype_meta)]
    M1_gene_exprSum = (sum(rowSums(M1_count)))/(dim(M1_count)[1])
    M2_gene_exprSum = (sum(rowSums(M2_count)))/(dim(M2_count)[1])
    return(list(M1_expr=M1_gene_exprSum,M2_expr=M2_gene_exprSum)) 
}


#' Cluster Using Dynamic Tree Cut
#'
#' This function performs clustering on a correlation matrix using dynamic tree cut.
#'
#' @name second_cluster
#' @param cor A correlation matrix.
#' @return A data frame of clusters.
#' @export
second_cluster <- function(cor,minClusterSize=10) {
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


#' Identify Cell Type-Specific Modules Based on Expression Scores
#'
#' This function identifies the most relevant module for a specified cell type based on module-specific expression scores. 
#' It selects either module 1 or module 2 based on their average expression scores, and further clusters 
#' the associated genes for the chosen module to refine the cell type-specific network.
#'
#' @name CelltypeClusterRes
#' @param M1_expr Numeric; the average expression score for module 1 within the specified cell type.
#' @param M2_expr Numeric; the average expression score for module 2 within the specified cell type.
#' @param Cluster_res Data frame; initial clustering results for genes associated with the specified cell type. 
#' This data frame includes gene identifiers and their assigned clusters.
#' @param celltype_cor Matrix; a correlation matrix of gene expression specific to the selected cell type. 
#' Each entry represents the correlation between gene pairs.
#'
#' @return A refined clustering result (`Cluster_res`) for the selected module, which represents the cell type-specific gene network.
#' @export
#'
CelltypeClusterRes <- function(M1_expr, M2_expr,Cluster_res,celltype_cor){
    # Determine which module (1 or 2) is more relevant based on average expression scores
    if (M1_expr > M2_expr) {
        celltype_module <- 1
    } else {
        celltype_module <- 2
    }
    # Extract genes belonging to the selected module
    celltype_relative_gene <- rownames(subset(Cluster_res, cluster == celltype_module))
    # Subset the correlation matrix to only include genes in the selected module
    celltype_relative_cormat <- celltype_cor[celltype_relative_gene, celltype_relative_gene]
    # Perform secondary clustering on the refined set of genes
    Cluster_res <- second_cluster(celltype_relative_cormat)
    
    return(Cluster_res)
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
        if (dim(GOres)[1] >= 3) {
            go_resdf <- data.frame(GOres)
            top3go <- go_resdf[order(go_resdf$p.adjust), ][1:3, ]
            top3go$module <- rep(g, nrow(top3go))
            df1 <- rbind(df1, top3go)            
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


#' Cell-type-specific network
#'
#' This function obtains cell type specific clustering results. 
#' @param celltype_feat_path A character string specifying the file path to the cell type-specific feature data in CSV format.
#' Each row represents a gene, and the columns represent features or conditions relevant to the cell type.
#' @param celltye_name A character string specifying the name of the cell type being analyzed.
#' @param count A matrix of gene expression counts with genes as rows and samples as columns.
#' This matrix provides the expression data for calculating module-specific scores.
#' @param meta A data frame containing metadata information for each sample, such as cell type labels.
#'
#' @return A data frame representing the refined clustering result for the selected module within the specified cell type.
#' Each row corresponds to a gene and includes clustering information specific to the cell type.
#' @export
#'
#' @examples
#' # Example usage:
#' # CelltypeNet("data/celltype_features.csv", "B cell", count_matrix, meta_data)
#'
CelltypeNet <- function(celltype_feat_path, celltye_name,count,meta){
    celltype_feat <- read.csv(celltype_feat_path, row.names = 1)
    Cluster_res = selectCellModules(celltype_feat) 

    score = calcModuleScore(count,Cluster_res$celltype_cluster,celltye_name,meta)
    M1_expr = score$M1_expr
    M2_expr = score$M2_expr
    # print(sprintf("module1 score is %s; module2 score is %s", M1_expr, M2_expr))

    celltype_net = CelltypeClusterRes(M1_expr, M2_expr,Cluster_res$celltype_cluster,Cluster_res$celltype_cor)
    return(celltype_net)
}
