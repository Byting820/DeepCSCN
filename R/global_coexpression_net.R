#' Gene Co-expression Network Construction
#'
#' This module contains functions for constructing gene co-expression networks.
#'
#' @title Gene Co-expression Network
#' @author Byting
#' @import Seurat
#' @import DT
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @export

# Function to perform correlation clustering
#'
#' This function performs clustering on a correlation matrix using dynamic tree cut.
#'
#' @name global_net
#' @param feat_matrix A frature matrix from GeneCluster model.
#' @return Data frames of cormatrix and global cluster.
#' @export
global_net <- function(feat_matrix,minClusterSize=10) {
    cor <- cosine(as.matrix(t(feat_matrix)))
    dist_cor = 1 - cor
    hclust_dist = hclust(as.dist(dist_cor), method = "average") 
    memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                          distM = dist_cor, 
                                          deepSplit = 2,
                                          pamRespectsDendro = FALSE,
                                          minClusterSize = minClusterSize)
    names(memb) = colnames(cor)
    memb_df = as.data.frame(memb)
    colnames(memb_df) = "cluster"
    sorted_clusters <- arrange(memb_df, cluster)
    return(list(cor_mat = cor,global_cluster=sorted_clusters))
}

# module visualization
module_visual <- function(cor_data,sorted_clusters){
    plot_data <- cor_data[rownames(sorted_clusters), rownames(sorted_clusters)]
    diag(plot_data) <- NA
    # Define heatmap annotations
    col_anno <- HeatmapAnnotation(module = as.factor(sorted_clusters$cluster),
                                col = list(module = c("1" = "#ffab66", "2" = "#ed6a6c", "3" = "#cfdd7b",
                                                        "4" = "#e1abd1", "5" = "#fac0bc", "6" = "#aadfc4",
                                                        "7" = "#edd99f", "8" = "#bbbbdfd5", "9" = "#f0ffbb",
                                                        "10" = "#c0c4c7", "11" = "#ffb01e", "12" = "#ed7d1c",
                                                        "13" = "#dc717c", "14" = "#b39869", "15" = "#4a96d0",
                                                        "16" = "#a2b3bf", "17" = "#c5b0d5", "18" = "#c49c94",
                                                        "19" = "#e377c2", "20" = "#7f7f7f", "21" = "#c7c7c7",
                                                        "22" = "#bcbd22", "23" = "#17becf", "24" = "#9edae5",
                                                        "25" = "#aec7e8")))

    ha2 <- rowAnnotation(module = as.factor(sorted_clusters$cluster),
                        col = list(module = c("1" = "#ffab66", "2" = "#ed6a6c", "3" = "#cfdd7b",
                                            "4" = "#e1abd1", "5" = "#fac0bc", "6" = "#aadfc4",
                                            "7" = "#edd99f", "8" = "#bbbbdfd5", "9" = "#f0ffbb",
                                            "10" = "#c0c4c7", "11" = "#ffb01e", "12" = "#ed7d1c",
                                            "13" = "#dc717c", "14" = "#b39869", "15" = "#4a96d0",
                                            "16" = "#a2b3bf", "17" = "#c5b0d5", "18" = "#c49c94",
                                            "19" = "#e377c2", "20" = "#7f7f7f", "21" = "#c7c7c7",
                                            "22" = "#bcbd22", "23" = "#17becf", "24" = "#9edae5",
                                            "25" = "#aec7e8")))

    # Create heatmap
    pdf("/Module_Heatmap.pdf", width = 6, height = 5)
    Heatmap(as.matrix(plot_data),
            col = gplots::colorpanel(250, 'lemonchiffon', "orange"),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_names = FALSE,
            show_column_names = FALSE,
            heatmap_legend_param = list(title = "cor", 
                                        title_position = "topleft", 
                                        at = c(0, 0.5, 1),
                                        legend_height = unit(1.5, "cm"), 
                                        legend_direction = "vertical"),
            left_annotation = ha2,
            top_annotation = col_anno)
    dev.off()    
}


# Identify Marker Genes Associated with Cell Types
marker_identify <- function(expr_path, metadata_path) {
    library(Seurat)
    library(SeuratData)

    expr_data <- read.table(expr_path, header = TRUE)
    meta_data <- read.table(metadata_path, sep = '\t', header = TRUE)
    meta <- meta_data[((meta_data$Experiment == 'pbmc1') & (meta_data$Method == 'Drop-seq')), ] 
    rownames(meta) <- meta$NAME_TYPE
    new_expr_data <- subset(expr_data, select = c(meta$NAME_TYPE)) 

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = new_expr_data, meta.data = meta) 
    Idents(seurat_obj) <- meta$CellType 
    markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(markers, 'markerGene.csv', row.names = FALSE)
}

#' Histogram Enrichment Plot
#'
#' This function creates a heatmap to visualize p-values and cross-table data.
#'
#' @name hist_enrichment_plot
#' @param pvalue A matrix of p-values.
#' @param cross_table A matrix of cross-table data.
#' @return A heatmap object.
#' @export
hist_enrichment_plot <- function(pvalue, cross_table) {
    ha1 = HeatmapAnnotation(
        dist1 = anno_barplot(
            colSums(cross_table),
            bar_width = 1,
            gp = gpar(col = "white", fill = "#FFE200"),
            border = FALSE,
            axis_param = list(at = c(0, 20, 40, 60), labels = c("0", "20", "40", "60")),
            height = unit(2.5, "cm")
        ), show_annotation_name = FALSE
    )

    ha2 = rowAnnotation(
        dist2 = anno_barplot(
            rowSums(cross_table),
            bar_width = 1,
            gp = gpar(col = "white", fill = "#FFE200"),
            border = FALSE,
            axis_param = list(at = c(0, 20, 40, 60, 80), labels = c("0", "20", "40", "60", "80")),
            width = unit(2.5, "cm")
        ), show_annotation_name = FALSE
    )

    col_fun = circlize::colorRamp2(c(0, 2, 10), c("white", "yellow", "red"))
    p <- Heatmap(pvalue,
                  name = "-log10P",
                  col = col_fun,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  top_annotation = ha1,
                  right_annotation = ha2,
                  width = ncol(pvalue) * unit(8, "mm"),
                  height = nrow(pvalue) * unit(8, "mm"),
                  row_names_side = "left",
                  column_labels = colnames(pvalue),
                  column_names_rot = 0,
                  border = "black",
                  heatmap_legend_param = list(
                      title = "-log10P", at = c(0, 2, 4, 6, 8, 10)
                  ),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      if (pvalue[i, j] > 2) {
                          grid.text(sprintf("%d", cross_table[i, j]), x, y, gp = gpar(fontsize = 12))
                      }
                      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
                  })
    return(p)
}


#' Function to associate modules with cell types
#'
#' This function associates global gene modules with cell types. 
#'
#' @name CelltypeModuleMap
#' @param clusterRes Global cluster res.
#' @param marker Celltype marker gene.
#' @return Celltype-module assocation matrix and P-value.
#' @export

CelltypeModuleMap <- function(clusterRes, marker){
    hclust_res <- read.csv(clusterRes)
    markergene <- read.csv(marker)
    colnames(hclust_res) <- c("gene", "module")
    
    # Merge data
    merged_genes <- merge(hclust_res, markergene, by = "gene", all = FALSE)
    cross_table <- table(merged_genes$module, merged_genes$cluster)
    pvalue <- cross_table

    # Calculate p-values
    for (i in rownames(cross_table)) {
        for (j in colnames(cross_table)) {
            p <- phyper((cross_table[i, j] - 1), 
                        sum(cross_table[, j]), 
                        (nrow(hclust_res) - sum(cross_table[, j])), 
                        sum(cross_table[i, ]), 
                        lower.tail = FALSE)
            q <- p.adjust(p, method = "fdr")
            pvalue[i, j] <- q
        }
    }
    pvalue <- -1 * log10(pvalue)
    return(list(crosstable = cross_table, pvalue = pvalue))  
}
