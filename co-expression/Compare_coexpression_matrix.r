# Description: Comparison of co-expression matrix of different methods.

library(GENIE3)
library(PIDC)
library(minet)
library(QUIC)
library(lsa)
library(scLink)
library(CSCORE)

comp_cor <- function(name,data){
    if (name == "GENIE3"){
        cov_mat <- GENIE3(t(data))  # requires gene x cell
    } else if (name == "minet"){
        cov_mat <- minet(data)
    } else if (name == "PIDC"){
        gene_sparsity <- colSums(data)
        non_zero_idx <- which(gene_sparsity != 0)
        cov_mat <- matrix(0.0, dim(data)[2], dim(data)[2])
        row.names(cov_mat) <- colnames(cov_mat) <- colnames(data)
        diag(cov_mat) <- 1.0
        tmp_data <- data[, non_zero_idx]
        tmp_cov_mat <- PIDC(t(tmp_data), verbose = FALSE) # requires gene x cell
        cov_mat[non_zero_idx, non_zero_idx] <- tmp_cov_mat
    } else if (name == "pearson"){
        cov_mat <- cor(data, method = "pearson")
    } else if (name == "spearman"){
        cov_mat <- cor(data, method = "spearman")
    } else if (name == "scLink"){
        cov_mat <- sclink_cor(data, ncores = 1, nthre = 10, dthre = 0.9)
    } else if (name == "GeneCluster"){
        cov_mat <- cosine(as.matrix(t(feat))) # requires feature matrix:feat x gene
    } else if (name == "glasso"){
        data <- apply(data, MARGIN = 2, FUN = function(x) { return(x - mean(x)) })
        cov_mat <- cov(data)
    } else if (name == "CS-CORE") {
    seurat_obj <- CreateSeuratObject(counts = t(data)) # requires gene x cell;
    CSCORE_result <- CSCORE(seurat_obj)
    # Obtain CS-CORE co-expression estimates
    CSCORE_coexp <- CSCORE_result$est
    # Obtain BH-adjusted p values
    CSCORE_p <- CSCORE_result$p_value
    p_matrix_BH = matrix(0, dim(t(data))[1], dim(t(data)))
    p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
    p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
    # Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
    CSCORE_coexp[p_matrix_BH > 0.05] <- 0
    cov_mat <- CSCORE_coexp
  }     
    return(cov_mat)
}


# input_data_foramt:cell*gene
# pbmc1 drop
raw_count <- read.csv('count.csv',sep='\t', header = 1, row.names = 1,check.names=F)
feat <- read.csv('features.csv',header = 1, row.names = 1,check.names=F)
# all_models_list <- c("pearson","spearman","minet","PIDC","GENIE3","scLink","glasso","GeneCluster")

# input data require cell*gene
count_t <- t(raw_count) 
genie3_cor <- comp_cor("GENIE3",count_t)
minet_cor <- comp_cor("minet",count_t)
pidc_cor <- comp_cor("PIDC",count_t)
pearson_cor <- comp_cor("pearson",count_t)
spearman_cor <- comp_cor("spearman",count_t)
sclink_cor <- comp_cor("scLink",count_t)
glasso_cor <- comp_cor("glasso",count_t)
GeneCluster_cor <- comp_cor("GeneCluster",feat)

count_t <- t(raw_count) 
cscore_cor <- comp_cor("CS-CORE",count_t)


# =====================================
#        Cor matrix heatmap
# =====================================
library(ComplexHeatmap)
library(circlize)
library(patchwork)
heatmap_plot = function(cor_matrix, col_fun, title){
    p <- Heatmap(as.matrix(cor_matrix), col = col_fun,
                cluster_rows = T,       # turn off row clustering
                cluster_columns = T,
                show_column_dend = FALSE,   # hide column dendrogram
                # show_row_dend = FALSE,   # hide row dendrogram
                column_title = title,
                column_title_gp = gpar(fontsize = 12,fontface = "bold"),
                show_row_names = F,   
                show_column_names = F,  
                name = "cor",
                # heatmap_legend_param = list(#title= "cor", 
                #                             title_position = "topcenter", 
                #                             # title_gp = gpar(fontsize = 8),
                #                             labels_gp = gpar(fontsize = 5),
                #                             legend_height=unit(2,"cm"), 
                #                             legend_direction="vertical"),
                width = 0.5, height = 0.5)
    return(as.ggplot(p))
}
