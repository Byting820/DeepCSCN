# Main Analysis Script
# 
# This script performs main analysis including global and cell type-specific gene co-expression networks.
# 

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(Seurat)
library(DT)
library(WGCNA)
library(lsa)

# Source custom functions
# source("R/global_coexpression_net.R")
# source("R/cell-type-specific_net.R")

# Load data
expr_data_path <- 'data/processed/pbmc1-Drop-1000hvg.csv'
features_path <- 'data/processed/features.csv'

if (file.exists(expr_data_path) & file.exists(features_path)) {
    expr_data <- read.csv(expr_data_path, sep = '\t', row.names = 1, check.names = FALSE)
    features <- read.csv(features_path, row.names = 1)
} else {
    stop("Error: Required data files not found.")
}

#1. global network
global_res = global_net(features)
write.csv(global_res$global_cluster,"data/res/global_cluster_res.csv")

# Function to associate modules with cell types
clusres <- "data/res/global_cluster_res.csv"
marker <- "data/processed/markerGene.csv"
result <- CelltypeModuleMap(clusres, marker)
cross_table <- result$crosstable 
pvalue <- result$pvalue

pdf("plot/hist_enrichment.pdf", width = 25, height = 10)
p <- hist_enrichment_plot(t(pvalue), t(cross_table))
draw(p)   
dev.off()


#2. cell-type-spscific network
count <- read.csv('data/raw/pbmc1-Drop-1000hvg.csv', sep = '\t', row.names = 1)
feat <- read.csv('data/processed/features.csv', row.names = 1, header = TRUE)
meta_data <- read.table('data/processed/meta_human.txt', sep = '\t', header = TRUE)
meta <- meta_data[((meta_data$Experiment == 'pbmc1') & (meta_data$Method == 'Drop-seq')),]
rownames(meta) <- meta$NAME_TYPE
meta$groupNo <- as.numeric(factor(meta$CellType))
sort_meta <- arrange(meta, meta$groupNo)

# Screening cell type features
DEG_path <- 'data/processed/markerGene.csv'
CelltypeFeat(DEG_path, feat, count, sort_meta)
celltype_feat_path <- "data/res/B cell_feat.csv"
celltye_name = "B cell"
celltype_net = CelltypeNet(celltype_feat_path, celltye_name,count,meta)
write.csv(celltype_net,"data/res/celltype_net.csv")

GOres <- Run_GO(celltype_net)
GOres$celltype = rep(celltye_name, nrow(GOres))
# print(GOres)
# write.csv(GOres,"data/res/celltype_GO_res.csv")


# markergene <- read.csv('data/processed/markerGene.csv', row.names = 1)
# celltype_marker <- markergene[markergene$cluster == "Bcell", ]
# # Plot heatmap for cell type correlations
# if (nrow(celltype_cor) > 1) {
#     pdf("plot/celltype_cor_heatmap.pdf", width = 8, height = 8)
#     cor_pic <- heatmap_plot(celltype_cor, celltype_marker)
#     draw(cor_pic)
#     dev.off()
# } else {
#     message("Warning: Insufficient data for heatmap plot.")
# }

