# DeepCSCN

## Overview

DeepCSCN(Cell-type-specific co-expression network inference based on deep learning) is an R package designed for constructing cell-type-specific co-expression networks. It provides a series of functions for data preprocessing, generation of cell-type-specific features, identification of cell-type-specific modules, and perform GO enrichment analyses for specific cell types.  DeepCSCN offers a novel tool for researchers to dissect the intricate gene co-expression networks within distinct cell types. 

## Installation

To install the DeepCSCN package, you can use the following commands in R:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install DeepCSCN from GitHub
devtools::install_github("Byting820/DeepCSCN")
```

## Data preparation

- expr: Gene expression matrix
- feat: Gene feature matrix extracted by GeneCluster model(https://github.com/Byting820/GeneCluster)
- meta: Celltype meta info
- marker: Marker genes for different cell types


## Usage

### Load the package

```r
library(DeepCSCN)
library(clusterProfiler)
library(org.Hs.eg.db)

# Data preprocess
hESC = read.csv('data/hESC_Raw.csv',row.names=1, check.names=FALSE)
plot_df1 = data.frame(mean=rowMeans(hESC))
p1 <- expr_density_plot(plot_df1)
p <- p1$plot
cutoff <- p1$peak_x

meta_path <- "data/hESC_time_info.csv"
filter_expr <- filter_expr(hESC, meta_path, round(cutoff,2))

hvg <- findHVG_py(filter_expr, num_HVG = 1000,
                  python_env = "your_python_env")
dim(hvg)
write.csv(hvg, "data/hESC_1000HVG.csv")

# cell-type-spscific network
# 1. Identify Marker Genes
output_dir <- 'data/hESC_marker.csv'
markers <- marker_identify(hESC, meta_path, output_dir) 

# 2. cell type feature
meta <- read.csv(meta_path, row.names = 1)
feat <- read.csv('data/features.csv',sep='\t',row.names=1)
DEG <- markers[markers$gene %in% intersect(markers$gene,rownames(hvg)),]
CelltypeFeat(DEG, feat, hvg, meta, "time_info") 

# 3. cell type cormat and first cluster
# 00h
feat_00h <- read.csv("data/res/feat/00h_feat.csv",row.names=1)
res <- FirstCluster(feat_00h)
cormat_00h <- res$celltype_cor
clusterRes_00h <- res$celltype_cluster
write.csv(clusterRes_00h,"data/res/00h_FirstClusterRes.csv")

# cell-type-specific score
res2 <- CellTypeSpecificScore(hvg, meta, clusterRes_00h, "00h")
M1_score <- res2$M1_score
M2_score <- res2$M2_score

# specific cluster
specific_module <- clusterRes_00h[clusterRes_00h$cluster == 2, , drop=FALSE]
cormat_ <- cormat_00h[rownames(specific_module),rownames(specific_module)]
specific_clusterRes <- specific_cluster(cormat_,minClusterSize=50)
write.csv(specific_clusterRes,"data/res/00h_specific_clusterRes.csv")

# GO
GOres <- Run_GO(specific_clusterRes)

```


## Functions

## The DeepCSCN package provides the following key functions:

- CelltypeFeat(): Get celltype related feratures. 
- CellTypeSpecificScore(): Calculates cell-type-specific scores and identifies cell-type-specific modules.


## Contact
For questions, please contact: yutingya820@163.com
