# DeepCSCN

## Overview

DeepCSCN(Cell-type-specific co-expression network inference based on deep learning) is an R package designed for constructing global gene co-expression networks and cell-type-specific co-expression networks. It provides a series of functions for data preprocessing, identify co-expressed gene modules, and perform enrichment analyses for specific cell types.  DeepCSCN offers a novel tool for researchers to dissect the intricate gene co-expression networks within distinct cell types. 

## Installation

To install the DeepCSCN package, you can use the following commands in R:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install DeepCSCN from GitHub
devtools::install_github("Byting820/DeepCSCN")
```

## Usage

### Load the package

```r
library(DeepCSCN)
# Load your expression data and features
# features from GeneCluster model(https://github.com/Byting820/GeneCluster)
expr_data_path <- 'data/processed/pbmc1-Drop-1000hvg.csv'
features_path <- 'data/processed/features.csv'

expr_data <- read.csv(expr_data_path, sep = '\t', row.names = 1, check.names = FALSE)
features <- read.csv(features_path, row.names = 1)

# 1. Generate global gene co-expression network
global_res <- global_net(features)
global_cluster_res <- global_res$sorted_clusters
write.csv(global_cluster_res, "data/res/global_cluster_res.csv")

# 2. Associate modules with cell types
clusres <- "data/res/global_cluster_res.csv"
marker <- "data/processed/markerGene.csv"
result <- CelltypeModuleMap(clusres, marker)
pdf("plot/hist_enrichment.pdf", width = 25, height = 10)
p <- hist_enrichment_plot(t(result$pvalue), t(result$crosstable ))
draw(p)   
dev.off()

# 3. Analyze cell type-specific networks
count <- read.csv('data/raw/pbmc1-Drop-1000hvg.csv', sep = '\t', row.names = 1)
meta_data <- read.table('data/processed/meta_human.txt', sep = '\t', header = TRUE)
meta <- meta_data[(meta_data$Experiment == 'pbmc1') & (meta_data$Method == 'Drop-seq'),]
rownames(meta) <- meta$NAME_TYPE
# Perform analysis for a specific cell type
DEG_path <- 'data/processed/markerGene.csv'
celltype_feat_path <- "data/res/B cell_feat.csv"
celltype_name <- "B cell"
CelltypeFeat(DEG_path, features, count, meta)
celltype_net <- CelltypeNet(celltype_feat_path, celltype_name, count, meta)
write.csv(celltype_net, "data/res/celltype_net.csv")

```


## Functions

## The DeepCSCN package provides the following key functions:

- global_net(): Generates the global gene co-expression network.
- CelltypeModuleMap(): Associates gene modules with specific cell types.
- CelltypeNet(): Calculates module-specific scores and identifies cell type-specific networks.
- Run_GO(): Performs Gene Ontology enrichment analysis.


## Contact
For questions, please contact: yutingya820@163.com
