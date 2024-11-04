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

## Data preparation

- count: Preprocessed gene expression matrix
- feat: Gene feature matrix extracted by GeneCluster model(https://github.com/Byting820/GeneCluster)
- meta: Celltype meta info
- marker: Marker genes for different cell types


## Usage

### Load the package

```r
library(DeepCSCN)
# Load your data
count <- read.csv('data/processed/pbmc1-Drop-1000hvg.csv', sep = '\t', row.names = 1, check.names = FALSE)
feat <- read.csv('data/processed/features.csv', row.names = 1)
meta <- read.table('data/processed/meta_human.txt', row.names = 1, check.names = FALSE)

# 1. Generate global gene co-expression network
global_res <- global_net(feat)
print(global_res$sorted_clusters)

# 2. Associate modules with cell types
clusres <- "data/res/global_cluster_res.csv"
marker <- "data/processed/markerGene.csv"
result <- CelltypeModuleMap(clusres, marker)
pdf("plot/hist_enrichment.pdf", width = 25, height = 10)
p <- hist_enrichment_plot(t(result$pvalue), t(result$crosstable ))
draw(p)   
dev.off()

# 3. Generate cell-type-specific networks
celltype_name <- "B cell"
celltype_feat_path <- "data/res/B cell_feat.csv"
CelltypeFeat('data/processed/markerGene.csv', features, count, meta)
celltype_net <- CelltypeNet(celltype_feat_path, celltype_name, count, meta)
print(celltype_net)

```


## Functions

## The DeepCSCN package provides the following key functions:

- global_net(): Generates the global gene co-expression network.
- CelltypeModuleMap(): Associates gene modules with specific cell types.
- CelltypeFeat(): Get celltype related feratures. 
- CelltypeNet(): Calculates module-specific scores and identifies cell-type-specific modules.


## Contact
For questions, please contact: yutingya820@163.com
