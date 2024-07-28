# construct cell-type-specific co-expression network

### =============================================================== ##
#    1.Calculate the correlation between the feature and the sample
### =============================================================== ##
library(dplyr)
library(circlize)
library(ComplexHeatmap)

corheatmap_plot <- function(plot_data,anno){
    # plot_data: Correlation matrix data used to draw heatmap
    col_anno = HeatmapAnnotation(celltype = anno,
                                col = list(celltype = c("B cell"="#1e77b4", "CD14+ monocyte"="#ed6a6c", "CD16+ monocyte"="#2da02b", "CD4+ T cell"="#ff9796", "Cytotoxic T cell"="#9167bd",
                                             "Dendritic cell" ="#bcbd22", "Megakaryocyte"="#c0c4c7",
                                             "Natural killer cell"="#bcbd22","Plasmacytoid dendritic cell"="#4a96d0"
                                            )))
    col_fun = colorRamp2(c(-0.2, 0, 0.2), c("#0084ff","white","#ff0404"))
    p <- Heatmap(plot_data, col = col_fun,
        cluster_rows = T,       # turn off row clustering
        cluster_columns = T,
        show_column_dend = FALSE,   # hide column dendrogram
        # column_title = "Heatmap of embedding importance",
        show_row_names = F, 
        show_column_names = F,
        name = "cor",
        # row_split = hclust_Res$x,
        column_split = sort_meta$CellType,
        # row_gap = unit(1, 'mm'),
        column_gap = unit(0.8, 'mm'),
        # row_km = 5,
        # column_km = 0,
        # row_title = NULL, 
        column_title = NULL,
        width = 0.85, height = 0.85,
        # right_annotation = row_anno,
        top_annotation = col_anno,
        )
    return(p)
}


### =============================================================== ##
#                 2.Screening cell type features
### =============================================================== ##
Get_top100feat <- function(FeatExprCor,feat){
    abs_FeatExprCor = abs(FeatExprCor)
    for (i in unique(sort_meta$CellType)){
        sample = sort_meta[sort_meta$CellType == i,]
        xcell_FeatExprCor = abs_FeatExprCor[,colnames(abs_FeatExprCor) %in% rownames(sample)]
        top100_feat = data.frame(cor=sort(rowMeans(xcell_FeatExprCor),decreasing=T)[1:100])
        xcell_feat = feat[,colnames(feat) %in% rownames(top100_feat)]
        write.csv(xcell_feat,paste("celltype_cor_heatmap/",i,"_feat.csv",sep=''))
    }
}


### =============================================================== ##
#                3.cell-type-specific co-expression
### =============================================================== ##
library(lsa)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

# Module cluster
cor_cluster <- function(cor_matrix,module_num){
    hclust_res <- cutree(hclust(as.dist(1-cor_matrix)),module_num)
    hclust_res <- as.data.frame(hclust_res)
    colnames(hclust_res) <- "cluster"
    sort_df <- arrange(hclust_res,cluster)
    return(sort_df)
}

heatmap_plot = function(cor_matrix,marker){
    # annotation info
    genelist = marker$gene
    index = which(rownames(cor_matrix) %in% genelist)
    labs = rownames(cor_matrix)[index]
    lab2 = rowAnnotation(foo = anno_mark(at = index,
                            labels = labs,
                            labels_gp = gpar(fontsize = 9),
                            lines_gp = gpar()))
    col_fun = colorRamp2(c(-1, 0, 1), c("#0084ff","white","#ff0404"))

    p <- Heatmap(cor_matrix, col = col_fun,
                cluster_rows = T,       # turn off row clustering
                cluster_columns = T,
                # show_column_dend = FALSE,   # hide column dendrogram
                # column_title = "sample",
                show_row_names = F, 
                show_column_names = F, 
                name = "Cor",
                width = 0.5, height = 0.5,
                right_annotation = lab2
                )
    return(p)
}

cor_cluster2 <- function(cor){
    dist_cor = 1-cor
    hclust_dist = hclust(as.dist(dist_cor), method = "average") 
    memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                        distM = dist_cor, 
                        deepSplit = 2,
                        pamRespectsDendro = FALSE,
                        minClusterSize = 10)
    names(memb) = colnames(cor)
    memb_df = as.data.frame(memb)
    colnames(memb_df) = "cluster"
    sort_df <- arrange(memb_df,cluster) 
    return(sort_df)
}

run_2cluster = function(cormat, clusterRes, module){
    # input:specific celltype cor matrix;
    # specific celltype first cluster res;
    # Select the module related to the cell type to be extracted (marker多的module)
    celltype_relative_gene = rownames(subset(clusterRes,cluster==module))
    celltype_relative_cormat = cormat[celltype_relative_gene,celltype_relative_gene]
    Cluster_res = cor_cluster2(celltype_relative_cormat)
    return(Cluster_res)
}

# GO
GO_analysis <- function(genelist,universe_gene){
    universe <- universe_gene
    GO <- enrichGO(genelist,
             OrgDb = org.Hs.eg.db,
             keyType = 'SYMBOL',
             ont = "ALL",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
            #  universe = universe
             )
    return(GO)
}

# Perform GO enrichment for each module
Run_GO <- function(Cluster_res,universe){
    df1 = data.frame()
    for (g in unique(Cluster_res$cluster)){
        module_gene = subset(Cluster_res,cluster==g)
        genelist = rownames(module_gene)
        # GO
        GOres <- GO_analysis(genelist,universe)  
        if (dim(GOres)[1] >=3 ){
            go_resdf = data.frame(GOres)
            top3go = go_resdf[order(go_resdf$p.adjust),][1:3,]
            top3go$module = rep(g, nrow(top3go))
            df1 = rbind(df1,top3go)            
        }
    }
    return(df1)
}

GO_heatmep = function(plot_df){
    my_palette <- colorRampPalette(colors = c("#7b67ed","white"))(length(bk <- seq(0,1e-3,by=1e-4)))
    p = pheatmap(t(plot_df),
                cellwidth = 20, cellheight = 20,
                cluster_rows=F,cluster_cols=F, 
                display_numbers = F, 
                number_color="white",
                fontsize_row=15,fontsize_col=15,
                angle_col=45,
                fontsize_number=15,
                color=my_palette,
                breaks=bk,
                legend_breaks=seq(0,1e-3,1e-4))
    return(p) 
}


if (TRUE){
    # dataset:pbmc1 drop-seq
    #1. Calculate the correlation between the feature and the sample
    count <- read.csv('count.csv',sep='\t',row.names=1)
    feat <- read.csv('features.csv',row.names=1,header=TRUE)
    meta_data <- read.table('meta_human.txt',sep='\t',header=TRUE)
    meta <- meta_data[((meta_data$Experiment=='pbmc1') & (meta_data$Method=='Drop-seq')),]
    rownames(meta) <- meta$NAME_TYPE
    meta$groupNo <- as.numeric(factor(meta$CellType))
    sort_meta <- arrange(meta,meta$groupNo)
    # celltype_list = c("B cell", "CD16+ monocyte", "CD4+ T cell", "Cytotoxic T cell",  "Natural killer cell")

    DEG = read.csv('DEG.csv',row.names=1)
    new_feat = feat[rownames(DEG),]
    new_count <- count[rownames(DEG), rownames(sort_meta)]
    FeatExprCor = cor(new_feat, new_count, use = "p")
    plot <- FeatExprCor[, rownames(sort_meta)]

    pdf("FeatureSampleCor_heatmap.pdf",width = 5,height = 3)
    anno = as.factor(sort_meta$CellType)
    p = corheatmap_plot(as.matrix(FeatExprCor), anno)
    draw(p)
    dev.off() 

    #2.Screening cell type features
    Get_top100feat(FeatExprCor,feat)

    #3.cell-type-specific co-expression
    celltype_feat = read.csv(celltype_feat,row.names=1)
    celltype_cor = cosine(as.matrix(t(celltype_feat))) #cor matrix
    markergene = read.csv('markerGene.csv',row.names=1)
    celltype_marker = markergene[markergene$cluster == "Bcell", ]
    pdf("celltype_cor_heatmap.pdf",width = 8,height = 8)
    cor_pic = heatmap_plot(celltype_cor,celltype_marker)
    draw(cor_pic)
    dev.off()

    module_num = 2
    Cluster_res = cor_cluster(celltype_cor,module_num)
    celltype_2Cluster_res = run_2cluster(celltype_cor,Cluster_res,2)

    GOres = Run_GO(celltype_2Cluster_res,rownames(celltype_2Cluster_res))
}


# compare to sclink and CS-CORE
library(scLink)
library(stringr)

Bcell_meta = meta[meta$CellType == "B cell",] 
Bcell_count = count[,rownames(Bcell_meta)] 
Bcell_top1000HEG <- Bcell_count[order(rowSums(Bcell_count),decreasing=T)[1:1000],]

# CS-CORE comp co-express
# https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html
library(devtools)
install_github("ChangSuBiostats/CS-CORE")
library(CSCORE)
library(Seurat)
library(WGCNA)
library(clusterProfiler)
# CSCORE input:raw is gene, col is cell
seurat_obj <- CreateSeuratObject(counts = Bcell_top1000HEG)
CSCORE_result <- CSCORE(seurat_obj)
# Obtain CS-CORE co-expression estimates
CSCORE_coexp <- CSCORE_result$est
# Obtain BH-adjusted p values
CSCORE_p <- CSCORE_result$p_value
p_matrix_BH = matrix(0, dim(Bcell_top1000HEG)[1], dim(Bcell_top1000HEG)[1])
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp[p_matrix_BH > 0.05] <- 0
# Compute the adjacency matrix based on the co-expression matrix
adj = WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 1)
# Compute the topological overlap matrix
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- rownames(Bcell_top1000HEG)
hclust_dist = hclust(as.dist(dissTOM), method = "average") 
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                     distM = dissTOM, 
                     deepSplit = 2,
                     pamRespectsDendro = FALSE,
                     minClusterSize = 10)
names(memb) = rownames(Bcell_top1000HEG)
# memb_tab <- table(memb)
memb_df = as.data.frame(memb)
colnames(memb_df) = "cluster"

# sclink
library(scLink)
networks = sclink_net(expr = t(Bcell_top1000HEG), ncores = 1, lda = seq(0.5, 0.1, -0.05))
networks$cor[1:3,1:3]
sclink_cor <- networks$cor
sclink_hclust = hclust(as.dist(1-sclink_cor), method = "average") 
clusters <- cutree(sclink_hclust, h = 0.85)
# memb_tab2 <- table(clusters)
memb3 = dynamicTreeCut::cutreeDynamic(dendro = sclink_hclust, 
                     distM = 1-sclink_cor, 
                     deepSplit = 2,
                     pamRespectsDendro = FALSE,
                     minClusterSize = 10)
names(memb3) = colnames(sclink_cor)
memb_df3 = as.data.frame(memb3)
colnames(memb_df3) = "cluster"

CSCORE_GOres = Run_GO(memb_df,rownames(memb_df))
sclink_GOres = Run_GO(memb_df3,rownames(memb_df3))

