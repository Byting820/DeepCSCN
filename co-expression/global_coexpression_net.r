# construct gene co-expression network
library(ComplexHeatmap)
library(circlize)
library(dplyr)

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)

# =====================================
#        Cor matrix heatmap
# =====================================

heatmap_plot = function(cor_matrix){
    # annotation info
    col_fun = colorRamp2(c(0.2,0.6, 1), c("#0084ff","white","#ff0404"))

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


# =====================================
#        Module cluster
# =====================================

cor_cluster <- function(cor_matrix,module_num){
    hclust_res <- cutree(hclust(as.dist(1-cor_matrix)),module_num)
    hclust_res <- as.data.frame(hclust_res)
    colnames(hclust_res) <- "cluster"
    # sort_df <- arrange(hclust_res,cluster)
    return(hclust_res)
}

clus_eval <- function(clusRes,data){
    library(cluster)
    # comp silhouette
    silh <- silhouette(clusRes,dist(data))  
    silh_avg <- summary(silh)$si.summary[3] 
    return(silh_avg)
}


# ==========================================================================
#               Associate modules with cell types
#   (Identification of cell types based on the number of marker genes)
# ==========================================================================

marker_identify <- function(expr_path,metadata_path){
    library(Seurat)
    library(SeuratData)
    library(DT)

    expr_data <- read.table(expr_path)
    meta_data <- read.table(metadata_path,sep='\t',header=TRUE)
    meta <- meta_data[((meta_data$Experiment=='pbmc1') & (meta_data$Method=='Drop-seq')),] 
    rownames(meta) <- meta$NAME_TYPE
    new_expr_data <- subset(expr_data,select=c(meta$NAME_TYPE)) 

    #create seurat
     seurat_obj <- CreateSeuratObject(counts = new_expr_data, meta.data = meta) 
    # str(seurat_obj) 
    Idents(object = seurat_obj)  
    Idents(seurat_obj) <- meta$CellType 
    markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(markers,'markerGene.csv')
}


hist_enrichment_plot <- function(pvalue,cross_table){
    ha1 = HeatmapAnnotation(
        dist1 = anno_barplot(
            colSums(cross_table),
            bar_width = 1,
            gp = gpar(col = "white", fill = "#FFE200"),
            border = FALSE,
            axis_param = list(at = c(0, 20, 40, 60),
                labels = c("0",  "20", "40", "60")),
            height = unit(2.5, "cm")
        ), show_annotation_name = FALSE)

    ha2 = rowAnnotation(
        dist2 = anno_barplot(
            rowSums(cross_table),
            bar_width = 1,
            gp = gpar(col = "white", fill = "#FFE200"),
            border = FALSE,
            axis_param = list(at = c(0, 20, 40, 60, 80),
                labels = c("0", "20", "40", "60","80")),
            width = unit(2.5, "cm")
        ), show_annotation_name = FALSE)

    col_fun = circlize::colorRamp2(c(0,2,10), c("white", "yellow", "red"))
    p <- Heatmap(pvalue,
                name = "-log10P",
                col = col_fun,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                top_annotation = ha1,
                right_annotation = ha2,
                width = ncol(pvalue)*unit(8, "mm"),
                height = nrow(pvalue)*unit(8, "mm"),
                row_names_side = "left",
                column_labels = colnames(pvalue),
                column_names_rot = 0,
                #  row_order = hist_order,
                border = "black",
                heatmap_legend_param = list(
            title = "-log10P", at = c(0, 2, 4, 6, 8, 10)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (pvalue[i,j]>2){
                grid.text(sprintf("%d", cross_table[i, j]), x, y, gp = gpar(fontsize = 12))
            }
            grid.rect(x = x, y = y, width = width, height = height,
                gp = gpar(col = "grey", fill = NA))
    })
    return(p)
}

cor_cluster <- function(cor){
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
    return(memb_df)
}

if (TRUE){
    # Dataset: pbmc1 drop
    cor = read.csv('cor.csv',row.names=1,check.names=F)
    feat = read.csv('features.csv',row.names=1)

    clusterRes = cor_cluster(cor)
    sort_df <- arrange(clusterRes,cluster)
    plot_data = cor[rownames(sort_df),rownames(sort_df)]
    diag(plot_data) <- NA

    # Correlation Heatmap between modules
    Cluster_res = cor_cluster(cor,optimal_module)
    sort_df <- arrange(Cluster_res,cluster)
    loc <- match(rownames(sort_df),rownames(cor))
    plot_data = cor[loc,loc]
    diag(plot_data) <- NA

    col_anno = HeatmapAnnotation(module = as.factor(sort_df$cluster),
                                col = list(module=c("1"="#ffab66", "2"="#ed6a6c","3"="#cfdd7b", "4"="#e1abd1", "5"="#fac0bc",
                                        "6"="#aadfc4", "7"="#edd99f", "8"="#bbbbdfd5", "9"="#f0ffbb", "10"="#c0c4c7",
                                        "11"="#ffb01e", "12"="#ed7d1c", "13"="#dc717c","14"="#b39869", "15"="#4a96d0",
                                        "16"="#a2b3bf","17"="#c5b0d5","18"="#c49c94","19"="#e377c2","20"="#7f7f7f","21"="#c7c7c7",
                                        "22"="#bcbd22","23"="#17becf","24"="#9edae5","25"="#aec7e8")))  

    ha2 = rowAnnotation(module = as.factor(sort_df$cluster),
                        col = list(module=c("1"="#ffab66", "2"="#ed6a6c","3"="#cfdd7b", "4"="#e1abd1", "5"="#fac0bc",
                                "6"="#aadfc4", "7"="#edd99f", "8"="#bbbbdfd5", "9"="#f0ffbb", "10"="#c0c4c7",
                                "11"="#ffb01e", "12"="#ed7d1c", "13"="#dc717c","14"="#b39869", "15"="#4a96d0",
                                "16"="#a2b3bf","17"="#c5b0d5","18"="#c49c94","19"="#e377c2","20"="#7f7f7f","21"="#c7c7c7",
                                "22"="#bcbd22","23"="#17becf","24"="#9edae5","25"="#aec7e8")))
    # col_fun = colorRamp2(c(0, 0.5, 1), c("#0084ff","white","#ff0404"))
    pdf("GeneCluster_Heatmap.pdf", width=6, height=5)
    Heatmap(as.matrix(plot_data), col = gplots::colorpanel(250,'lemonchiffon',"orange"),
            cluster_rows = F,       # turn off row clustering
            cluster_columns = F,
            show_row_names = F, 
            show_column_names = F,
            # name = "cor",
            heatmap_legend_param = list(title= "cor", 
                            title_position = "topleft", 
                            at = c(0,0.5,1),
                            legend_height=unit(1.5,"cm"), 
                            legend_direction="vertical"),
            width = 0.5, height = 0.5,
            left_annotation = ha2,
            top_annotation = col_anno)
    dev.off()


    # Associate modules with cell types
    # 4 dataset:pbmc1drop/pbmc1indrops/pmbc2drop/pbmc2indrops
    enrich_propress <- function(clusterRes,marker){
        hclust_res <- read.csv(clusterRes)
        # hclust_res <- read.csv(clusterRes,row.names=1) # if wgcna res,run this code
        markergene <- read.csv(marker)
        colnames(hclust_res) <- c("gene","module")
        merge_gene <- merge(hclust_res,markergene,by="gene",all=F)
        cross_table <- table(merge_gene$module,merge_gene$cluster)
        pvalue <- cross_table

        for (i in rownames(cross_table)){
            for (j  in colnames(cross_table)){
                p <- phyper((cross_table[i,j]-1), 
                            sum(cross_table[,j]), 
                            (nrow(hclust_res)-sum(cross_table[,j])), 
                            sum(cross_table[i,]), 
                            lower.tail=F)
                q <- p.adjust(p, method="fdr")
                pvalue[i,j] <- q
            }
        }
        pvalue <- -1 * log10(pvalue)

        # pvalue <- t(pvalue)
        # cross_table <- t(cross_table)

        return(list(crosstable = cross_table, pvalue = pvalue))      
    }

    # pbmc1drop
    pbmc1drop_clusres = "Evaluate/03CoexpressionNet_contruct/module_cluster_res/new_pbmc1drop_GeneCluster_clusterRes.csv"
    pbmc1drop_marker = "Evaluate/03CoexpressionNet_contruct/data/pbmc1drop_markerGene.csv"
    result <- enrich_propress(pbmc1drop_clusres, pbmc1drop_marker)
    pbmc1drop_crosstable <- result$crosstable #(14,9)
    pbmc1drop_pvalue <- result$pvalue

    pdf("hist_enrichment.pdf", width=25, height=10)
    p = hist_enrichment_plot(t(pbmc1drop_pvalue), t(pbmc1drop_crosstable))
    draw(p)   
    dev.off()


    #############################################
    # Associate modules with cell types:method2
    library(WGCNA)
    expr_data = read.csv('count.csv',sep='\t',row.names=1,check.names=F)
    expr_data = t(expr_data)
    # meta data
    traitData <- read.table("meta_human.txt",sep='\t',header=TRUE)
    meta <- traitData[((traitData$Experiment=='pbmc1') & (traitData$Method=='Drop-seq')),]
    rownames(meta) <- meta$NAME_TYPE
    df <- table(meta$NAME_TYPE,meta$CellType)
    select_count = expr_data[rownames(df),]
    GeneCluster_module <- read.csv('GeneCluster_res.csv',row.names=1) 
    GeneCluster_module = as.matrix(GeneCluster_module)
    # GeneCluster_moduleColors <- labels2colors(GeneCluster_module)
    GeneCluster_MEs = moduleEigengenes(select_count,GeneCluster_module)$eigengenes
    GeneCluster_MET = orderMEs(GeneCluster_MEs)

    sort_traits <- df[match(rownames(GeneCluster_MET),rownames(df)),]
    design <- model.matrix(~0+sort_traits)
    colnames(design) <- colnames(df)
    moduleTraitCor = cor(GeneCluster_MET, design, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(GeneCluster_MET))
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    pdf("Module-trait_associations.pdf",width = 25, height=25)
    par(mar = c(18, 10, 5, 5))
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design), yLabels = names(GeneCluster_MET),
                    ySymbols = names(GeneCluster_MET), colorLabels = FALSE, colors = blueWhiteRed(50), 
                    textMatrix = textMatrix, setStdMargins = FALSE, 
                    cex.lab.x = 2.5,cex.lab.y = 2.0,
                    cex.text = 2.0, zlim = c(-1,1), 
                    main = paste("Module-trait relationships"))
    dev.off()
}

