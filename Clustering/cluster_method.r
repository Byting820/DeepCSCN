# differ cluster method evaluation
library(umap)
library(mclust)
library(ggplot2)
library(cluster)
library(reticulate)

my_mclust <- function(dim2_data,group1,group2){
    mclust_data = Mclust(umap_data$layout,G=group1:group2)
    return(mclust_data$classification)
}

mykmeans <- function(data,k){
    kmeans_res <- kmeans(data,centers=k)
    return(list(cluster = kmeans_res$cluster,center = kmeans_res$centers))
}

my_hclust <- function(data,method,k){
    hclust_res <- hclust(dist(data),method)
    return(cutree(hclust_res, k))
}

my_fanny <- function(data,k,metric,memb.exp){
    library(cluster)
    dist_matrix <- dist(data,method=metric)
    fanny_res <- fanny(dist_matrix,k,memb.exp=memb.exp)
    return(fanny_res$clustering)
}


scatter_plot <- function(plot_data,xlab,ylab,title){
    p = ggplot(data=plot_data,aes(x=umap1,y=umap2,colour=factor(clus_id)))+
            geom_point(size=1)+
            # scale_x_continuous(breaks=c(-5, 0, 5, 10))+
            # scale_y_continuous(breaks=c(-5, 0, 5, 10))+
            theme_bw()+
            labs(title=title,col="clus_id")+
            theme(legend.position="none",
                  axis.title = element_blank(),
                  axis.text = element_text(size = 16),
                  panel.border = element_rect(fill=NA,color="black"),
                  panel.grid=element_blank() 
                  )
    return(p)
}


expr_data <- read.csv('count.csv',sep='\t',header=1,row.names=1,check.names=F)
umap_data = umap(expr_data)

# Since the dataset has no real labels,
# in order to have the same number of clusters for different clustering methods,
# we first determine the number of clusters using the mclust clustering method,
# and the other clustering methods use the same number of clusters for clustering.

# mclust cluster
mclust_res <- my_mclust(umap_data$layout,5,40) 
plot_df <- data.frame(umap1=umap_data$layout[,1],umap2=umap_data$layout[,2],clus_id=mclust_res)
p1 = scatter_plot(plot_df,"UMAP1","UMAP2","Mclust")

# kmenas cluster
num_cluster = max(mclust_res)  #Num of clusters determined by mclust
kmeans_res <- mykmeans(umap_data$layout,num_cluster)
plot_df <- data.frame(umap1=umap_data$layout[,1],umap2=umap_data$layout[,2],clus_id=kmeans_res$cluster)
p2 = scatter_plot(plot_df,"UMAP1","UMAP2","KMeans")

# hclust cluster
hclust_res <- my_hclust(umap_data$layout,'complete',num_cluster)
plot_df <- data.frame(umap1=umap_data$layout[,1],umap2=umap_data$layout[,2],clus_id=hclust_res)
p3 = scatter_plot(plot_df,"UMAP1","UMAP2","Hclust")

# fanny cluster
fanny_res <- my_fanny(umap_data$layout, num_cluster, "euclidean", 1.15)
plot_df <- data.frame(umap1=umap_data$layout[,1],umap2=umap_data$layout[,2],clus_id=fanny_res)
p4 = scatter_plot(plot_df,"UMAP1","UMAP2","Fanny")


# deepcluster
features <- read.csv('features.csv',sep=',',row.names=1)
umap_feat = umap(features)
deep_res <- mykmeans(umap_feat$layout,num_cluster)
plot_df <- data.frame(umap1=umap_feat$layout[,1],umap2=umap_feat$layout[,2],clus_id=deep_res$cluster)
p5 = scatter_plot(plot_df,"UMAP1","UMAP2","GeneCluster")


library(patchwork)
p = (p2/p1)|(p3/p4)
ggsave("all_scatter.pdf",width = 6, height = 6)