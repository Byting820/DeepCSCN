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

