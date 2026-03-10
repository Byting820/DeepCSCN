############################################################
# Function: overlap_marker
# Description: Calculate overlap between marker genes and 
#              clustering modules (from FirstCluster results)
############################################################
overlap_marker = function(marker_path, expr, FirstClusterRes_path, time){
    df = data.frame()
    # Load marker genes and keep only those present in expression matrix
    marker = read.csv(marker_path)
    marker = marker[marker$gene %in% rownames(expr), ]
    i_marker = marker[marker$cluster == time,]
    
    # Load clustering result
    clusterRes = read.csv(FirstClusterRes_path, row.names=1)
    module1 = subset(clusterRes,cluster==1)
    module2 = subset(clusterRes,cluster==2) 
    # Calculate overlap between markers and clustering modules
    module1_overlap_num = length(intersect(i_marker$gene, rownames(module1)))
    module2_overlap_num = length(intersect(i_marker$gene, rownames(module2)))

    tmp_df = data.frame(time = time, 
                        module1_overlap=module1_overlap_num, 
                        module2_overlap = module2_overlap_num)
    df = rbind(df, tmp_df)

    return(df)
}


############################################################
# Function: Comp_DifferTime_tau
# Description: Compute tissue-specificity index (Tau) across 
#              different time points for each gene
############################################################
Comp_DifferTime_tau = function(expr, meta_path){
    ## expr
    #           H9_00hb4s_001 H9_00hb4s_002 H9_00hb4s_003 H9_00hb4s_004
    # A1BG       2.510133      0.000000      1.702251      4.259372
    # A1CF       1.424896      1.208699      0.000000      0.000000
    ## meta
    #                   time_info
    # H9_00hb4s_001       00h
    # H9_00hb4s_002       00h
    library(dplyr)
    meta = read.csv(meta_path, row.names=1)
    merged_data <- merge(t(expr), meta, by = "row.names")
    rownames(merged_data) <- merged_data$Row.names
    merged_data$Row.names <- NULL
    # Calculate mean expression per time point
    mean_expression <- merged_data %>%
                        group_by(time_info) %>%
                        summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
                        column_to_rownames("time_info") %>%
                        t()
    
    df = data.frame()
    for (time in unique(meta$time_info)){
        time_meta = meta[meta$time_info == time, , drop=FALSE]
        time_path = paste("Data/time_clusterRes/",'/',time,'_clusterRes.csv',sep='')
        # Load clustering result at this time point
        time_clusterRes = read.csv(time_path, row.names=1, check.names=FALSE)
        time_expr = expr[rownames(time_clusterRes), rownames(time_meta)]

        tau_values <- numeric(nrow(time_expr))
        count=1
        # Compute tau for each gene
        for (i in rownames(time_expr)) {
            gene_expression <- as.numeric(as.character(mean_expression[i, ])) 
            max_expression <- max(gene_expression)
            normalized_values <- gene_expression / max_expression 
            # sum(1 - x_i)
            numerator <- sum(1 - normalized_values)
            tau <- numerator / (length(gene_expression) - 1)            
            
            tau_values[count] <- tau
            count = count+1
            if (count > nrow(time_expr)){
                break
            }
        }
        tmp_df = data.frame(time=time, gene=rownames(time_expr), tau=(tau_values))
        df = rbind(df,tmp_df)
    }
    return(df)
}

############################################################
# Function: stat_specific_num
# Description: Count number of time-specific and housekeeping 
#              genes based on tau threshold
############################################################
stat_specific_num = function(hHep_tau){
    hHep_tau_res = data.frame()

    for (i in unique(hHep_tau$time)){
        time_tau = hHep_tau[hHep_tau$time == i ,]
        # Define thresholds for specificity
        specific_gene_num = nrow(time_tau[time_tau$tau>=0.85,])
        home_gene_num = nrow(time_tau[time_tau$tau<=0.15,])
        
        df = data.frame(time = i, 
                        specific_gene_num = specific_gene_num, 
                        home_gene_num = home_gene_num)
        hHep_tau_res = rbind(hHep_tau_res, df)
    }
    return(hHep_tau_res)
}


############################################################
# Example usage
############################################################

# Input paths
marker_path = "data/hESC_marker.csv"
FirstClusterRes_00h_path = "data/res/00h_FirstClusterRes.csv"
expr = read.csv('data/hESC_1000HVG.csv', sep='\t', row.names=1, check.names=FALSE)

# Overlap between marker genes and cluster modules
overlap_marker_00h = overlap_marker(marker_path, expr, FirstClusterRes_00h_path, "00h")

# Tissue-specific index (Tau) across time points
hESC_DeepCSCN_tau <- Comp_DifferTime_tau(expr, meta_path)

# Count number of specific and housekeeping genes
hESC_DeepCSCN_tau_res <- stat_specific_num(hESC_DeepCSCN_tau)

