# Description: Benchmark of exisiting models on experimental datasets (with 1000HVGs).

suppressPackageStartupMessages({
  library(generics)
  library(purrr)
  library(Matrix)
  library(MASS) # for covariance estimation
  library(igraph) # network plotting
  library(cvms) # for evaluation metrics
  # library(prg) # precision-recall-gain curve
  # library(zoo) # rolling means
  library(argparse) # argument parsing
  library(caret) # confusion matrix

  library(QUIC) # solve GLasso model
  library(scLink)
  library(GENIE3)
  library(PIDC)
  library(minet)
  library(ZILGM)
  library(CSCORE)
  library(Seurat)

})


# =====================================
#        COVARIANCE COMPUTATION
# =====================================

computeCov <- function(data, type) {
  if (type == "pearson") {
    cov_mat <- cor(data, method = "pearson")
  } else if (type == "spearman") {
    cov_mat <- cor(data, method = "spearman")
  } else if (type == "covariance") {
    data <- apply(data, MARGIN = 2, FUN = function(x) { return(x - mean(x)) })
    cov_mat <- cov(data)
  } else if (type == "scLink") {
    # ========================================
    cov_mat <- sclink_cor(data, ncores = 1, nthre = 10, dthre = 0.9)
  } else if (type == "scDesign2") {
    cov_mat <- fit_Gaussian_copula(t(data), jitter = FALSE, zp_cutoff = 1.1)$cov_mat # requires gene x cell
  } else if (type == "GENIE3") {
    cov_mat <- GENIE3(t(data)) # requires gene x cell
  } else if (type == "PIDC") {
    gene_sparsity <- colSums(data)
    non_zero_idx <- which(gene_sparsity != 0)
    cov_mat <- matrix(0.0, dim(data)[2], dim(data)[2])
    row.names(cov_mat) <- colnames(cov_mat) <- colnames(data)
    diag(cov_mat) <- 1.0
    tmp_data <- data[, non_zero_idx]
    tmp_cov_mat <- PIDC(t(tmp_data), verbose = FALSE) # requires gene x cell
    cov_mat[non_zero_idx, non_zero_idx] <- tmp_cov_mat
  } else if (type == "minet") {
    cov_mat <- minet(data) # requires cell x gene; default parameters method="mrnet", estimator="spearman"
  } else if (type == "ZILGM") {
    # ZILGM not requires to explicitly compute the covariance, so just return the raw expression matrix
    return(data)
  } else if (type == "CS-CORE") {
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
  } else {
    stop(sprintf("Unkown covariance type %s!", type))
  }
  cov_mat[which(is.na(cov_mat))] <- 0.0
  return(cov_mat)
}


# =====================================
#           GRAPHICAL MODEL
# =====================================

thresholdingMat <- function(mat, thr) {
  # mat_dim <- dim(mat)
  # for (i in 1:(mat_dim[1] - 1)) {
  #   for (j in (i + 1):mat_dim[2]) {
  #     if (abs(mat[i, j]) < thr) {
  #       mat[i, j] <- 0.0
  #       mat[j, i] <- 0.0
  #     }
  #   }
  # }
  # return(mat)
  diag_val <- diag(mat) 
  mat <- as.matrix(forceSymmetric(mat)) 
  mat[abs(mat) < thr] <- 0.0
  diag(mat) <- diag_val
  return(mat)
}

estimateNet <- function(cov_mat, model_args, type) {
  if (type == "thresholding") {
    model_args$mat <- cov_mat
    est_network <- do.call("thresholdingMat", model_args)  
  } else if (type == "glasso") {
    model_args$S <- cov_mat
    model_args$msg <- 0
    est_network <- do.call("QUIC", model_args)$X
  } else if (type == "ZILGM") {
    # for ZILGM, varaible "cov_mat" is actually the raw expression matrix
    est_network <- zilgm(
      X = cov_mat, lambda = c(model_args), nlambda = 1,
      family = "NBII", update_type = "IRLS",
      do_boot = FALSE,
      sym = "OR", verbose = 2
    )
    est_network <- as.matrix(est_network$network[[1]])
    diag(est_network) <- 1.0
  } else {
    stop(sprintf("Unknown graphical model %s!", type))
  }
  est_network[which(is.na(est_network))] <- 0.0
  return(est_network)
}


# =====================================
#             EVALUATION
# =====================================

preprocessNet <- function(net) {
  net_triu <- net[upper.tri(net, diag = FALSE)] 
  adj_net_triu <- 1 * (net_triu != 0)
  sign_net_triu <- sign(net_triu) 
  return(list(original = net_triu, adj = adj_net_triu, sign = sign_net_triu))
}

networkMetric <- function(est_network, true_network) {
  # Pre-process networks
  metric_ind <- c("Balanced Accuracy", "F1", "Sensitivity", "Specificity", "Kappa", "MCC")
  true_object <- preprocessNet(true_network)
  est_object <- preprocessNet(est_network)
  # Metrics of adjacent matrix
  data_table <- as.data.frame(list(true = true_object$adj, prediction = est_object$adj))
  
  cvms::evaluate
  metrics <- evaluate(data_table, target_col = "true", prediction_cols = "prediction", type = "binomial") 
  confusion_mat <- metrics$"Confusion Matrix"[[1]]
  TN <- confusion_mat$N[1]
  FP <- confusion_mat$N[2]
  FN <- confusion_mat$N[3]
  TP <- confusion_mat$N[4]
  adj_metrics <- cbind(list(TN = TN, FP = FP, FN = FN, TP = TP), metrics[metric_ind])
  # Summary
  all_metrics <- adj_metrics
  return(list(
    metric = all_metrics,
    sparsity = mean(est_object$sign)
  ))
}


# =====================================
#           MODEL RUNNING
# =====================================

space.sampling <- function(model_name, upper, num_steps) {
  if (model_name == "thresholding") {
    lower <- 1e-3
    upper <- upper - 1e-3
    step <- (upper - lower) / num_steps 
    thr_seq <- seq(lower, upper, by = step)  
    random_thr_seq <- sample(num_steps)  
    model_args <- lapply(1:num_steps, function(i) { return(list(thr = thr_seq[random_thr_seq[i]])) }) 
  } else if (model_name == "glasso" | model_name == "quic") {
    lower <- 1e-3
    upper <- 1.0
    step <- (upper - lower) / num_steps
    rho_seq <- seq(lower, upper, by = step)
    random_rho_seq <- sample(num_steps)
    model_args <- lapply(1:num_steps, function(i) { return(list(rho = rho_seq[random_rho_seq[i]])) })
  } else if (model_name == "ZILGM"){
    lower <- 1e-3 * upper
    upper <- upper
    model_args <- seq(lower, upper, length.out = num_steps)
  } else {
    stop(sprintf("Unknown model name %s!", model_name))
  }
  return(model_args)
}


runModel <- function(name, data, num_steps) {
  # pearson, spearman, glasso, scLink, GENIE3, PIDC, scDesign2, minet
  if (name == "pearson") {
    cov_type <- "pearson"
    estimate_type <- "thresholding"
  } else if (name == "spearman") {
    cov_type <- "spearman"
    estimate_type <- "thresholding"
  } else if (name == "glasso") {
    cov_type <- "covariance"
    estimate_type <- "glasso"
  } else if (name == "scLink") {
    cov_type <- "scLink"
    estimate_type <- "glasso"
  } else if (name == "GENIE3") {
    cov_type <- "GENIE3"
    estimate_type <- "thresholding"
  } else if (name == "PIDC") {
    cov_type <- "PIDC"
    estimate_type <- "thresholding"
  } else if (name == "scDesign2") {
    cov_type <- "scDesign2"
    estimate_type <- "glasso"
  } else if (name == "minet") {
    cov_type <- "minet"
    estimate_type <- "thresholding"
  } else if (name == "ZILGM") {
    cov_type <- "ZILGM"
    estimate_type <- "ZILGM"
  } else if (name == "CS-CORE") {
    cov_type <- "CS-CORE"
    estimate_type <- "thresholding"
  } else {
    stop(sprintf("Unknown model name %s!", name))
  }
  # -----
  cov_mat <- computeCov(data, type = cov_type) 
  par_space <- space.sampling(estimate_type, max(abs(cov_mat)), num_steps) 
  est_network_list <- list()
  for (t in 1:length(par_space)) {
    if (t %% 5 == 0) { print(sprintf("%d iteration done", t)) } 
    step_pars <- par_space[[t]]  #[] extracts a list, [[]] extracts elements within the list
    est_network_list[[t]] <- estimateNet(cov_mat, model_args = step_pars, type = estimate_type)  
    if (estimate_type == "glasso"){
      rownames(est_network_list[[t]]) <- rownames(cov_mat)
      colnames(est_network_list[[t]]) <- rownames(cov_mat)
    } 
  }
  res <- list(name = name, cov_type = cov_type, estimate_type = estimate_type, est_network_list = est_network_list, par_list = par_space)
  return(res) 
}

seqEstEvaluation <- function(network, est_list) {
  metric_list <- list()
  conf_mat_list <- list()
  sparsity_list <- list()
  mse_list <- list()
  for (i in 1:range(length(est_list))) {
    iter_summary <- networkMetric(est_list[[i]], network)
    metric_list[[i]] <- iter_summary$metric
    conf_mat_list[[i]] <- iter_summary$conf_mat
    sparsity_list[[i]] <- iter_summary$sparsity
    mse_list[[i]] <- iter_summary$MSE
  }
  # -----
  metric_table <- as.data.frame(reduce(metric_list, rbind))
  metric_table["Sparsity"] <- sparsity_list
  metric_table["MSE"] <- mse_list
  return(list(metric = metric_table, conf_mat = conf_mat_list))
}

