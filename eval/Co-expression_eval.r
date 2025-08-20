#Description:GeneCluster co-expression performance eval.
#Datasets: BEELINE 7 datasets 1000hvg
#Author:Yuting Bai <yutingya820@163.com>

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  library(stringr)
  source("DeepCSCN/eval/BenchmarkUtils.r")
})

loadRealData <- function(cor_filename, net_filename) {
  data <- read.csv(cor_filename, header = 1, row.names = 1,check.names=F)
  network <- read.csv(net_filename, header = 1, row.names = 1,check.names=F)
  return(list(data = as.matrix(data), network = as.matrix(network)))
}

all_models_list <- c("scWGCNA","hdWGCNA","sclink","CSCORE","SingleCellGGM","GeneCluster")
dataset_name = 'hESC'

if (TRUE){
  defaultW <- getOption("warn")
  options(warn = -1)  #options环境设置参数 warn = -1 可以忽视任何警告
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = dataset_name, help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  args <- parser$parse_args()
  species_type <- dataset_name
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, species_type))

  # Model running
  model_metrics <- list()
  for (i in 1:length(all_models_list)) {
    model_name <- all_models_list[i]

    TFlist = read.csv("data/true_network/human-tfs.csv")
    data_dir_path <- "data/res/cormat"
    net_dir_path <- "data/true_network"
    data_filename <- sprintf("%s/%s/%s_cormat.csv", data_dir_path, dataset_name, model_name)
    #data_filename <- sprintf("%s/%s/Rawexpr_cormat.csv", data_dir_path, dataset_name, model_name)
    net_filename <- sprintf("%s/%s_net.csv", net_dir_path, data_name)
    data_obj <- loadRealData(data_filename, net_filename)
    data <- data_obj$data
    true_network <- data_obj$network
    true_network <- true_network[intersect(rownames(true_network),rownames(data)),intersect(colnames(true_network),colnames(data))]
    true_network <- true_network[intersect(rownames(true_network),TFlist$TF),intersect(colnames(true_network),TFlist$TF)]

    # TF_list <- rownames(true_network)
    # TF_list <- intersect(TF_list,colnames(data))
    # true_network <- true_network[TF_list,TF_list]
    print(sprintf("Data shape : #genes=%d", dim(data)[1]))

    start_time <- Sys.time() 

    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
  
    estimate_type <- "thresholding"
    cov_mat <- data
    name <- model_name
    cov_type <- "GeneCluster"
    par_space <- space.sampling(estimate_type, max(abs(cov_mat)), num_steps)
    est_network_list <- list()
    for (t in 1:length(par_space)) {
      if (t %% 5 == 0) { print(sprintf("%d iteration done", t)) } 
      step_pars <- par_space[[t]]  #[] extracts a list, [[]] extracts elements within the list
      est_network_list[[t]] <- estimateNet(cov_mat, model_args = step_pars, type = estimate_type)  
    }
    est_res <- list(name = name, cov_type = cov_type, estimate_type = estimate_type, est_network_list = est_network_list, par_list = par_space)

    # Evaluation
    est_sub_network_list <- lapply(est_res$est_network_list, function (net){
      return (net[rownames(true_network), colnames(true_network)])
    }) 

    metric_summary <- seqEstEvaluation(true_network, est_sub_network_list)
    metric_table <- metric_summary$metric
    metric_table["Pars"] <- unlist(est_res$par_list)
    metric_table["Model"] <- name
    # -----
    model_metrics[[i]] <- metric_table

    end_time <- Sys.time() # 记录终止时间
    final_time <- end_time - start_time
    print(final_time)
    print(sprintf("[ %s ] Total run time = [ %s ]", model_name, final_time))
  }

  all_model_metrics <- reduce(model_metrics, rbind)
  options(warn = defaultW)
  # Save data
  metric_save_path <- "data/res/eval_res"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s_CoexpressPrec.csv", metric_save_path, data_name))
}

