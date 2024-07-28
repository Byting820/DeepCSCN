#Description:GeneCluster co-expression performance eval.
#Four Datasets: Human PBMC1drop/PBMC1indrops/PBMC2drop/PBMC2indrops

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  library(stringr)
  source("BenchmarkUtils.r")
})

loadRealData <- function(cor_filename, net_filename) {
  data <- read.csv(cor_filename, header = 1, row.names = 1,check.names=F)
  network <- read.csv(net_filename, header = 1, row.names = 1,check.names=F)
  return(list(data = as.matrix(data), network = as.matrix(network)))
}


if (TRUE){
  defaultW <- getOption("warn")
  options(warn = -1)  
  # Parameters
  parser <- ArgumentParser()
  parser$add_argument("-d", "--data_name", default = 'pbmc2-Drop', help = "scRNA-seq data name.")
  parser$add_argument("-n", "--num_steps", default = '50', help = "Number of steps.")
  args <- parser$parse_args()
  species_type <- "PBMC"
  data_name <- args$data_name
  num_steps <- as.integer(args$num_steps)
  # Load data
  print(paste(rep("*", 70), collapse = ""))
  print(sprintf("[ %s | %s ] Loading data...", data_name, species_type))
  data_dir_path <- "cor_matrix"
  net_dir_path <- "Data/Benchmark_dataset/reference_net"
  data_filename <- sprintf("%s/%s500.csv", data_dir_path, data_name)
  net_filename <- sprintf("%s/%s-human_TF_similarity-sub_net_mat.csv", net_dir_path, data_name)
  data_obj <- loadRealData(data_filename, net_filename)
  data <- data_obj$data
  true_network <- data_obj$network
  TF_list <- rownames(true_network)
  TF_list <- intersect(TF_list,colnames(data))
  true_network <- true_network[TF_list,TF_list]
  print(sprintf("Data shape : #genes=%d", dim(data)[1]))

  # Model running
  model_name = "GeneCluster"
  start_time <- Sys.time() 
  print(paste(rep("=", 70), collapse = ""))
  print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
  
  estimate_type <- "thresholding"
  cov_mat <- data
  name <- "GeneCluster"
  cov_type <- "GeneCluster"
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
  est_res <- list(name = name, cov_type = cov_type, estimate_type = estimate_type, est_network_list = est_network_list, par_list = par_space)

  # Evaluation
  est_sub_network_list <- lapply(est_res$est_network_list, function (net){
    return (net[TF_list, TF_list])
  }) 

  metric_summary <- seqEstEvaluation(true_network, est_sub_network_list)
  metric_table <- metric_summary$metric
  metric_table["Pars"] <- unlist(est_res$par_list)
  metric_table["Model"] <- name
  # -----

  end_time <- Sys.time() 
  final_time <- end_time - start_time
  print(final_time)
  
  options(warn = defaultW)
  # Save data
  metric_save_path <- "Coexpress_Eva_Res"
  write.csv(apply(metric_table, 2, as.character), sprintf("%s/%s-metrics.csv", metric_save_path, data_name))
}


