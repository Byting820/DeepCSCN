# Description: Benchmark of exisiting models on experimental datasets (with 1000HVGs).

suppressPackageStartupMessages({
  library(Matrix)
  library(argparse) # argument parsing
  library(stringr)
  source("BenchmarkUtils.r")
})

loadRealData <- function(data_filename, net_filename) {
  data <- read.csv(data_filename, header = 1, row.names = 1,check.names=F)
  gene_names <- colnames(data)
  gene_names <- str_split(gene_names, "_")
  gene_names <- unlist(lapply(gene_names, function(x){return (x[[2]])}))
  colnames(data) <- gene_names
  network <- read.csv(net_filename, header = 1, row.names = 1,check.names=F)
  return(list(data = as.matrix(data), network = as.matrix(network)))
}

# Test nine exisiting methods
all_models_list <- c("pearson","spearman","minet","PIDC","GENIE3","scLink","glasso")


# Model running on real experimental dataset
# Experimental dataset names:
#   Cortex1-10xChromium, Cortex2-10xChromium, Cortex2-Smart_seq2, Cortex2-Smart_seq2
#   pbmc1-drop, pbmc1-indrops, pbmc2-drop, pbmc2-indrops
if (TRUE) {
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
  data_dir_path <- "/Dataset"
  net_dir_path <- "/reference_net"
  data_filename <- sprintf("%s/%s300.csv", data_dir_path, data_name)
  net_filename <- sprintf("%s/%s-human_TF_similarity-sub_net_mat.csv", net_dir_path, data_name)
  data_obj <- loadRealData(data_filename, net_filename)
  data <- data_obj$data 
  true_network <- data_obj$network
  TF_list <- rownames(true_network)
  TF_list <- intersect(TF_list,colnames(data))
  true_network <- true_network[TF_list,TF_list]
  print(sprintf("Data shape : #cells=%d, #genes=%d", dim(data)[1], dim(data)[2]))
  # Model running
  model_metrics <- list()
  for (i in 1:range(length(all_models_list))) {
    start_time <- Sys.time() 

    model_name <- all_models_list[i]
    print(paste(rep("=", 70), collapse = ""))
    print(sprintf("[ %s ] Total num of steps = %d", model_name, num_steps))
    est_res <- runModel(name = model_name, data = data, num_steps = num_steps)
    
    # Evaluation
    est_sub_network_list <- lapply(est_res$est_network_list, function (net){
      return (net[TF_list, TF_list])
    }) 

    metric_summary <- seqEstEvaluation(true_network, est_sub_network_list)
    
    metric_table <- metric_summary$metric
    metric_table["Pars"] <- unlist(est_res$par_list)
    metric_table["Model"] <- model_name
    # -----
    model_metrics[[i]] <- metric_table

    end_time <- Sys.time()
    final_time <- end_time - start_time
    print(final_time)
    print(sprintf("[ %s ] Total run time = [ %s ]", model_name, final_time))
  }
  all_model_metrics <- reduce(model_metrics, rbind) 
  options(warn = defaultW)
  # Save data
  metric_save_path <- "data/resCoexpress_EvaRes"
  write.csv(apply(all_model_metrics, 2, as.character), sprintf("%s/%s-CSCORE-metrics300.csv", metric_save_path, data_name))
}
