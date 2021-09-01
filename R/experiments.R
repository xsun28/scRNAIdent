source("R/config.R")
library(reticulate)
options(reticulate.conda_binary = conda_home)
options(connectionObserver = NULL)
use_condaenv('r-reticulate')
library(stringr)
library(tidyverse)
library(tensorflow)
tf$constant("Tensorflow test")
source("R/dataset_config.R")
source("R/experiments_config.R")
source("R/experiments_analysis.R")
source("R/utils.R")
source("R/data_constructor.R")
source("R/methods.R")
source("R/analysis.R")
source("R/output.R")
source("R/preprocess_data.R")
source("R/plot.R")
source("R/summary.R")
library(SingleCellExperiment)
library(R.methodsS3)
logger <- utils.get_logger("DEBUG",log_file)
##Wrapper
runExperiments <- function(experiment=c("simple_accuracy","cell_number","celltype_number","sequencing_depth","celltype_structure",
                                        "batch_effects","sample_bias","celltype_complexity","inter_species",
                                        "random_noise","inter_protocol")){
  switch(experiment,
         simple_accuracy = experiments.simple_accuracy(experiment),
         cell_number = experiments.cell_number(experiment),
         celltype_number = experiments.celltype_number(experiment),
         sequencing_depth = experiments.sequencing_depth(experiment),
         celltype_structure = experiments.celltype_structure(experiment),
         batch_effects = experiments.batch_effects(experiment),
         sample_bias = experiments.sample_bias(experiment),
         celltype_complexity = experiments.celltype_complexity(experiment),
         inter_species = experiments.inter_species(experiment),
         random_noise = experiments.random_noise(experiment),
         inter_protocol = experiments.inter_protocol(experiment),
         imbalance_impacts = experiments.imbalance_impacts(experiment),
         unknown_types = experiments.unknown_types(experiment),
         stop("Unkown experiments")
         )
}

####base function for constructing sampled datasets for all methods
experiments.base.data_constructor <- function(experiment, exp_config){
  if_cv <- exp_config$cv
  if(if_cv){###if intra-dataset using cross-validation
    data <- utils.load_datasets(experiments.assign.data$train_dataset[[experiment]]) %>%
      constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
    # colData(data)$unique_id <- 1:dim(colData(data))[[1]]
    colData(data)$barcode <- colnames(data)
    colnames(data) <- 1:length(colnames(data)) ###make the sample id unique
    return(list(train_data=data,test_data=NULL))
  }
  ###if inter-dataset or fixed train/test dataset
  all_train_data <- utils.load_datasets(experiments.assign.data$train_dataset[[experiment]])
  if(!purrr::is_null(exp_config$fixed_train)&&exp_config$fixed_train){
    print("train sample is fixed")
    sample_seed <- if(purrr::is_null(exp_config$sample_seed)) 10 else exp_config$sample_seed
    train_data <- constructor.data_constructor(all_train_data, config=exp_config,experiment=experiment,if_train = TRUE, sample_seed)
  }else{
    print("train sample is not fixed")
    train_data <- constructor.data_constructor(all_train_data, config=exp_config,experiment=experiment,if_train = TRUE, sample_seed = exp_config$sample_seed)
  }
  # colData(train_data)$unique_id <- 1:dim(colData(train_data))[[1]]
  if(purrr::is_null(colData(train_data)$barcode)) colData(train_data)$barcode <- colnames(train_data)
  colnames(train_data) <- 1:length(colnames(train_data))
  
  ####construct test data
  if(!purrr::is_null(exp_config$fixed_test)&&exp_config$fixed_test){
    sample_seed <- if(purrr::is_null(exp_config$sample_seed)) 10 else exp_config$sample_seed
    print("test sample is fixed")
    if(experiments.assign.data$train_dataset[[experiment]]==experiments.assign.data$test_dataset[[experiment]]){
      print(str_glue('train_dataset = test_dataset'))
      test_data <- all_train_data[,setdiff(colnames(all_train_data),colnames(train_data))]
    }else{
      print(str_glue('train_dataset != test_dataset'))
      test_data <- utils.load_datasets(experiments.assign.data$test_dataset[[experiment]])
    }
    test_data <- constructor.data_constructor(test_data, config=exp_config, 
                                              experiment = experiment,if_train = FALSE, sample_seed)
  }else{
    print("test sample is not fixed")
    if(experiments.assign.data$train_dataset[[experiment]]==experiments.assign.data$test_dataset[[experiment]]){
      #####same train/test dataset 
      print(str_glue('train_dataset = test_dataset'))
      test_data <- all_train_data[,setdiff(colnames(all_train_data),colnames(train_data))]
    }else{#####interdataset
      print(str_glue('train_dataset != test_dataset'))
      test_data <- utils.load_datasets(experiments.assign.data$test_dataset[[experiment]]) 
    }
    test_data <- constructor.data_constructor(test_data, config=exp_config, experiment = experiment,if_train = FALSE,sample_seed = exp_config$sample_seed)
  }
  # colData(test_data)$unique_id <- (dim(colData(train_data))[[1]]+1):(dim(colData(train_data))[[1]]+dim(colData(test_data))[[1]])
  if(purrr::is_null(colData(test_data)$barcode)) colData(test_data)$barcode <- colnames(test_data)
  colnames(test_data) <- (dim(colData(train_data))[[1]]+1):(dim(colData(train_data))[[1]]+dim(colData(test_data))[[1]])    
  intersected_gene_train_test <- utils.intersect_on_genes(train_data,test_data)
  train_data <- intersected_gene_train_test$data1
  test_data <- intersected_gene_train_test$data2
  if(experiment %in% c("batch_effects") & exp_config$batch_free){
    batch_corrected_data <- preprocess.remove_batch_effects(list(train_data,test_data),exp_config$batch_correct_method)
    train_data <- batch_corrected_data[[1]]
    test_data <- batch_corrected_data[[2]]
  }
  return(list(train_data=train_data,test_data=test_data))
}


###base function for assigning methods for all experiments
experiments.base.assign <- function(experiment, train_dataset, test_dataset, exp_config){
  require(stringr)
  require(readr)
  if(missing(exp_config)){
    stop("missing exp_config in base assign function")
  }
  print(str_glue("experiment {experiment} configuration: {exp_config}"))
  methods <- experiments.methods[[experiment]]
  ##assigning methods
  assign_methods <- methods$assign
  assign_results <- experiments.run_assign(assign_methods,train_dataset,test_dataset,exp_config)
  print("finish prediction for assign methods")
  assign_results
}

  
###base function for assigning methods for all experiments
experiments.base.cluster <- function(experiment, exp_config,data){
  require(stringr)
  require(readr)
  if(missing(exp_config)){
    stop("missing exp_config in base cluster function")
  }
  methods <- experiments.methods[[experiment]]
  ##clustering methods
  cluster_methods <- methods$cluster
  cluster_results <- experiments.run_cluster(cluster_methods,data,exp_config)
  print("finish prediction for cluster methods")
  cluster_results
}


###base function for all experiments except batch effects
experiments.base <- function(experiment, exp_config){
  train_test_datasets <- experiments.base.data_constructor(experiment, exp_config)
  train_dataset <- train_test_datasets$train_data
  test_dataset <- train_test_datasets$test_data
  exp_config$have_test <- if(purrr::is_null(exp_config$have_test)) purrr::is_null(test_dataset) else exp_config$have_test
  if(exp_config$have_test){
    test_samples <- colnames(test_dataset)
    total_dataset <- utils.append_sce(train_dataset,test_dataset)
  }else{#### using cv and test_dataset=train_dataset
    total_dataset <- train_dataset
    test_dataset <- train_dataset
  }
  ####assign methods
  assign_results <- experiments.base.assign(experiment,train_dataset, test_dataset, exp_config)
  train_cell_types <- unique(colData(train_dataset)$label)
  assign_results <- utils.label_unassigned(train_cell_types,assign_results,T)
  ####cluster methods
  if(purrr::is_null(exp_config$fixed_test)||(!exp_config$fixed_test)||(!exp_config$clustered)){
    cluster_results <- if(experiment %in% c("batch_effects")) experiments.base.cluster(experiment,exp_config,total_dataset) else experiments.base.cluster(experiment,exp_config,test_dataset)
    if(experiment %in% c("batch_effects")) cluster_results <- cluster_results[test_samples,]
    cluster_results <- utils.label_unassigned(NULL,cluster_results,F)
  }else{
    print("test dataset is fixed and cluster results already generated")
    cluster_results <- NULL
  }
  return(list(experiment=experiment,assign_results=assign_results,cluster_results=cluster_results,exp_config=exp_config,train_dataset=train_dataset,test_dataset=test_dataset,total_dataset=total_dataset))
}

###simple accuracy experiment
experiments.simple_accuracy <- function(experiment){
  base_exp_config <- experiments.parameters[[experiment]]
  experiments.config.check_config(base_exp_config)
  for(i in seq_along(base_exp_config$intra_dataset)){
    train_dataset <- base_exp_config$intra_dataset[[i]]
    test_dataset <- base_exp_config$intra_dataset[[i]]
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init.simple_accuracy(train_dataset, test_dataset, base_exp_config)
    exp_config <- experiments.config.update(experiment, train_dataset, test_dataset,init_exp_config)
    base_results <- experiments.base(experiment,exp_config)
    report_results <- do.call(experiments.analysis,base_results)
    analy_results <- bind_rows(report_results$analy_results)
    pred_results <- bind_rows(report_results$pred_results)
    output.sink(experiment,pred_results,analy_results,exp_config)
    plot.plot(experiment,analy_results,pred_results,exp_config)
    utils.clean_marker_files()
    print(str_glue("{experiment} train_dataset={train_dataset}"))
    print(analy_results)
  }
  summarize_experiments(experiment,init_exp_config)
}


####experiment with different number of cells per type
experiments.cell_number <- function(experiment){
  base_exp_config <- experiments.parameters[[experiment]]
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$train_dataset
      test_dataset <- used_dataset[[j]]$test_dataset
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
    test_num <- init_exp_config$test_num
    combined_cluster_results <- vector('list',test_num)
    combined_assign_results <- vector('list',test_num)
    combined_raw_results <- vector('list',test_num)
    for(i in 1:test_num){
      exp_config <- init_exp_config
      exp_config$current_increment_index <- i-1
      exp_config <- experiments.config.update(experiment, train_dataset, test_dataset,exp_config)
      val <- if(exp_config$fixed_train) exp_config$target_test_num else exp_config$target_train_num
      print(str_glue('target sample num={val}'))
      base_results <- experiments.base(experiment,exp_config)
      if(exp_config$fixed_test){
        if(!purrr::is_null(base_results$cluster_results)){
          clustered_results <- base_results$cluster_results
        }else{
          base_results$cluster_results <- clustered_results
        }
      }
      report_results <- do.call(experiments.analysis,base_results)
      results <- report_results$analy_results%>% 
        purrr::map(~{
         .$sample_num <- val
         return(.)})
      raw_results <- report_results$pred_results
      raw_results$sample_num <- val
      combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
      combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
      combined_raw_results[[i]] <- raw_results
    }
    collapsed_combined_assign_results <- bind_rows(combined_assign_results)
    collapsed_combined_cluster_results <- bind_rows(combined_cluster_results)
    collapsed_combined_raw_results <- bind_rows(combined_raw_results)
    final_results <- bind_rows(collapsed_combined_assign_results,collapsed_combined_cluster_results)
    output.sink(experiment,collapsed_combined_raw_results,final_results,exp_config)
    plot.plot(experiment,final_results,collapsed_combined_raw_results,exp_config)
    utils.clean_marker_files()
    final_results
    print(str_glue("{experiment} train_dataset={train_dataset}, test_dataset={test_dataset}"))
    print(final_results)
  }
  summarize_experiments(experiment,init_exp_config)
}

###cell type number experiment



###experiments with different sequencing depth
experiments.sequencing_depth <- function(experiment){
  base_exp_config <- experiments.parameters[[experiment]]
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$train_dataset
      test_dataset <- used_dataset[[j]]$test_dataset
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
    test_num <- init_exp_config$test_num
    combined_cluster_results <- vector('list',test_num)
    combined_assign_results <- vector('list',test_num)
    combined_raw_results <- vector('list',test_num)
    for(i in 1:test_num){
      exp_config <- init_exp_config
      exp_config$current_increment_index <- i-1
      exp_config <-  experiments.config.update(experiment, train_dataset, test_dataset,exp_config)
      low_quant <- exp_config$low_quantile
      high_quant <- exp_config$high_quantile
      base_results <- experiments.base(experiment,exp_config)
      if(exp_config$fixed_test){
        if(!purrr::is_null(base_results$cluster_results)){
          clustered_results <- base_results$cluster_results
        }else{
          base_results$cluster_results <- clustered_results
        }
      }
      report_results <- do.call(experiments.analysis,base_results)
      results <- report_results$analy_results %>%
        purrr::map(~{
          .$quantile <- str_glue("{round(low_quant,2)}-{round(high_quant,2)}")
          return(.)})
      raw_results <- report_results$pred_results
      raw_results$quantile <- str_glue("{round(low_quant,2)}-{round(high_quant,2)}")
      combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
      combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
      combined_raw_results[[i]] <- raw_results
    }
    collapsed_combined_assign_results <- bind_rows(combined_assign_results)
    collapsed_combined_cluster_results <- bind_rows(combined_cluster_results)
    collapsed_combined_raw_results <- bind_rows(combined_raw_results)
    final_results <- bind_rows(combined_assign_results,collapsed_combined_cluster_results)
    output.sink(experiment,collapsed_combined_raw_results,final_results,exp_config)
    plot.plot(experiment,final_results,collapsed_combined_raw_results, exp_config)
    utils.clean_marker_files()
    print(str_glue("{experiment} train_dataset={train_dataset}, test_dataset={test_dataset}"))
    print(final_results)
  }
  summarize_experiments(experiment,init_exp_config)
}


#####experiments with different cell types
experiments.celltype_structure <- function(experiment){
  source("R/experiments_config.R") ####reset global environment to default value
  exp_config <- experiments.parameters[[experiment]]
  ####assigning
  assign_train_dataset <- experiments.assign.data$train_dataset[[experiment]]
  assign_test_dataset <- experiments.assign.data$test_dataset[[experiment]]
  data <- constructor.type_architecturer(exp_config, assign_train_dataset)
  levels <- colnames(colData(data))[grepl("Level*",colnames(colData(data)))]
  write_rds(data, str_glue("{data_home}{metadata(data)$study}_celltype_structure_dataset.RDS"))
  experiments.assign.data$train_dataset[[experiment]] <<- str_glue("{metadata(data)$study}_celltype_structure_dataset.RDS")
  if(assign_test_dataset!=assign_train_dataset){
    data <- constructor.type_architecturer(exp_config, assign_test_dataset)
    write_rds(data, str_glue("{data_home}{metadata(data)$study}_celltype_structure_test_dataset.RDS"))
    experiments.assign.data$test_dataset[[experiment]] <<- str_glue("{metadata(data)$study}_celltype_structure_test_dataset.RDS")
  }else{
    experiments.assign.data$test_dataset[[experiment]] <<- str_glue("{metadata(data)$study}_celltype_structure_dataset.RDS")
  }
  combined_cluster_results <- vector('list',length(levels))
  combined_assign_results <- vector('list',length(levels))
  combined_raw_results <- vector('list',length(levels))
  methods <- experiments.methods[[experiment]]
  
  for(i in seq_along(levels)){
    level <- levels[[i]]
    print(str_glue('starting cell type level:{level}'))
    config <- exp_config
    config$level <- level
    config$sample_num <- (length(levels)-i+1)*exp_config$sample_num
    base_results <- experiments.base(experiment,config) 
    current_celltype_hierarchy <<- utils.createCellTypeHierarchy(base_results$total_dataset,colData(base_results$total_dataset)$label)
    current_celltype_weights <<- utils.createCellTypeWeights(base_results$total_dataset,colData(base_results$total_dataset)$label)
    report_results <- do.call(experiments.analysis,base_results)
    
    if(assign_test_dataset!=assign_train_dataset){
      results <- report_results$analy_results%>% 
        purrr::map(~{
        .$level<-level
        return(.)})
    }else{
      results <- report_results$analy_results%>% 
        purrr::map(~{.$level<-level
        return(.)})
    }
    raw_results <- report_results$pred_results
    raw_results$level <- level
    combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
    combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
    combined_raw_results[[i]] <- raw_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  combined_raw_results <- bind_rows(combined_raw_results)
  final_results <- bind_rows(combined_assign_results,combined_cluster_results)
  output.sink(experiment,combined_raw_results,final_results,exp_config)
  plot.plot(experiment,final_results,combined_raw_results)
  utils.clean_marker_files()
  final_results
}

###############
experiments.celltype_complexity <- function(experiment){
  
}

############

experiments.batch_effects <- function(experiment){
  require(gtools)
  base_exp_config <- experiments.parameters[[experiment]]
  # datasets_perm2 <- gtools::permutations(n=length(experiments.data$batch_effects),r=2,v=unlist(experiments.data$batch_effects),repeats.allowed = F)
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$train_dataset
      test_dataset <- used_dataset[[j]]$test_dataset
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
    
    exp_config <-  experiments.config.update(experiment, train_dataset, test_dataset,init_exp_config)
    data <- utils.load_datasets(list(train_dataset,test_dataset))
    data_intersected <- list(str_glue("{metadata(data[[1]])$study_name}_{metadata(data[[2]])$study_name}_intersected.RDS"),
                             str_glue("{metadata(data[[2]])$study_name}_{metadata(data[[1]])$study_name}_intersected.RDS")
                              )
    preprocess.intersect_sces(data,unlist(data_intersected))
    
    experiments.config.update.train_test_datasets(experiment, data_intersected[[1]], data_intersected[[2]])
    exp_config$batch_free <- F
    base_results <- experiments.base(experiment,exp_config)
    report_results <- do.call(experiments.analysis,base_results)
    results <- report_results$analy_results%>% 
      purrr::map(~{.$batch_effects_removed <- FALSE
      return(.)}
      )
    raw_results <- report_results$pred_results
    # raw_results$train_dataset <- experiments.assign.data$train_dataset[[experiment]]
    # raw_results$test_dataset <- experiments.assign.data$test_dataset[[experiment]]
    raw_results$batch_effects <- T
    utils.clean_marker_files()
  
  ###remove batch effects from dataset and save
  
    if(exp_config$remove_batch){
      print("in batch effects experiment, starting removing batch effects")
      exp_config$batch_free <- T
      base_results_no_be <- experiments.base(experiment,exp_config)
      report_results_no_be <- do.call(experiments.analysis,base_results_no_be)
      results_no_be <- report_results_no_be$analy_results%>% 
        purrr::map(~{.$batch_effects_removed <- T
        return(.)}
        )
      raw_results_no_be <- report_results_no_be$pred_results
      raw_results_no_be$batch_effects <- F
      utils.clean_marker_files()
      utils.remove_files(unlist(data_intersected))
    }
    if(!purrr::is_null(results_no_be)){
      total_results <- bind_rows(bind_rows(results), bind_rows(results_no_be))
      total_raw_results <- bind_rows(raw_results,raw_results_no_be)
    }else{
      total_results <- bind_rows(results)
      total_raw_results <- raw_results
    }
    output.sink(experiment,total_raw_results,total_results,exp_config)
    plot.plot(experiment,total_results,total_raw_results,exp_config)
    utils.clean_marker_files()
    print(str_glue("finished batch effect dataset pair {experiments.assign.data$train_dataset[[experiment]]}-{experiments.assign.data$test_dataset[[experiment]]}"))
  }
  summarize_experiments(experiment,init_exp_config)
}


############
experiments.sample_bias <- function(experiment){
  require(gtools)
  base_exp_config <- experiments.parameters[[experiment]]
  # datasets_perm2 <- gtools::permutations(n=length(experiments.data[[experiment]]),r=2,v=unlist(experiments.data[[experiment]]),repeats.allowed = F)
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
  exp_config <- init_exp_config
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$dataset$train_dataset
      train_ind <- used_dataset[[j]]$inds$train_ind
      test_dataset <- used_dataset[[j]]$dataset$test_dataset
      test_ind <- used_dataset[[j]]$inds$test_ind
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    
    exp_config <-  experiments.config.update(experiment, train_dataset, test_dataset,exp_config,train_ind=train_ind,test_ind=test_ind)
    base_results <- experiments.base(experiment,exp_config)
    if(exp_config$fixed_test){
      if(!purrr::is_null(base_results$cluster_results)){
        exp_config$cluster_results[[test_dataset]] <- base_results$cluster_results
      }else{
        base_results$cluster_results <- exp_config$cluster_results[[test_dataset]]
      }
    }
    report_results <- do.call(experiments.analysis,base_results)
    results <- report_results$analy_results
    raw_results <- report_results$pred_results
    utils.clean_marker_files()
    results <- bind_rows(results)
    raw_results <- bind_rows(raw_results)
    output.sink(experiment,raw_results,results,exp_config)
    plot.plot(experiment,results,raw_results,exp_config)
    print(str_glue("finished sample-bias dataset pair {experiments.assign.data$train_dataset[[experiment]]}-{experiments.assign.data$test_dataset[[experiment]]}"))
  }
  summarize_experiments(experiment,init_exp_config)
}


###########
experiments.imbalance_impacts <- function(experiment){
  base_exp_config <- experiments.parameters[[experiment]]
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$train_dataset
      test_dataset <- used_dataset[[j]]$test_dataset
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
    type_pctgs <- init_exp_config$type_pctgs
    test_num <- init_exp_config$test_num
    combined_cluster_results <- vector('list',test_num)
    combined_assign_results <- vector('list',test_num)
    combined_raw_results <- vector('list',test_num)
    single_method_results <- vector('list',test_num)
    for(i in 1:test_num){
      exp_config <- init_exp_config
      exp_config$current_increment_index <- i-1
      exp_config <- experiments.config.update(experiment, train_dataset, test_dataset,exp_config)
      entropy <- round(utils.calc_entropy(type_pctgs[[i]]),2)
      base_results <- experiments.base(experiment,exp_config)
      if(exp_config$fixed_test){
        if(!purrr::is_null(base_results$cluster_results)){
          clustered_results <- base_results$cluster_results
        }else{
          base_results$cluster_results <- clustered_results
        }
      }
      report_results <- do.call(experiments.analysis,base_results)
      results <- report_results$analy_results %>%
        purrr::map(~{
          .$imbl_entropy <- entropy
          return(.)})
      raw_results <- report_results$pred_results
      raw_results$imbl_entropy <- entropy
      single_method_result <- report_results$single_method_result
      single_method_result$imbl_entropy <- entropy
      combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
      combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
      combined_raw_results[[i]] <- raw_results
      single_method_results[[i]] <- single_method_result
    }
    collapsed_combined_assign_results <- bind_rows(combined_assign_results)
    collapsed_combined_cluster_results <- bind_rows(combined_cluster_results)
    collapsed_combined_raw_results <- bind_rows(combined_raw_results)
    collapsed_combined_single_method_results <- bind_rows(single_method_results)
    final_results <- bind_rows(collapsed_combined_assign_results,collapsed_combined_cluster_results)
    output.sink(experiment,collapsed_combined_raw_results,final_results,exp_config,single_method_results=collapsed_combined_single_method_results)
    plot.plot(experiment,final_results,collapsed_combined_raw_results, exp_config,single_method_results=collapsed_combined_single_method_results)
    utils.clean_marker_files()
    print(str_glue("{experiment} train_dataset={train_dataset}, test_dataset={test_dataset}"))
    print(final_results)
  }
  summarize_experiments(experiment,init_exp_config)
}


###cell type number experiment

experiments.celltype_number <- function(experiment){
  base_exp_config <- experiments.parameters[[experiment]]
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$train_dataset
      test_dataset <- used_dataset[[j]]$test_dataset
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
    test_num <- init_exp_config$test_num
    combined_cluster_results <- vector('list',test_num)
    combined_assign_results <- vector('list',test_num)
    combined_raw_results <- vector('list',test_num)
    single_method_results <- vector('list',test_num)
    exp_config <- init_exp_config
    for(i in 1:test_num){
      exp_config$current_increment_index <- i
      exp_config <- experiments.config.update(experiment, train_dataset, test_dataset,exp_config)
      val <- if(exp_config$fixed_train) exp_config$test_number else exp_config$train_number   
      print(str_glue('cell type num={val}'))
      base_results <- experiments.base(experiment,exp_config)
      if(exp_config$fixed_test){
        exp_config$train_data <- base_results$train_dataset
        if(!purrr::is_null(base_results$cluster_results)){
          clustered_results <- base_results$cluster_results
        }else{
          base_results$cluster_results <- clustered_results
        }
      }else if(exp_config$fixed_train){
        exp_config$test_data <- base_results$test_dataset
      }
      report_results <- do.call(experiments.analysis,base_results)
      results <- report_results$analy_results%>% 
        purrr::map(~{
          .$type_num <- val
          return(.)})
      raw_results <- report_results$pred_results
      raw_results$type_num <- val
      single_method_result <- report_results$single_method_result
      single_method_result$type_num <- val
      
      combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
      combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
      combined_raw_results[[i]] <- raw_results
      single_method_results[[i]] <- single_method_result
      
    }
    collapsed_combined_assign_results <- bind_rows(combined_assign_results)
    collapsed_combined_cluster_results <- bind_rows(combined_cluster_results)
    collapsed_combined_raw_results <- bind_rows(combined_raw_results)
    collapsed_combined_single_method_results <- bind_rows(single_method_results)
    
    final_results <- bind_rows(collapsed_combined_assign_results,collapsed_combined_cluster_results)
    output.sink(experiment,collapsed_combined_raw_results,final_results,exp_config,single_method_results=collapsed_combined_single_method_results)
    plot.plot(experiment,final_results,collapsed_combined_raw_results,exp_config,single_method_results=collapsed_combined_single_method_results)
    
    utils.clean_marker_files()
    final_results
    print(str_glue("{experiment} train_dataset={train_dataset}, test_dataset={test_dataset}"))
    print(final_results)
  }
  summarize_experiments(experiment,init_exp_config)
}



##############
experiments.unknown_types <- function(experiment){
  base_exp_config <- experiments.parameters[[experiment]]
  experiments.config.check_config(base_exp_config)
  used_dataset <- if(base_exp_config$use_intra_dataset) base_exp_config$intra_dataset else base_exp_config$inter_dataset
  for(j in seq_along(used_dataset)){
    if(base_exp_config$use_intra_dataset){
      train_dataset <- used_dataset[[j]]
      test_dataset <- used_dataset[[j]]
    }else{
      train_dataset <- used_dataset[[j]]$train_dataset
      test_dataset <- used_dataset[[j]]$test_dataset
    }
    experiments.config.update.train_test_datasets(experiment, train_dataset, test_dataset)
    init_exp_config <- experiments.config.init(experiment, train_dataset, test_dataset,base_exp_config)
    
    test_num <- init_exp_config$test_num
    combined_cluster_results <- vector('list',test_num)
    combined_assign_results <- vector('list',test_num)
    combined_raw_results <- vector('list',test_num)
    single_method_results <- vector('list',test_num)
    for(i in 1:test_num){
      exp_config <- init_exp_config
      exp_config$current_increment_index <- i
      exp_config <- experiments.config.update(experiment, train_dataset, test_dataset,exp_config)
      print(str_glue('unknown type={exp_config$unknown_type}'))
      base_results <- experiments.base(experiment,exp_config)
      if(exp_config$fixed_test){
        if(!purrr::is_null(base_results$cluster_results)){
          clustered_results <- base_results$cluster_results
        }else{
          base_results$cluster_results <- clustered_results
        }
      }
      report_results <- do.call(experiments.analysis,base_results)
      results <- report_results$analy_results%>% 
        purrr::map(~{
          .$unknown_type <- exp_config$unknown_type
          return(.)})
      raw_results <- report_results$pred_results
      raw_results$unknown_type <- exp_config$unknown_type
      single_method_result <- report_results$single_method_result
      single_method_result$unknown_type <- exp_config$unknown_type
      combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
      combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
      combined_raw_results[[i]] <- raw_results
      single_method_results[[i]] <- single_method_result
    }
    collapsed_combined_assign_results <- bind_rows(combined_assign_results)
    collapsed_combined_cluster_results <- bind_rows(combined_cluster_results)
    collapsed_combined_raw_results <- bind_rows(combined_raw_results)
    collapsed_combined_single_method_results <- bind_rows(single_method_results)
    final_results <- bind_rows(collapsed_combined_assign_results,collapsed_combined_cluster_results)
    output.sink(experiment,collapsed_combined_raw_results,final_results,exp_config,single_method_results=collapsed_combined_single_method_results)
    plot.plot(experiment,final_results,collapsed_combined_raw_results, exp_config,single_method_results=collapsed_combined_single_method_results)
    utils.clean_marker_files()
    print(str_glue("{experiment} train_dataset={train_dataset}, test_dataset={test_dataset}"))
    print(final_results)
  }
  summarize_experiments(experiment,init_exp_config)
  
}





############
experiments.inter_species <- function(experiment){
  
}
############
experiments.random_noise <- function(experiment){
  
}
############
experiments.inter_protocol <- function(experiment){
  
}

#######

experiments.run_assign <- function(methods, train_data, test_data=NULL, exp_config){
  if_cv <- exp_config$cv
  if(if_cv){
    print('start cross validation')
    results <- experiments.run_cv(methods,train_data,exp_config)
  }else{
    results <- data.frame(matrix(nrow=dim(colData(test_data))[[1]],ncol=length(methods)))
    colnames(results) <- methods
    results <- mutate(results,label=as.character(colData(test_data)$label))
    # rownames(results) <- colData(test_data)$unique_id
    rownames(results) <- colnames(test_data)
    for(m in methods){
      print(str_glue('start assign method {m}'))
      m_result <- utils.try_catch_method_error(run_assign_methods(m,train_data,test_data,exp_config))
      if(inherits(m_result,"try-error")){
        print(str_glue("error occurs in {m}:{m_result}"))
        error(logger, str_glue("error occurs in train={experiments.assign.data$train_dataset[[experiment]]}, test={experiments.assign.data$test_dataset[[experiment]]},{m}:{m_result}"))
        # results <- dplyr::select(results,-m)
      }else{
        print(str_glue("{m} finished correctly"))
        if(is.factor(m_result)){
          m_result <- as.character(m_result)
        }
        results[[m]] <- m_result
      }
    }
  }
  results
}


experiments.run_cluster <- function(methods,data,exp_config){
  results <- data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods)))
  colnames(results) <- methods
  results <- mutate(results,label=as.character(colData(data)$label))
  # rownames(results) <- colData(data)$unique_id
  rownames(results) <- colnames(data)
  
  for(m in methods){
    print(str_glue('start cluster method {m}'))
    m_result <- utils.try_catch_method_error(run_cluster_methods(m,data,exp_config))
    if(inherits(m_result,"try-error")){
      print(str_glue("error occurs in {m}:{m_result}"))
      error(logger, str_glue("error occurs in train={experiments.assign.data$train_dataset[[experiment]]}, test={experiments.assign.data$test_dataset[[experiment]]}, {m}:{m_result}"))
      # results <- dplyr::select(results,-m)
    }else{
      print(str_glue("{m} finished correctly"))
      results[[m]] <- m_result
    }
  }
  results
}


###cross validation function
experiments.run_cv <- function(methods, data,exp_config){
  require(caret)
  require(tidyverse)
  stopifnot(is(data,"SingleCellExperiment"))
  cv_folds <- exp_config$cv_fold
  folds <- createFolds(factor(colData(data)$label),cv_folds,returnTrain = TRUE)
  combined_results <- data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods)))
  colnames(combined_results) <- methods
  combined_results <- mutate(combined_results,label=as.character(colData(data)$label))
  rownames(combined_results) <- colnames(data)
  for(i in seq_along(folds)){
    train_data <- data[,folds[[i]]]
    test_data <- data[,-folds[[i]]]
    for(m in methods){
      print(str_glue('start assign method {m} in {i}th fold of CV'))
      m_result <- utils.try_catch_method_error(run_assign_methods(m,train_data,test_data,exp_config))
      if(inherits(m_result,"try-error")){
        # combined_results[[m]][-folds[[i]]] <- NULL
        print(str_glue("error occurs in {m}, cv folds:{i}:{m_result}"))
        error(logger, str_glue("error occurs in {output.dataset_name[[experiment]]}, {m}, cv folds:{i}:{m_result}"))
        next
      }else{
        print(str_glue("{m} in cv folds:{i} finished correctly"))
        if(is.factor(m_result)){
          m_result <- as.character(m_result)
        }
        combined_results[[m]][-folds[[i]]] <- m_result
      }
    }
  }
  combined_results
}
