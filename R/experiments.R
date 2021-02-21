source("R/config.R")
library(reticulate)
options(reticulate.conda_binary = conda_home)
use_condaenv('r-reticulate')
library(stringr)
library(tidyverse)
library(tensorflow)
tf$constant("Tensorflow test")
source("R/dataset_config.R")
source("R/experiments_config.R")
source("R/output_config.R")
source("R/utils.R")
source("R/data_constructor.R")
source("R/methods.R")
source("R/analysis.R")
source("R/output.R")
source("R/preprocess_data.R")
source("R/plot.R")
library(SingleCellExperiment)
##Wrapper
runExperiments <- function(experiment=c("simple_accuracy","cell_number", "sequencing_depth","celltype_structure",
                                        "batch_effects","inter_diseases","celltype_complexity","inter_species",
                                        "random_noise","inter_protocol")){
  switch(experiment,
         simple_accuracy = experiments.simple_accuracy(experiment),
         cell_number = experiments.cell_number(experiment),
         sequencing_depth = experiments.sequencing_depth(experiment),
         celltype_structure = experiments.celltype_structure(experiment),
         batch_effects = experiments.batch_effects(experiment),
         inter_diseases = experiments.inter_diseases(experiment),
         celltype_complexity = experiments.celltype_complexity(experiment),
         inter_species = experiments.inter_species(experiment),
         random_noise = experiments.random_noise(experiment),
         inter_protocol = experiments.inter_protocol(experiment),
         stop("Unkown experiments")
         )
}

###base function for assigning methods for all experiments
experiments.base.assign <- function(experiment, exp_config){
  require(stringr)
  require(readr)
  if(missing(exp_config)){
    stop("missing exp_config in base assign function")
  }
  print(str_glue("experiment {experiment} configuration: {exp_config}"))
  methods <- experiments.methods[[experiment]]
  
  ##assigning methods
  assign_methods <- methods$assign
  if_cv <- exp_config$cv
  
  if(if_cv){###if intra-dataset using cross-validation
    data <- utils.load_datasets(experiments.assign.data$train_dataset[[experiment]]) %>%
      constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
    colData(data)$unique_id <- 1:dim(colData(data))[[1]]
    assign_results <- experiments.run_assign(assign_methods,data,NA,exp_config)
    print("finish prediction for assign methods")
    return(list(assign_results=assign_results,train_data=data,test_data=NULL))
  }else{###if inter-dataset or fixed train/test dataset
    all_data <- utils.load_datasets(experiments.assign.data$train_dataset[[experiment]])
    if(experiment=="cell_number"){
      saved_train_data_path <- str_glue("{pretrained_home}/{experiment}_saved_train_data.RData")
      if(exp_config$test_sampling && exp_config$train_saved){
        load(saved_train_data_path)
      }else{
        train_data <- constructor.data_constructor(all_data, config=exp_config,experiment=experiment,if_train = TRUE)
        colData(train_data)$unique_id <- 1:dim(colData(train_data))[[1]]
        if(exp_config$test_sampling) save(train_data, file = saved_train_data_path)
      }
    }else{
      train_data <- constructor.data_constructor(all_data, config=exp_config,experiment=experiment,if_train = TRUE)
      colData(train_data)$unique_id <- 1:dim(colData(train_data))[[1]]
    }
    
    if(experiment=="cell_number"){
      saved_test_data_path <- str_glue("{pretrained_home}/{experiment}_saved_test_data.RData")
      if(exp_config$train_sampling&&exp_config$test_saved){
        load(saved_test_data_path)
      }else{
        if(experiments.assign.data$train_dataset[[experiment]]==experiments.assign.data$test_dataset[[experiment]]){
          #####fixed train/test dataset 
          print(str_glue('train_dataset = test_dataset'))
          test_sample_data <- all_data[,setdiff(colnames(all_data),colnames(train_data))]
          test_data <- constructor.data_constructor(test_sample_data, config=exp_config,experiment = experiment,if_train = FALSE)
        }else{#####interdataset
          print(str_glue('train_dataset != test_dataset'))
          test_data <- utils.load_datasets(experiments.assign.data$test_dataset[[experiment]]) %>%
            constructor.data_constructor(config=exp_config,experiment = experiment,if_train = FALSE)
        }
        colData(test_data)$unique_id <- (dim(colData(train_data))[[1]]+1):(dim(colData(train_data))[[1]]+dim(colData(test_data))[[1]])
        if(exp_config$train_sampling) save(test_data, file = saved_test_data_path)
      }
    }else{
      if(experiments.assign.data$train_dataset[[experiment]]==experiments.assign.data$test_dataset[[experiment]]){
          #####fixed train/test dataset
          print(str_glue('train_dataset = test_dataset'))
          test_sample_data <- all_data[,setdiff(colnames(all_data),colnames(train_data))]
          test_data <- constructor.data_constructor(test_sample_data, config=exp_config,experiment = experiment,if_train = FALSE)
      }else{#####interdataset
        print(str_glue('train_dataset != test_dataset'))
        test_data <- utils.load_datasets(experiments.assign.data$test_dataset[[experiment]]) %>%
          constructor.data_constructor(config=exp_config,experiment = experiment,if_train = FALSE)
      }
      colData(test_data)$unique_id <- (dim(colData(train_data))[[1]]+1):(dim(colData(train_data))[[1]]+dim(colData(test_data))[[1]])
    }
    assign_results <- experiments.run_assign(assign_methods,train_data,test_data,exp_config)
    print("finish prediction for assign methods")
    return(list(assign_results=assign_results,train_data=train_data,test_data=test_data))
  }
}

experiments.base.marker_gene_assign <- function(experiment, exp_config, data){
  require(stringr)
  require(readr)
  if(missing(exp_config)){
    exp_config <- experiments.parameters[[experiment]]
  }
  methods <- experiments.methods[[experiment]]
  ##clustering methods
  maker_gene_assign_methods <- methods$marker_gene_assign
  maker_gene_assign_results <- experiments.run_marker_gene_assign(maker_gene_assign_methods,data,exp_config)
  print("finish prediction for marker gene assign methods")
  maker_gene_assign_results
}
  
###base function for assigning methods for all experiments
experiments.base.cluster <- function(experiment, exp_config,data){
  require(stringr)
  require(readr)
  if(missing(exp_config)){
    exp_config <- experiments.parameters[[experiment]]
  }
  methods <- experiments.methods[[experiment]]
  ##clustering methods
  cluster_methods <- methods$cluster
  cluster_results <- experiments.run_cluster(cluster_methods,data,exp_config)
  print("finish prediction for cluster methods")
  cluster_results
}

###base function 
experiments.base.analyze <- function(assign_results,cluster_results,exp_config){
  print("start analyzing assigned/unassigned assigning/clustering methods")
  if(missing(exp_config)){
    exp_config <- experiments.parameters[[experiment]]
  }
  # methods <- experiments.methods[[experiment]]
  # cluster_methods <- methods$cluster
  # if(!purrr::is_null(methods$marker_gene_assign)){
  #   assign_methods <- c(methods$assign,methods$marker_gene_assign)
  # }else{
  #   assign_methods <- methods$assign
  # }
  cluster_methods <- colnames(cluster_results)[colnames(cluster_results)!="label"]
  assign_methods <- colnames(assign_results)[colnames(assign_results)!="label"]
  unlabeled_pctg_results <- analysis.run(assign_results,assign_methods,c("unlabeled_pctg"))
  results <- list(assign_results=assign_results,cluster_results=cluster_results)
  cluster_analysis_results <- analysis.run(results$cluster_results,cluster_methods,exp_config$metrics)
  cluster_analysis_results$supervised <- F
  assign_analysis_results <- analysis.run(results$assign_results,assign_methods,exp_config$metrics) %>%
    dplyr::bind_cols(unlabeled_pctg_results)
  assign_analysis_results$supervised <- T
  
  ####
  assigned_results <- utils.select_assigned(results)
  
  cluster_analysis_assigned_results <- analysis.run(assigned_results$cluster_results,cluster_methods,exp_config$metrics)
  cluster_analysis_assigned_results$assigned <- TRUE
  cluster_analysis_assigned_results$supervised <- F
  assign_analysis_assigned_results <- analysis.run(assigned_results$assign_results,assign_methods,exp_config$metrics)
  assign_analysis_assigned_results$assigned <- TRUE
  assign_analysis_assigned_results$supervised <- T
  
  # if(experiment %in% c('celltype_structure')){
  #   report_results <- list(all_assign_results=assign_analysis_results,
  #                          all_cluster_results=cluster_analysis_results,
  #                          assigned_assign_results=assign_analysis_assigned_results,
  #                          assigned_cluster_results=cluster_analysis_assigned_results)
  #   return(report_results)
  # }
  unassigned_results <- utils.select_unassigned(results)
  cluster_analysis_unassigned_results <- analysis.run(unassigned_results$cluster_results,cluster_methods,exp_config$metrics)
  cluster_analysis_unassigned_results$assigned <- FALSE
  cluster_analysis_unassigned_results$supervised <- F
  assign_analysis_unassigned_results <- analysis.run(unassigned_results$assign_results,assign_methods,exp_config$metrics)
  assign_analysis_unassigned_results$assigned <- FALSE
  assign_analysis_unassigned_results$supervised <- T
  report_results <- list(all_assign_results=assign_analysis_results,
                         all_cluster_results=cluster_analysis_results,
                         assigned_assign_results=assign_analysis_assigned_results,
                         assigned_cluster_results=cluster_analysis_assigned_results,
                         unassigned_assign_results=assign_analysis_unassigned_results,
                         unassigned_cluster_results=cluster_analysis_unassigned_results
                         )
  report_results
}

###base function for all experiments except batch effects
experiments.base <- function(experiment, exp_config){
  methods <- experiments.methods[[experiment]]
  assign_data_results <- experiments.base.assign(experiment,exp_config)
  assign_results <- assign_data_results$assign_results
  assign_data <- assign_data_results$train_data
  test_only <- !purrr::is_null(assign_data_results$test_data)
  if(test_only){
    test_samples <- colData(assign_data_results$test_data)$unique_id
    assign_data <- utils.append_sce(assign_data,assign_data_results$test_data)
  }
  if(length(methods$marker_gene_assign)>=1){
    marker_gene_assign_results <- experiments.base.marker_gene_assign(experiment,exp_config,assign_data)
    if(test_only){
      marker_gene_assign_results <- marker_gene_assign_results[test_samples,]
    }
    assign_results <- bind_cols(assign_results,dplyr::select(marker_gene_assign_results,-label))
  }
  assign_results <- utils.label_unassigned(assign_results,T)
  cluster_results <- experiments.base.cluster(experiment,exp_config,assign_data)
  if(test_only){
    cluster_results <- cluster_results[test_samples,]
  }
  cluster_results <- utils.label_unassigned(cluster_results,F)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  if(experiment %in% c("celltype_structure")){
    current_celltype_hierarchy <<- utils.createCellTypeHierarchy(assign_data,colData(assign_data)$label)
    current_celltype_weights <<- utils.createCellTypeWeights(assign_data,colData(assign_data)$label)
  }
  report_results <- experiments.base.analyze(assign_results,cluster_results)
  list(pred_results=combined_results,analy_results=report_results)
}


###simple accuracy experiment
experiments.simple_accuracy <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  base_results <- experiments.base(experiment,exp_config)
  report_results <- base_results$analy_results
  pred_results <- base_results$pred_results
  output.sink(experiment,pred_results,report_results)
  plot.plot(experiment,report_results,pred_results)
  report_results
}


####experiment with different number of cells per type
experiments.cell_number <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  cell_numbers <- exp_config$sample_num
  cell_pctgs <- exp_config$sample_pctg
  combined_cluster_results <- vector('list',length(cell_numbers))
  combined_assign_results <- vector('list',length(cell_numbers))
  combined_raw_results <- vector('list',length(cell_numbers))
  methods <- experiments.methods[[experiment]]
  sampling_by_pctg <- purrr::is_null(cell_numbers)||length(cell_numbers)==0
  cell_sampling_plan <- if(sampling_by_pctg) cell_pctgs else cell_numbers
  for(i in seq_along(cell_sampling_plan)){
    val <- cell_sampling_plan[[i]]
    print(str_glue('starting sampling :{val}'))
    config <- exp_config
    config$train_saved <- if(i==1) F else T
    config$test_saved <- if(i==1) F else T
    config$trained <- F
    if(!sampling_by_pctg){
      if(config$train_sampling){
        config$train_sample_num <- val
      }
      if(config$test_sampling){
        config$test_sample_num <- val
        config$trained <- if(i==1) F else T 
      }
    }else{
      if(config$train_sampling){
        config$train_sample_num <- NULL
        config$train_sample_pctg <- val
      }
      if(config$test_sampling){
        config$test_sample_num <- NULL
        config$test_sample_pctg <- val
        config$trained <- if(i==1) F else T 
      }
    }
    if(config$train_sampling&&!config$test_sampling){
      config$fixed_test <- T
    }else if(!exp_config$train_sampling&&exp_config$test_sampling){
      config$fixed_train <- T
    }
    base_results <- experiments.base(experiment,config)
    if(sampling_by_pctg){
      results <- base_results$analy_results%>% 
        purrr::map(~{.$sample_pctg <- val
        return(.)})
      raw_results <- base_results$pred_results
      raw_results$sample_pctg <- val
    }
    else{
      results <- base_results$analy_results%>% 
        purrr::map(~{.$sample_num <- val
                return(.)})
      raw_results <- base_results$pred_results
      raw_results$sample_num <- val
    }
    combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
    combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
    combined_raw_results[[i]] <- raw_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  combined_raw_results <- bind_rows(combined_raw_results)
  final_results <- bind_rows(combined_assign_results,combined_cluster_results)
  output.sink(experiment,combined_raw_results,final_results)
  plot.plot(experiment,final_results,combined_raw_results)
  final_results
}

###experiments with different sequencing depth
experiments.sequencing_depth <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  shallow_quant <- exp_config$quantile$low
  shallow_config <- exp_config
  shallow_config$quantile <- shallow_quant
  shallow_config$right <- FALSE
  shallow_results_base <- experiments.base(experiment,shallow_config)
  shallow_results <- shallow_results_base$analy_results %>%
    purrr::map(~{.$quantile<-shallow_quant 
            return(.)})
  shallow_raw_results <- shallow_results_base$pred_results
  shallow_raw_results$quantile <- shallow_quant
  
  deep_quant <- exp_config$quantile$high
  deep_config <- exp_config
  deep_config$quantile <- deep_quant
  deep_config$right <- TRUE
  deep_results_base <- experiments.base(experiment,deep_config)
  deep_results <- deep_results_base$analy_results %>%
    purrr::map(~{.$quantile<-deep_quant 
    return(.)})
  deep_raw_results <- deep_results_base$pred_results
  deep_raw_results$quantile <- deep_quant
  
  combined_assign_results <- bind_rows(bind_rows(shallow_results[grepl(".*_assign_.*",names(shallow_results))]),
                                       bind_rows(deep_results[grepl(".*_assign_.*",names(deep_results))]))
  combined_cluster_results <- bind_rows(bind_rows(shallow_results[grepl(".*_cluster_.*",names(shallow_results))]),
                                        bind_rows(deep_results[grepl(".*_cluster_.*",names(deep_results))]))
  combined_raw_results <- bind_rows(shallow_raw_results,deep_raw_results)
  final_results <- bind_rows(combined_assign_results,combined_cluster_results)
  output.sink(experiment,combined_raw_results,final_results)
  plot.plot(experiment,final_results,combined_raw_results)
  final_results
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
    
    results <- base_results$analy_results%>% 
      purrr::map(~{.$level<-level
      return(.)})
    raw_results <- base_results$pred_results
    raw_results$level <- level
    combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
    combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
    combined_raw_results[[i]] <- raw_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  combined_raw_results <- bind_rows(combined_raw_results)
  final_results <- bind_rows(combined_assign_results,combined_cluster_results)
  output.sink(experiment,combined_raw_results,final_results)
  plot.plot(experiment,final_results,combined_raw_results)
  final_results
}

###############
experiments.celltype_complexity <- function(experiment){
  
}

############
experiments.batch_effects <- function(experiment){
  require(gtools)
  exp_config <- experiments.parameters[[experiment]]
  # raw_data <- utils.load_datasets(experiments.data$batch_effects_no_free)
  methods <- experiments.methods[[experiment]]
  ###not remove batch effects
  # experiments.data[[experiment]] <<- purrr::map(raw_data,~{str_glue("{metadata(.)$study_name}_intersected.RDS")})
  # print(str_glue("after global experiments.data {experiment} is {experiments.data[[experiment]]}"))

  datasets_perm2 <- gtools::permutations(n=length(experiments.data$batch_effects_no_free),r=2,v=unlist(experiments.data$batch_effects_no_free),repeats.allowed = F)

  combined_cluster_results <- vector('list',dim(datasets_perm2)[[1]])
  combined_assign_results <- vector('list',dim(datasets_perm2)[[1]])
  combined_raw_results <- vector('list',dim(datasets_perm2)[[1]])
  intersected_datasets <- vector('list',dim(datasets_perm2)[[1]])
  for(i in 1:dim(datasets_perm2)[1]){
    # experiments.data[[experiment]] <<- unlist(datasets_comb2[,i])
    # print(str_glue("experiment data is {experiments.data[[experiment]]}"))
    data <- utils.load_datasets(datasets_perm2[i,])
    data_intersected <- list(str_glue("{metadata(data[[1]])$study_name}_{metadata(data[[2]])$study_name}_intersected.RDS"),
                             str_glue("{metadata(data[[2]])$study_name}_{metadata(data[[1]])$study_name}_intersected.RDS")
                            )
    intersected_datasets[[i]] <- data_intersected
    preprocess.intersect_sces(data,unlist(data_intersected))
    
    experiments.assign.data$train_dataset[[experiment]] <<- data_intersected[[1]]
    print(str_glue("assign train data is {experiments.assign.data$train_dataset[[experiment]]}"))
    experiments.assign.data$test_dataset[[experiment]] <<- data_intersected[[2]]
    print(str_glue("assign test data is {experiments.assign.data$test_dataset[[experiment]]}"))
    base_results <- experiments.base(experiment,exp_config)
    results <- base_results$analy_results%>% 
      purrr::map(~{.$train_dataset <- str_split(datasets_perm2[i,1],"\\.")[[1]][[1]]
      .$test_dataset <- str_split(datasets_perm2[i,2],"\\.")[[1]][[1]]
      .$batch_effects_removed <- FALSE
      return(.)}
      )
    raw_results <- base_results$pred_results
    raw_results$train_dataset <- experiments.assign.data$train_dataset[[experiment]]
    raw_results$test_dataset <- experiments.assign.data$test_dataset[[experiment]]
    combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
    combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
    combined_raw_results[[i]] <- raw_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  combined_raw_results <- bind_rows(combined_raw_results)
    
  ###remove batch effects from dataset and save

  if(exp_config$remove_batch){
    # datasets_comb2_no_be <- combn(experiments.data[[experiment]],2)
    combined_cluster_results_no_be <- vector('list',dim(datasets_perm2)[[1]])
    combined_assign_results_no_be <- vector('list',dim(datasets_perm2)[[1]])
    combined_raw_results_no_be <- vector('list',dim(datasets_perm2)[[1]])
    
    utils.update_batch_effects_free_config(experiment)
    
    for(i in 1:dim(datasets_perm2)[1]){
      # experiments.data[[experiment]] <<- unlist(datasets_comb2_no_be[,i])
      # print(str_glue("experiment data is {experiments.data[[experiment]]}"))
      
      data <- utils.load_datasets(unlist(intersected_datasets[[i]])) 
      data_no_be <- list(str_glue("{metadata(data[[1]])$study_name}_{metadata(data[[2]])$study_name}_batch_effects_removed.RDS"),
                         str_glue("{metadata(data[[2]])$study_name}_{metadata(data[[1]])$study_name}_batch_effects_removed.RDS")
                         )
      preprocess.remove_batch_effects(data,unlist(data_no_be))
      
      experiments.assign.data$train_dataset[[experiment]] <<- data_no_be[[1]]
      print(str_glue("assign train data no batch is {experiments.assign.data$train_dataset[[experiment]]}"))
      experiments.assign.data$test_dataset[[experiment]] <<- data_no_be[[2]]
      print(str_glue("assign test data no batch is {experiments.assign.data$test_dataset[[experiment]]}"))
      
      
      base_results_no_be <- experiments.base(experiment,exp_config)
      results_no_be <- base_results_no_be$analy_results%>% 
        purrr::map(~{.$train_dataset <- str_split(datasets_perm2[i,1],"\\.")[[1]][[1]]
        .$test_dataset <- str_split(datasets_perm2[i,2],"\\.")[[1]][[1]]
        .$batch_effects_removed <- T
        return(.)}
        )
      raw_results_no_be <- base_results_no_be$pred_results
      raw_results_no_be$train_dataset <- experiments.assign.data$train_dataset[[experiment]]
      raw_results_no_be$test_dataset <- experiments.assign.data$test_dataset[[experiment]]
      combined_assign_results_no_be[[i]] <- bind_rows(results_no_be[grepl(".*_assign_.*",names(results_no_be))])
      combined_cluster_results_no_be[[i]] <- bind_rows(results_no_be[grepl(".*cluster.*",names(results_no_be))])
      combined_raw_results_no_be[[i]] <- raw_results_no_be
    }
    
    combined_assign_results_no_be <- bind_rows(combined_assign_results_no_be)
    combined_cluster_results_no_be <- bind_rows(combined_cluster_results_no_be)
    combined_raw_results_no_be <- bind_rows(combined_raw_results_no_be)
  }
  if(!purrr::is_null(combined_assign_results_no_be)){
    
    total_combined_assign_results <- bind_rows(combined_assign_results, combined_assign_results_no_be)
    total_combined_cluster_results <- bind_rows(combined_cluster_results, combined_cluster_results_no_be)
    final_results <- bind_rows(total_combined_assign_results,total_combined_cluster_results)
    # report_results <- bind_rows(bind_rows(report_results),bind_rows(report_results_no_be))
    total_raw_results <- bind_rows(combined_raw_results,combined_raw_results_no_be)
  }else{
    final_results <- bind_rows(combined_assign_results,combined_cluster_results)
    total_raw_results <- combined_raw_results
  }
  output.sink(experiment,total_raw_results,final_results)
  plot.plot(experiment,final_results,total_raw_results)
  final_results
}

############
experiments.inter_diseases <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  datasets_perm2 <- gtools::permutations(n=length(experiments.data[[experiment]]),r=2,v=unlist(experiments.data[[experiment]]),repeats.allowed = F)
  combined_cluster_results <- vector('list',dim(datasets_perm2)[[1]])
  combined_assign_results <- vector('list',dim(datasets_perm2)[[1]])
  combined_raw_results <- vector('list',dim(datasets_perm2)[[1]])
  for(i in 1:dim(datasets_perm2)[[1]]){
    # experiments.data[[experiment]] <<- unlist(datasets_comb2[,i])
    print(str_glue("experiment data is {experiments.data[[experiment]]}"))
    experiments.assign.data$train_dataset[[experiment]] <<- unlist(datasets_perm2[i,])[1]
    print(str_glue("assign train data is {experiments.assign.data$train_dataset[[experiment]]}"))
    experiments.assign.data$test_dataset[[experiment]] <<- unlist(datasets_perm2[i,])[2]
    print(str_glue("assign test data is {experiments.assign.data$test_dataset[[experiment]]}"))
    base_results <- experiments.base(experiment,exp_config)
    results <- base_results$analy_results%>% 
      purrr::map(~{.$train_dataset <- str_split(datasets_perm2[i,1],"\\.")[[1]][[1]]
                   .$test_dataset <- str_split(datasets_perm2[i,2],"\\.")[[1]][[1]]
                   return(.)}
                 )
    raw_results <- base_results$pred_results
    raw_results$train_dataset <- experiments.assign.data$train_dataset[[experiment]]
    raw_results$test_dataset <- experiments.assign.data$test_dataset[[experiment]]
    combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
    combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
    combined_raw_results[[i]] <- raw_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  combined_raw_results <- bind_rows(combined_raw_results)
  final_results <- bind_rows(combined_assign_results,combined_cluster_results)
  output.sink(experiment,combined_raw_results,final_results)
  plot.plot(experiment,final_results,combined_raw_results)
  final_results
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

experiments.run_assign <- function(methods, train_data, test_data=NA, exp_config){
  if_cv <- exp_config$cv
  if(if_cv){
    print('start cross validation')
    results <- experiments.run_cv(methods,train_data,exp_config)
  }else{
    results <- data.frame(matrix(nrow=dim(colData(test_data))[[1]],ncol=length(methods)))
    colnames(results) <- methods
    results <- mutate(results,label=as.character(colData(test_data)$label))
    rownames(results) <- colData(test_data)$unique_id
    for(m in methods){
      print(str_glue('start assign method {m}'))
      m_result <- utils.try_catch_method_error(run_assign_methods(m,train_data,test_data,exp_config))
      if(inherits(m_result,"try-error")){
        results <- dplyr::select(results,-m)
      }else{
        if(is.factor(m_result)){
          m_result <- as.character(m_result)
        }
        results[[m]] <- m_result
      }
    }
  }
  results
}

experiments.run_marker_gene_assign <- function(methods,data,exp_config){
  results <- data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods)))
  colnames(results) <- methods
  results <- mutate(results,label=as.character(colData(data)$label))
  rownames(results) <- colData(data)$unique_id
  for(m in methods){
    print(str_glue('start marker gene assign method {m}'))
    m_result <- utils.try_catch_method_error(run_assign_methods(m,data,NULL,exp_config))
    if(inherits(m_result,"try-error")){
      results <- dplyr::select(results,-m)
    }else{
      if(is.factor(m_result)){
        m_result <- as.character(m_result)
      }
      results[[m]] <- m_result
    }
  }
  results
}

experiments.run_cluster <- function(methods,data,exp_config){
  results <- data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods)))
  colnames(results) <- methods
  results <- mutate(results,label=as.character(colData(data)$label))
  rownames(results) <- colData(data)$unique_id
  for(m in methods){
    print(str_glue('start cluster method {m}'))
    m_result <- utils.try_catch_method_error(run_cluster_methods(m,data))
    if(inherits(m_result,"try-error")){
      results <- dplyr::select(results,-m)
    }else{
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
  rownames(combined_results) <- colData(data)$unique_id
  for(i in seq_along(folds)){
    train_data <- data[,folds[[i]]]
    test_data <- data[,-folds[[i]]]
    for(m in methods){
      print(str_glue('start assign method {m} in {i}th fold of CV'))
      m_result <- utils.try_catch_method_error(run_assign_methods(m,train_data,test_data,exp_config))
      if(inherits(m_result,"try-error")){
        # combined_results[[m]][-folds[[i]]] <- NULL
        next
      }else{
        if(is.factor(m_result)){
          m_result <- as.character(m_result)
        }
        combined_results[[m]][-folds[[i]]] <- m_result
      }
    }
  }
  combined_results
}




