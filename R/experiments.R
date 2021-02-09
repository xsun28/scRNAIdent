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
    exp_config <- experiments.parameters[[experiment]]
  }
  print(str_glue("experiment {experiment} configuration: {exp_config}"))
  methods <- experiments.methods[[experiment]]
  
  ##assigning methods
  assign_methods <- methods$assign
  if_cv <- exp_config$cv
  
  if(if_cv){###if intra-dataset using cross-validation
    data <- utils.load_datasets(experiments.assign.data$train_dataset[[experiment]]) %>%
      constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
    assign_results <- experiments.run_assign(assign_methods,data,NA,exp_config)
  }else{###if inter-dataset 
    train_data <- utils.load_datasets(experiments.assign.data$train_dataset[[experiment]]) %>%
      constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
    test_data <- utils.load_datasets(experiments.assign.data$test_dataset[[experiment]]) %>%
      constructor.data_constructor(config=exp_config,experiment = experiment,if_train = FALSE)
    if(exp_config$use_test_only){
      data <- test_data ##using same dataset for clustering and marker gene based  below for other inter-dataset experiments
    }else{
      data <- utils.append_sce(train_data,test_data) ##using same dataset for clustering and marker gene based  below for batch effects
    }
    assign_results <- experiments.run_assign(assign_methods,train_data,test_data,exp_config)
  }
  print("finish prediction for assign methods")
  list(assign_results=assign_results,assign_data=data)
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
  # data <- utils.load_datasets(experiments.cluster.data[[experiment]]) %>%
  #   constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE,seed)
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
  assign_data <- assign_data_results$assign_data
  if(length(methods$marker_gene_assign)>=1){
    marker_gene_assign_results <- experiments.base.marker_gene_assign(experiment,exp_config,assign_data)
    assign_results <- bind_cols(assign_results,dplyr::select(marker_gene_assign_results,-label))
  }
  cluster_results <- experiments.base.cluster(experiment,exp_config,assign_data)
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
  combined_cluster_results <- vector('list',length(cell_numbers))
  combined_assign_results <- vector('list',length(cell_numbers))
  combined_raw_results <- vector('list',length(cell_numbers))
  methods <- experiments.methods[[experiment]]
  for(i in seq_along(cell_numbers)){
    num <- cell_numbers[[i]]
    print(str_glue('starting sample num:{num}'))
    config <- list(sample_num=num, cv=exp_config$cv, cv_fold=exp_config$cv_fold, metrics=exp_config$metrics)
    base_results <- experiments.base(experiment,config) 
    results <- base_results$analy_results%>% 
      purrr::map(~{.$sample_num<-num
              return(.)})
    raw_results <- base_results$pred_results
    raw_results$sample_num <- num
    combined_assign_results[[i]] <- bind_rows(results[grepl(".*_assign_.*",names(results))])
    combined_cluster_results[[i]] <- bind_rows(results[grepl(".*cluster.*",names(results))])
    combined_raw_results[[i]] <- raw_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  combined_raw_results <- bind_rows(combined_raw_results)
  final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
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
  final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
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
  final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
  output.sink(experiment,combined_raw_results,final_results)
  plot.plot(experiment,final_results,combined_raw_results)
  final_results
}

###############
experiments.celltype_complexity <- function(experiment){
  
}

############
experiments.inter_diseases <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  datasets_comb2 <- combn(experiments.cluster.data[[experiment]],2)
  combined_cluster_results <- vector('list',dim(datasets_comb2)[[2]])
  combined_assign_results <- vector('list',dim(datasets_comb2)[[2]])
  combined_raw_results <- vector('list',dim(datasets_comb2)[[2]])
  for(i in dim(datasets_comb2)[[2]]){
    experiments.cluster.data[[experiment]] <<- unlist(datasets_comb2[,i])
    print(str_glue("cluster data is {experiments.cluster.data[[experiment]]}"))
    experiments.assign.data$train_dataset[[experiment]] <<- unlist(datasets_comb2[,i])[1]
    print(str_glue("assign train data is {experiments.assign.data$train_dataset[[experiment]]}"))
    experiments.assign.data$test_dataset[[experiment]] <<- unlist(datasets_comb2[,i])[2]
    print(str_glue("assign test data is {experiments.assign.data$test_dataset[[experiment]]}"))
    base_results <- experiments.base(experiment,exp_config)
    results <- base_results$analy_results%>% 
      purrr::map(~{.$train_dataset <- experiments.assign.data$train_dataset[[experiment]]
                   .$test_dataset <- experiments.assign.data$test_dataset[[experiment]]
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
  final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
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

####experiments of inter-datasets with batch effects removed or not
experiments.batch_effects <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  raw_data <- utils.load_datasets(experiments.cluster.data$batch_effects_no_free)
  methods <- experiments.methods[[experiment]]
  ###not remove batch effects
  # print(str_glue("before global experiments.cluster.data {experiment} is {experiments.cluster.data[[experiment]]}"))
  
  experiments.cluster.data[[experiment]] <<- purrr::map(raw_data,~{str_glue("{metadata(.)$study_name}_intersected.RDS")})
  print(str_glue("after global experiments.cluster.data {experiment} is {experiments.cluster.data[[experiment]]}"))
  preprocess.intersect_sces(raw_data,experiments.cluster.data[[experiment]])
  train_datasets_combinations <- combn(experiments.cluster.data[[experiment]],length(experiments.cluster.data[[experiment]])-1)
  total_assign_results <- vector('list',dim(train_datasets_combinations)[2])
  assign_data <- NULL
  for(i in 1:dim(train_datasets_combinations)[2]){
    experiments.assign.data$train_dataset[[experiment]] <<- train_datasets_combinations[,i]
    print(str_glue('train dataset is {experiments.assign.data$train_dataset[[experiment]]}'))
    experiments.assign.data$test_dataset[[experiment]] <<- setdiff(experiments.cluster.data$batch_effects,train_datasets_combinations[,i]) 
    print(str_glue('test dataset is {experiments.assign.data$test_dataset[[experiment]]}'))
    assign_data_results <- experiments.base.assign(experiment,exp_config)
    assign_results <- assign_data_results$assign_results
    total_assign_results[[i]] <- assign_results
    if(purrr::is_null(assign_data)){
      assign_data <- assign_data_results$assign_data
    }
  }
  total_assign_results <- bind_rows(total_assign_results)
  if(length(methods$marker_gene_assign)>=1){
    marker_gene_assign_results <- experiments.base.marker_gene_assign(experiment,exp_config,assign_data)
    total_assign_results <- bind_cols(total_assign_results,dplyr::select(marker_gene_assign_results,-label))
  }
  cluster_results <- experiments.base.cluster(experiment,exp_config,assign_data)
  report_results <- experiments.base.analyze(total_assign_results,cluster_results) %>%
    purrr::map(~{ .$batch_effects_removed <- FALSE
      return(.)})
  total_raw_results <- bind_cols(total_assign_results,dplyr::select(cluster_results,-label))
  total_raw_results$batch_effects_removed <- FALSE
  
  ###remove batch effects from dataset and save
  report_results_no_be <- NULL
  if(exp_config$remove_batch){
    data <- utils.load_datasets(experiments.cluster.data[[experiment]]) 
    experiments.cluster.data[[experiment]] <<- purrr::map(data,~{str_glue("{metadata(.)$study_name}_batch_effects_removed.RDS")})
    preprocess.remove_batch_effects(data,experiments.cluster.data[[experiment]])
    train_datasets_combinations_no_be <- combn(experiments.cluster.data[[experiment]],length(experiments.cluster.data[[experiment]])-1)
    total_assign_results_no_be <- vector('list',dim(train_datasets_combinations_no_be)[2])
    assign_data_no_be <- NULL
    utils.update_batch_effects_free_config(experiment)
    for(i in 1:dim(train_datasets_combinations_no_be)[2]){
      experiments.assign.data$train_dataset[[experiment]] <<- train_datasets_combinations_no_be[,i]
      print(str_glue('train dataset is {experiments.assign.data$train_dataset[[experiment]]}'))
      experiments.assign.data$test_dataset[[experiment]] <<- setdiff(experiments.cluster.data[[experiment]],train_datasets_combinations_no_be[,i]) 
      print(str_glue('test dataset is {experiments.assign.data$test_dataset[[experiment]]}'))

      
      assign_data_results_no_be <- experiments.base.assign(experiment,exp_config)
      assign_results_no_be <- assign_data_results_no_be$assign_results
      total_assign_results_no_be[[i]] <- assign_results_no_be
      if(purrr::is_null(assign_data_no_be)){
        assign_data_no_be <- assign_data_results_no_be$assign_data
      }
    }
    total_assign_results_no_be <- bind_rows(total_assign_results_no_be)
    if(length(methods$marker_gene_assign_batch_free)>=1){
      experiments.methods[[experiment]]$marker_gene_assign <<- experiments.methods[[experiment]]$marker_gene_assign_batch_free
      marker_gene_assign_results_no_be <- experiments.base.marker_gene_assign(experiment,exp_config,assign_data_no_be)
      total_assign_results_no_be <- bind_cols(total_assign_results_no_be,dplyr::select(marker_gene_assign_results_no_be,-label))
    }
    cluster_results_no_be <- experiments.base.cluster(experiment,exp_config,assign_data_no_be)
    report_results_no_be <- experiments.base.analyze(total_assign_results_no_be,cluster_results_no_be) %>%
      purrr::map(~{ .$batch_effects_removed <- TRUE
      return(.)})
    total_raw_results_no_be <- bind_cols(total_assign_results_no_be,dplyr::select(cluster_results_no_be,-label))
    total_raw_results_no_be$batch_effects_removed <- TRUE
  }
  if(!purrr::is_null(report_results_no_be)){
    
    combined_assign_results <- bind_rows(bind_rows(report_results[grepl(".*_assign_.*",names(report_results))]),
                                         bind_rows(report_results_no_be[grepl(".*_assign_.*",names(report_results_no_be))]))
    combined_cluster_results <- bind_rows(bind_rows(report_results[grepl(".*_cluster_.*",names(report_results))]),
                                          bind_rows(report_results_no_be[grepl(".*_cluster_.*",names(report_results_no_be))]))
    final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
    # report_results <- bind_rows(bind_rows(report_results),bind_rows(report_results_no_be))
    total_raw_results <- bind_rows(total_raw_results,total_raw_results_no_be)
  }else{
    report_results <- bind_rows(report_results)
  }
  output.sink(experiment,total_raw_results,final_results)
  plot.plot(experiment,final_results,total_raw_results)
  final_results
}

#######

experiments.run_assign <- function(methods, train_data, test_data=NA, exp_config){
  if_cv <- exp_config$cv
  if(if_cv){
    print('start cross validation')
    results <- experiments.run_cv(methods,train_data,exp_config)
  }else{
    results <- as_tibble(data.frame(matrix(nrow=dim(colData(test_data))[[1]],ncol=length(methods))))
    colnames(results) <- methods
    results <- mutate(results,label=as.character(colData(test_data)$label))
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
  results <- as_tibble(data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods))))
  colnames(results) <- methods
  results <- mutate(results,label=as.character(colData(data)$label))
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
  results <- as_tibble(data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods))))
  colnames(results) <- methods
  results <- mutate(results,label=as.character(colData(data)$label))
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
  combined_results <- as_tibble(data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods))))
  colnames(combined_results) <- methods
  combined_results <- mutate(combined_results,label=as.character(colData(data)$label))
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




