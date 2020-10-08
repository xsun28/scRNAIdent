source("R/config.R")
source("R/experiments_config.R")
source("R/utils.R")
source("R/data_constructor.R")
source("R/methods.R")
source("R/analysis.R")
source("R/output.R")
##Wrapper
runExperiments <- function(experiment=c("simple_accuracy","cell_number", "sequencing_depth","cell_types", "batch_effects")){
  switch(experiment,
         simple_accuracy = experiments.simple_accuracy(experiment),
         cell_number = experiments.cell_number(experiment),
         sequencing_depth = experiments.equencing_depth(experiment),
         cell_types = experiments.cell_types(experiment),
         batch_effects = experiments.batch_effects(experiment),
         stop("Unkown experiments")
         )
}

###base function for all experiments
experiments.base <- function(experiment, exp_config){
  if(missing(exp_config)){
    exp_config <- experiments.parameters[[experiment]]
  }
  print(str_glue("experiment {experiment} configuration: {exp_config}"))
  methods <- experiments.methods[[experiment]]
  ##clustering methods
  cluster_methods <- methods$cluster
  data <- utils.load_dataset(experiments.cluster.data[[experiment]]) %>%
    constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
  cluster_results <- experiments.run_cluster(cluster_methods,data,exp_config)
  
  print("finish prediction for cluster methods, start analyzing")
  cluster_analysis_results <- analysis.run(cluster_results, cluster_methods,exp_config$metrics)
  ##assigning methods
  assign_methods <- methods$assign
  if_cv <- exp_config$cv
  if(if_cv){###if intra-dataset using cross-validation
    data <- utils.load_dataset(experiments.assign.data$train_dataset[[experiment]]) %>%
              constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
    assign_results <- experiments.run_asssign(assign_methods,data,NA,exp_config)
  }else{###if inter-dataset 
    train_data <- utils.load_dataset(experiments.assign.data$train_dataset[[experiment]]) %>%
              constructor.data_constructor(config=exp_config,experiment = experiment,if_train = TRUE)
    test_data <- utils.load_dataset(experiments.assign.data$test_dataset[[experiment]]) %>%
              constructor.data_constructor(config=exp_config,experiment = experiment,if_train = FALSE)
    assign_results <- experiments.run_asssign(assign_methods,train_data,test_data,exp_config)
  }
  print("finish prediction for assign methods, start analyzing")
  assign_analysis_results <- analysis.run(assign_results,assign_methods,exp_config$metrics)
  list(assign_results=assign_analysis_results,cluster_results=cluster_analysis_results)
}

###simple accuracy experiment
experiments.simple_accuracy <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  results <- experiments.base(experiment,exp_config)
  output.sink(experiment,results)
  results
}


####experiment with different number of cells per type
experiments.cell_number <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  cell_numbers <- exp_config$sample_num
  combined_cluster_results <- vector('list',length(cell_numbers))
  combined_assign_results <- vector('list',length(cell_numbers))
  methods <- experiments.methods[[experiment]]
  for(i in seq_along(cell_numbers)){
    num <- cell_numbers[[i]]
    print(str_glue('starting sample num:{num}'))
    config <- list(sample_num=num, cv=exp_config$cv, cv_fold=exp_config$cv_fold, metrics=exp_config$metrics)
    results <- experiments.base(experiment,config)
    cluster_results <- results$cluster_results
    assign_results <- results$assign_results
    cluster_results$sample_num <- num
    assign_results$sample_num <- num
    combined_assign_results[[i]] <- assign_results
    combined_cluster_results[[i]] <- cluster_results
  }
  combined_assign_results <- bind_rows(combined_assign_results)
  combined_cluster_results <- bind_rows(combined_cluster_results)
  final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
  output.sink(experiment,final_results)
  final_results
}

###experiments with different sequencing depth
experiments.sequencing_depth <- function(experiment){
  exp_config <- experiments.parameters[[experiment]]
  shallow_quant <- exp_config$quantile$low
  shallow_config <- exp_config
  shallow_config$quantile <- shallow_quant
  shallow_config$right <- FALSE
  shallow_results <- experiments.base(experiment,shallow_config)%>%
    map(~{.$quantile<-shallow_quant 
            return(.)})
  deep_quant <- exp_config$quantile$high
  deep_config <- exp_config
  deep_config$quantile <- deep_quant
  deep_config$right <- TRUE
  deep_results <- experiments.base(experiment,deep_config)%>%
    map(~{.$quantile<-deep_quant 
    return(.)})
  
  combined_assign_results <- bind_rows(shallow_results$assign_results,deep_results$assign_results)
  combined_cluster_results <- bind_rows(shallow_results$cluster_results,deep_results$cluster_results)
  final_results <- list(assign_results=combined_assign_results,cluster_results=combined_cluster_results)
  output.sink(experiment,final_results)
  final_results
}



#######

experiments.run_asssign <- function(methods, train_data, test_data=NA, exp_config){
  assign_results <- vector("list",length(methods))
  if_cv <- exp_config$cv
  if(if_cv){
    cv_folds <- exp_config$cv_fold
    print('start cross validation')
    results <- experiments.run_cv(methods,train_data,cv_folds)
  }else{
    results <- as_tibble(data.frame(matrix(nrow=dim(colData(test_data))[[1]],ncol=length(methods))))
    colnames(results) <- methods
    results <- mutate(results,label=colData(test_data)$label)
    for(m in methods){
      print(str_glue('start assign method {m}'))
      results[[m]] <- run_assign_methods(m,train_data,test_data)
    }
  }
  results
}


experiments.run_cluster <- function(methods,data,exp_config){
  results <- as_tibble(data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods))))
  colnames(results) <- methods
  results <- mutate(results,label=colData(data)$label)
  for(m in methods){
    print(str_glue('start cluster method {m}'))
    results[[m]] <- run_cluster_methods(m,data)
  }
  results
}


###cross validation function
experiments.run_cv <- function(methods, data,cv_folds){
  require(caret)
  require(tidyverse)
  stopifnot(is(data,"SingleCellExperiment"))
  folds <- createFolds(factor(colData(data)$label),cv_folds,returnTrain = TRUE)
  combined_results <- as_tibble(data.frame(matrix(nrow=dim(colData(data))[[1]],ncol=length(methods))))
  colnames(combined_results) <- methods
  combined_results <- mutate(combined_results,label=colData(data)$label)
  for(i in seq_along(folds)){
    train_data <- data[,folds[[i]]]
    test_data <- data[,-folds[[i]]]
    for(m in methods){
      print(str_glue('start assign method {m} in {i}th fold of CV'))
      pred_labels <- run_assign_methods(m,train_data,test_data)
      combined_results[[m]][-folds[[i]]] <- pred_labels
    }
  }
  combined_results
}




