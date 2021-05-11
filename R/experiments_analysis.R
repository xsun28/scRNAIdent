####Wrapper
experiments.analysis <- function(experiment, assign_results,cluster_results,exp_config,...){
  switch(experiment,
         simple_accuracy = experiments.analysis.simple_accuracy(assign_results,cluster_results,exp_config,...),
         cell_number = experiments.analysis.cell_number(assign_results,cluster_results,exp_config,...),
         celltype_number = experiments.analysis.celltype_number(assign_results,cluster_results,exp_config,...),
         sequencing_depth = experiments.analysis.sequencing_depth(assign_results,cluster_results,exp_config,...),
         celltype_structure = experiments.analysis.celltype_structure(assign_results,cluster_results,exp_config,...),
         batch_effects = experiments.analysis.batch_effects(assign_results,cluster_results,exp_config,...),
         inter_diseases = experiments.analysis.inter_diseases(assign_results,cluster_results,exp_config,...),
         celltype_complexity = experiments.analysis.celltype_complexity(assign_results,cluster_results,exp_config,...),
         inter_species = experiments.analysis.inter_species(assign_results,cluster_results,exp_config,...),
         random_noise = experiments.analysis.random_noise(assign_results,cluster_results,exp_config,...),
         inter_protocol = experiments.analysis.inter_protocol(assign_results,cluster_results,exp_config,...),
         imbalance_impacts = experiments.analysis.imbalance_impacts(assign_results,cluster_results,exp_config,...),
         unknown_types = experiments.analysis.unknown_types(assign_results,cluster_results,exp_config,...),
         stop("Wrong experiment")
  )
}


###base function 
experiments.analysis.base <- function(assign_results,cluster_results,exp_config){
  print("start analyzing assigned/unassigned assigning/clustering methods")
  if(missing(exp_config)){
    stop("missing exp_config in base analysis function")
  }

  cluster_methods <- colnames(cluster_results)[colnames(cluster_results)!="label"]
  assign_methods <- colnames(assign_results)[colnames(assign_results)!="label"]
  unlabeled_pctg_results <- analysis.run(assign_results,assign_methods,c("unlabeled_pctg"))
  supervised_cluster_num <- analysis.run(assign_results,assign_methods,c("cluster_num"))
  supervised_pred_type_max_pctg <- analysis.run(assign_results,assign_methods,c("pred_type_max_pctg"))
  assign_analysis_results <- analysis.run(assign_results,assign_methods,exp_config$metrics) %>%
    dplyr::bind_cols(unlabeled_pctg_results) %>% dplyr::bind_cols(supervised_cluster_num) %>% 
    dplyr::bind_cols(supervised_pred_type_max_pctg)
  assign_analysis_results$supervised <- T
  
  unsupervised_cluster_num <- analysis.run(cluster_results,cluster_methods,c("cluster_num"))
  unsupervised_pred_type_max_pctg <- analysis.run(cluster_results,cluster_methods,c("pred_type_max_pctg"))
  cluster_analysis_results <- analysis.run(cluster_results,cluster_methods,exp_config$metrics) %>%
    dplyr::bind_cols(unsupervised_cluster_num) %>% 
    dplyr::bind_cols(unsupervised_pred_type_max_pctg)
  cluster_analysis_results$supervised <- F
  
  
  
  results <- list(assign_results=assign_results,cluster_results=cluster_results)
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

experiments.analysis.simple_accuracy <- function(assign_results,cluster_results,exp_config,...){
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  experiments.analysis.attach_dataset_props("simple_accuracy", exp_config,report_results=report_results,
                                            combined_results=combined_results,...)
}

experiments.analysis.cell_number <- function(assign_results,cluster_results,exp_config,...){
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  experiments.analysis.attach_dataset_props("cell_number", exp_config,report_results=report_results,
                                            combined_results=combined_results,...)
}

experiments.analysis.celltype_number <- function(assign_results,cluster_results,exp_config,...){
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  experiments.analysis.attach_dataset_props("celltype_number", exp_config,report_results=report_results,
                                            combined_results=combined_results,...)
}

experiments.analysis.sequencing_depth <- function(assign_results,cluster_results,exp_config,...){
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  experiments.analysis.attach_dataset_props("sequencing_depth", exp_config,report_results=report_results,
                                            combined_results=combined_results,...)
}

experiments.analysis.batch_effects <- function(assign_results,cluster_results,exp_config,...){
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  experiments.analysis.attach_dataset_props("batch_effects", exp_config,report_results=report_results,
                                            combined_results=combined_results,...)
}

experiments.analysis.inter_diseases <- function(assign_results,cluster_results,exp_config,...){
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  experiments.analysis.attach_dataset_props("inter_diseases", exp_config,report_results=report_results,
                                            combined_results=combined_results,...)
}

experiments.analysis.imbalance_impacts <- function(assign_results,cluster_results,exp_config,...){
  require(tidyverse)
  extra_args <- list(...)
  train_dataset <- extra_args$train_dataset
  test_dataset <- extra_args$test_dataset
  
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  if(exp_config$fixed_train){
    type_pctg <- dplyr::group_by(as.data.frame(colData(test_dataset)), label) %>% dplyr::summarize(type_pctg=round(n()/length(colnames(test_dataset)),2))
  }else{
    type_pctg <- dplyr::group_by(as.data.frame(colData(train_dataset)), label) %>% dplyr::summarize(type_pctg=round(n()/length(colnames(train_dataset)),2))
  }
  cell_types <- type_pctg$label
  methods <- colnames(combined_results)[colnames(combined_results)!="label"]
  supervised_methods <- colnames(assign_results)[colnames(assign_results)!="label"]
  unsupervised_methods <- colnames(cluster_results)[colnames(cluster_results)!="label"]
  single_method_result <- bind_rows(purrr::map(methods, function(method){
                                                                      if(method %in% unsupervised_methods){
                                                                        type_accuracy_f1 <- bind_rows(purrr::map(cell_types, function(type){ method_pred_true <- dplyr::select(combined_results,method,label)
                                                                                                                                              if(all(sapply(method_pred_true[[method]],function(x){is_null(x)||is.na(x)}))){
                                                                                                                                                return(list(f1=NA,acc=NA))
                                                                                                                                              }
                                                                                                                                                type_pred_true <- dplyr::filter(method_pred_true,label==type)
                                                                                                                                                unique_clusters <- unique(type_pred_true[[method]])                                                
                                                                                                                                                cluster_f_betas <- purrr::map(unique_clusters,function(cluster_num){ 
                                                                                                                                                analysis.cluster.fbeta(method_pred_true$label,method_pred_true[[method]],1,type,cluster_num)
                                                                                                                                              })
                                                                                                                                              max_fscore_cluster <- unique_clusters[[which.max(cluster_f_betas)]]
                                                                                                                                              list(f1=max(unlist(cluster_f_betas)),acc=sum(type_pred_true[[method]]==max_fscore_cluster)/nrow(type_pred_true))
                                                                                                                                                  }))
                                                                        type_accuracy <- type_accuracy_f1$acc
                                                                        type_f1 <- type_accuracy_f1$f1
                                                                        type_unassigned_pctg <- NA
                                                                        supervised <- F
                                                                      }else if(method %in% supervised_methods){
                                                                        type_unassigned_pctg <- purrr::map_dbl(cell_types, function(type){ type_pred_true <- dplyr::select(combined_results,method,label) %>% dplyr::filter(label==type)
                                                                                                    analysis.assign.unlabeled_pctg(type_pred_true$label,type_pred_true[[method]])
                                                                                                    })
                                                                        type_accuracy <- purrr::map_dbl(cell_types, function(type){ type_pred_true <- dplyr::select(combined_results,method,label) %>% dplyr::filter(label==type)
                                                                                                    analysis.assign.accuracy(type_pred_true$label,type_pred_true[[method]])
                                                                                                    })
                                                                        type_f1 <- purrr::map_dbl(cell_types, function(type){ method_pred_true <- dplyr::select(combined_results,method,label) 
                                                                                                                              type_pred_true <- dplyr::filter(method_pred_true,label==type)
                                                                                                    analysis.cluster.fbeta(method_pred_true$label,method_pred_true[[method]],1,type,type)
                                                                                                    })
                                                                          
                                                                        supervised <- T
                                                                      }else{
                                                                        stop(str_glue("unkown method={method}"))
                                                                      }
                                                                      tibble(method=method,type=type_pctg$label,type_pctg=type_pctg$type_pctg,unlabeled_pctg=type_unassigned_pctg, type_accuracy=type_accuracy,type_f1=type_f1, supervised=supervised)
                                                                          }))
  

  experiments.analysis.attach_dataset_props("imbalance_impacts",exp_config,report_results=report_results,
                                            combined_results=combined_results,
                                            single_method_result=single_method_result,...)
}



experiments.analysis.unknown_types <- function(assign_results,cluster_results,exp_config,...){
  require(tidyverse)
  
  report_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  unknown_type <- exp_config$unknown_type
  methods <- colnames(combined_results)[colnames(combined_results)!="label"]
  supervised_methods <- colnames(assign_results)[colnames(assign_results)!="label"]
  unsupervised_methods <- colnames(cluster_results)[colnames(cluster_results)!="label"]
  single_method_result <- bind_rows(purrr::map(methods, function(method){
    if(method %in% unsupervised_methods){
      type_accuracy_f1 <- bind_rows(purrr::map(unknown_type, function(type){ method_pred_true <- dplyr::select(combined_results,method,label) 
      type_pred_true <- dplyr::filter(method_pred_true,label==type)
      unique_clusters <- unique(type_pred_true[[method]])                                                
      cluster_f_betas <- purrr::map(unique_clusters,function(cluster_num){ 
        analysis.cluster.fbeta(method_pred_true$label,method_pred_true[[method]],1,type,cluster_num)
      })
      max_fscore_cluster <- unique_clusters[[which.max(cluster_f_betas)]]
      list(f1=max(unlist(cluster_f_betas)),acc=sum(type_pred_true[[method]]==max_fscore_cluster)/nrow(type_pred_true))
      }))
      type_accuracy <- mean(type_accuracy_f1$acc)
      type_f1 <- mean(type_accuracy_f1$f1)
      type_unassigned_pctg <- NA
      supervised <- F
    }else if(method %in% supervised_methods){
      type_unassigned_pctg <- mean(purrr::map_dbl(unknown_type, function(type){ type_pred_true <- dplyr::select(combined_results,method,label) %>% dplyr::filter(label==type)
                                                                                 analysis.assign.unlabeled_pctg(type_pred_true$label,type_pred_true[[method]])
                                                                               }
                                                  )
                                   )
      
      type_accuracy <- mean(purrr::map_dbl(unknown_type, function(type){ type_pred <- dplyr::select(combined_results,method,label) %>% dplyr::filter(label==type) 
                                                                         sum(type_pred[,1]=="unassigned")/nrow(type_pred)
                                                                        }
                                           )
                            )
      type_f1 <- NA
      supervised <- T
    }else{
      stop(str_glue("unkown method={method}"))
    }
    tibble(method=method,unknown_celltype=unknown_type,unlabeled_pctg=type_unassigned_pctg, type_accuracy=type_accuracy,type_f1=type_f1, supervised=supervised)
  }))
  experiments.analysis.attach_dataset_props("celltype_detection",exp_config,report_results=report_results,
                                            combined_results=combined_results,
                                            single_method_result=single_method_result,...)
}


######calculate dataset properties and attach to results
####attached train and test dataset properties to experiment results
experiments.analysis.attach_dataset_props <- function(experiment,exp_config,...){
  extra_args <- list(...)
  train_dataset <- extra_args$train_dataset
  test_dataset <- extra_args$test_dataset
  total_dataset <- extra_args$total_dataset
  have_test <- exp_config$have_test
  report_results <- extra_args$report_results
  combined_results <- extra_args$combined_results
  if(experiment %in% c("imbalance_impacts","unknown_types")) single_method_result <-  extra_args$single_method_result
  if(have_test){
    if(experiment %in% c("batch_effects")) batch_effects_quant = analysis.dataset.batch_effects(train_dataset, test_dataset)
    train_data_props <- analysis.dataset.properties(train_dataset)
    test_data_props <- analysis.dataset.properties(test_dataset)
    train_test_corr <- analysis.interdata.correlation(train_dataset,test_dataset)
    total_data_props <- analysis.dataset.properties(total_dataset) %>% 
      purrr::list_modify("dataset"=NULL)
    for(prop in names(train_data_props)){
      report_results <- purrr::map(report_results,~{.[[str_glue("train_{prop}")]]=train_data_props[[prop]]
      return(.)
      })
    }
    for(prop in names(test_data_props)){
      report_results <- purrr::map(report_results,~{.[[str_glue("test_{prop}")]]=test_data_props[[prop]]
      return(.)
      })
    }
    for(prop in names(total_data_props)){
      report_results <- purrr::map(report_results,~{.[[str_glue("total_{prop}")]]=total_data_props[[prop]]
      return(.)
      })
    }
    
    report_results <- purrr::map(report_results,~{.$train_test_correlation=train_test_corr
    return(.)
    })
    if(experiment %in% c("batch_effects")){
      report_results <- purrr::map(report_results,~{.$batch_effects_amount=batch_effects_quant
      return(.)
      }) 
    }
    combined_results[["train_dataset"]] <- train_data_props[["dataset"]]
    combined_results[["test_dataset"]] <- test_data_props[["dataset"]]
    if(experiment %in% c("imbalance_impacts","unknown_types")){
      single_method_result[["train_dataset"]] <- train_data_props[["dataset"]]
      single_method_result[["test_dataset"]] <- test_data_props[["dataset"]]
    }
  }else{
    data_props <- analysis.dataset.properties(train_dataset)
    for(prop in names(data_props)){
      report_results <- purrr::map(report_results,~{.[[prop]]=data_props[[prop]]
      return(.)
      })
    }
    combined_results[["dataset"]] <- data_props[["dataset"]]
    if(experiment %in% c("imbalance_impacts","unknown_types")){
      single_method_result[["dataset"]] <- data_props[["dataset"]]
    }
  }
  if(experiment %in% c("imbalance_impacts","unknown_types")){
    return(list(pred_results=combined_results,analy_results=report_results,single_method_result=single_method_result))
  }else{
    return(list(pred_results=combined_results,analy_results=report_results))
  }
}  
  
  
