####Wrapper
experiments.analysis <- function(experiment, assign_results,cluster_results,exp_config){
  switch(experiment,
         simple_accuracy = experiments.analysis.simple_accuracy(assign_results,cluster_results,exp_config),
         cell_number = experiments.analysis.cell_number(assign_results,cluster_results,exp_config),
         celltype_number = experiments.analysis.celltype_number(assign_results,cluster_results,exp_config),
         sequencing_depth = experiments.analysis.sequencing_depth(assign_results,cluster_results,exp_config),
         celltype_structure = experiments.analysis.celltype_structure(assign_results,cluster_results,exp_config),
         batch_effects = experiments.analysis.batch_effects(assign_results,cluster_results,exp_config),
         inter_diseases = experiments.analysis.inter_diseases(assign_results,cluster_results,exp_config),
         celltype_complexity = experiments.analysis.celltype_complexity(assign_results,cluster_results,exp_config),
         inter_species = experiments.analysis.inter_species(assign_results,cluster_results,exp_config),
         random_noise = experiments.analysis.random_noise(assign_results,cluster_results,exp_config),
         inter_protocol = experiments.analysis.inter_protocol(assign_results,cluster_results,exp_config),
         imbalance_impacts = experiments.analysis.imbalance_impacts(assign_results,cluster_results,exp_config),
         stop("Wrong experiment")
  )
}


###base function 
experiments.analysis.base <- function(assign_results,cluster_results,exp_config){
  print("start analyzing assigned/unassigned assigning/clustering methods")
  if(missing(exp_config)){
    stop("missing exp_config in base analysis function")
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

experiments.analysis.simple_accuracy <- function(assign_results,cluster_results,exp_config){
  return(experiments.analysis.base(assign_results,cluster_results,exp_config))
}

experiments.analysis.cell_number <- function(assign_results,cluster_results,exp_config){
  return(experiments.analysis.base(assign_results,cluster_results,exp_config))
}

experiments.analysis.sequencing_depth <- function(assign_results,cluster_results,exp_config){
  return(experiments.analysis.base(assign_results,cluster_results,exp_config))
}

experiments.analysis.batch_effects <- function(assign_results,cluster_results,exp_config){
  return(experiments.analysis.base(assign_results,cluster_results,exp_config))
}

experiments.analysis.inter_diseases <- function(assign_results,cluster_results,exp_config){
  return(experiments.analysis.base(assign_results,cluster_results,exp_config))
}

experiments.analysis.imbalance_impacts <- function(assign_results,cluster_results,exp_config){
  require(tidyverse)

  
  
  base_results <- experiments.analysis.base(assign_results,cluster_results,exp_config)
  combined_results <- bind_cols(assign_results,dplyr::select(cluster_results,-label))
  type_pctg <- dplyr::group_by(combined_results, label) %>% dplyr::summarize(type_pctg=round(n()/length(combined_results$label),2))
  cell_types <- type_pctg$label
  methods <- colnames(combined_results)[colnames(combined_results)!="label"]
  supervised_methods <- colnames(assign_results)[colnames(assign_results)!="label"]
  unsupervised_methods <- colnames(cluster_results)[colnames(cluster_results)!="label"]
  single_method_result <- bind_rows(purrr::map(methods, function(method){
                                                                      if(method %in% unsupervised_methods){
                                                                        type_accuracy_f1 <- bind_rows(purrr::map(cell_types, function(type){ method_pred_true <- dplyr::select(combined_results,method,label) 
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
  base_results$single_method_result <- single_method_result
  base_results
}
