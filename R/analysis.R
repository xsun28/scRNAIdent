source('R/utils.R')
source('R/experiments_config.R')
###separate assigned and unassigned and analyze respectively

analysis.run <- function(results,methods,metrics){
  metrics_functions_mapping <- c(ARI="analysis.cluster.ARI",
                                AMI="analysis.cluster.AMI",
                                FMI="analysis.cluster.FMI",
                                wNMI="analysis.cluster.wNMI",
                                wRI="analysis.cluster.wRI",
                                unlabeled_pctg="analysis.assign.unlabeled_pctg")
  
  analysis_results <- as.data.frame(matrix(nrow = length(methods), ncol = length(metrics)))
  rownames(analysis_results) <- methods
  colnames(analysis_results) <- metrics
  if(dim(results)[1] == 0) return(analysis_results)
  for (m in methods){
    for(metric in metrics){
      f <- get(metrics_functions_mapping[[metric]])
      analysis_results[m,metric] <- f(results$label,results[[m]])
    }
  }
  analysis_results
}

###Adjusted rand index
analysis.cluster.ARI <- function(true,pred){
  require(aricode)
  pred[is.na(pred)] <- -1 ####replace na with -1 unclustered cells
  ARI(true,pred)
}

###Adjusted mutual information
analysis.cluster.AMI <- function(true,pred){
  require(aricode)
  pred[is.na(pred)] <- -1
  AMI(true,pred)
}

###Fowlkes Mallows index
analysis.cluster.FMI <- function(true,pred){
  require(dendextend)
  pred[is.na(pred)] <- -1
  FM_index(true,pred)[[1]]
}

#####calculate the percentage of unlabeled cells for assign methods
analysis.assign.unlabeled_pctg <- function(labels,pred){
  unique_labels <- unique(labels)
  return(sum(unlist(purrr::map(pred, ~{!. %in% unique_labels})))/length(pred))
}


###### create cell hierarchy for wNMI
analysis.cluster.createHierarchy <- function(data,labels){
  require(Wind)
  stopifnot(is(data,"SingleCellExperiment"))
  ge_matrix <- as.matrix(counts(data))
  createRef(ge_matrix,labels)
}


####### create weights for wRI
analysis.culster.createWeights <- function(data,labels){
  require(Wind)
  stopifnot(is(data,"SingleCellExperiment"))
  ge_matrix <- as.matrix(counts(data))
  createWeights(ge_matrix,labels)
}
  
#######calculate wNMI
analysis.cluster.wNMI <- function(data,labels,preds){
  analysis.cluster.createHierarchy(data,labels) %>%
    wNMI(labels,preds)
}

#######calculate wRI
analysis.cluster.wRI <- function(data,labels,preds){
  weights <- analysis.culster.createWeights(data,labels)
  wRI(labels,preds,weights$W0,weights$W1)[[1]]
}
