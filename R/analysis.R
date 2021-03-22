
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
      analysis_results[m,metric] <- f(as.character(results$label),results[[m]])
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
  if(is.null(pred)==length(pred)||is.na(pred)==length(pred)) return(NA)
  return(sum(unlist(purrr::map(pred, ~{!. %in% unique_labels})))/length(pred))
}

#######calculate wNMI
analysis.cluster.wNMI <- function(labels,preds){
  require(Wind)
  wNMI(current_celltype_hierarchy,labels,preds)
}

#######calculate wRI
analysis.cluster.wRI <- function(labels,preds){
  require(Wind)
  wRI(labels,preds,current_celltype_weights$W0,current_celltype_weights$W1)[[1]]
}

analysis.pivot_table <- function(data,row,col1,col2=NULL,func="n()"){
  library(pivottabler)
  pt <- PivotTable$new()
  pt$addData(data)
  pt$addColumnDataGroups(col1)
  if(!is.null(col2)) pt$addColumnDataGroups(col2, expandExistingTotals=TRUE) 
  pt$addRowDataGroups(row)
  pt$defineCalculation(calculationName="TotalCounts", summariseExpression=func)
  pt$renderPivot()
}

analysis.dataset.complexity <- function(data){
  stopifnot(is(data,"SingleCellExperiment"))
  stopifnot("label" %in% colnames(colData(data)))
  agg_sce <- aggregateAcrossCells(data,ids = data$label)
  agg_count <- assay(agg_sce,"counts")
  corr <- cor(agg_count)
  diag(corr) <- 0
  cor <- apply(corr,1,max)
  complexity <- mean(cor)
  complexity
}

analysis.dataset.entropy <- function(data){
  stopifnot(is(data,"SingleCellExperiment"))
  stopifnot("label" %in% colnames(colData(data)))
  entropy <- sum((dplyr::group_by(as.data.frame(colData(data)),label)%>%dplyr::summarize(prob=n()/length(colData(data)$label)))[['prob']]%>%
            map_dbl(~{-log(.)*.}))
  entropy
}

analysis.dataset.sequencing_depth <- function(data){
  require(scater)
  stopifnot(is(data,"SingleCellExperiment"))
  seq_dep_stats <- as.list(summary(perCellQCMetrics(data)$sum/dim(data)[1]))
  c(median=seq_dep_stats$Median,IQR=(seq_dep_stats$`3rd Qu.`-seq_dep_stats$`1st Qu.`))
}

analysis.dataset.cell_types <- function(data,threshold=0){
  require(scater)
  stopifnot(is(data,"SingleCellExperiment"))
  dim(dplyr::group_by(as.data.frame(colData(data)),label)%>%dplyr::summarize(num=n())%>%dplyr::filter(num>=threshold))[1]
}

analysis.dataset.sparsity <- function(data){
  require(scater)
  stopifnot(is(data,"SingleCellExperiment"))
  mean(1-nexprs(data)/nrow(data))#pctg of undetected genes per cell
}

analysis.dataset.sample_num <- function(data){
  stopifnot(is(data,"SingleCellExperiment"))
  length(unique(colnames(data)))
}

analysis.dataset.gene_num <- function(data){
  stopifnot(is(data,"SingleCellExperiment"))
  length(unique(rownames(data)))
}

analysis.dataset.properties <- function(data){
  stopifnot(is(data,"SingleCellExperiment"))
  complexity <- analysis.dataset.complexity(data)
  entropy <- analysis.dataset.entropy(data)
  sequencing_depth <- analysis.dataset.sequencing_depth(data)
  cell_types <- analysis.dataset.cell_types(data)
  sparsity <- analysis.dataset.sparsity(data)
  sample_num <- analysis.dataset.sample_num(data)
  gene_num <- analysis.dataset.gene_num(data)
  list(name=metadata(data)$study_name,complexity=complexity,entropy=entropy,seq_depth_med=sequencing_depth[['median']],
      seq_depth_IQR=sequencing_depth[['IQR']],cell_types=cell_types,sparsity=sparsity,sample_num=sample_num,
      gene_num=gene_num)
}

analysis.dataset.properties_table <- function(dataset_names){
  prop_table <- vector("list",length(dataset_names))
  for(i in seq_along(dataset_names)){
    dataset_ <- utils.load_datasets(dataset_names[[i]])
    dataset_prop <- analysis.dataset.properties(dataset_)
    prop_table[[i]] <- dataset_prop
  }
  bind_rows(prop_table)
}

analysis.dataset.batch_effects <- function(data1,data2){
  stopifnot(is(data1,"SingleCellExperiment"))
  stopifnot(is(data1,"SingleCellExperiment"))
  cell_types1 <- unique(colData(data1)$label)
  cell_types2 <- unique(colData(data2)$label)
  common_types <- intersect(cell_types1,cell_types2)
  
}
