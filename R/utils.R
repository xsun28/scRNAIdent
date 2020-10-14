source('R/config.R')
library(purrr)
###load a dataset
utils.load_dataset <- function(file_name){
  data_path <- paste(data_home,file_name,sep='')
  return(read_rds(data_path))    
}

###convert count matrix to singlecellexperiment object and add cell-,gene-level metadata
utils.convert_to_SingleCellExperiment <- function(data_matrix,genes,cellIds,coldata,meta){
  require(SingleCellExperiment)
  require(scater)
  colnames(data_matrix) <- cellIds
  rownames(data_matrix) <- genes
  sce <- SingleCellExperiment(list(counts=data_matrix),colData=coldata,metadata=meta)
  sce <- addPerCellQC(sce)
  rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE))
  sce
}

utils.get_dataset_paths <- function(data_home,dataset_names){
  require(tidyverse)
  paths <- map_chr(dataset_names,~paste(data_home,.,sep=''))
}

###filter out cells with all 0s reads and genes with all 0s reads
utils.filter <- function(data,filter_gene=TRUE, filter_cells=TRUE){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  if(filter_gene & filter_cells){
    return(data[rowData(data)$count>0,colData(data)$detected>0])
  }else if (filter_cell){
    return(data[,colData(data)$detected>0])
  }else if(filter_gene){
    return(data[rowData(data)$count>0,])
  }else{
    return(data)
  }
}

###sampling sample_num cells from each cell type
utils.sampler <- function(data,sample_num,types){
  sample_idx <- map(types,~ which(colData(data)$label==.)) %>% 
    map(~sample(.,sample_num,replace=FALSE)) 
  data[,unlist(sample_idx)]
}

###select different sequencing depth cells, right means take the deeper sequencing data
utils.seqDepthSelector <- function(data, quantile,right=TRUE){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  if(right){
    return(data[,which(percent_rank(colData(data)$sum)>=quantile)])
  }else{
    return(data[,which(percent_rank(colData(data)$sum)<=quantile)])
  }
}

###label unassigned cells in the tibble results
utils.label_unassigned <- function(assign_results){
  cell_types <- unique(assign_results$label)
  f <- function(x){
    x[which(!x%in% cell_types)] <- 'unassigned'
    x
  }
  map(assign_results,~f(.))
}

###select assigned cells by all methods from tibble results
utils.select_assigned <- function(assign_results){
    if("label" %in% names(assign_results)){
      assign_results1 <- select(assign_results,-label)
    }
    labeled_idx <- reduce(map(assign_results1,~which(.!='unassigned')),intersect)
    assign_results[labeled_idx,]
}

###select unassigned cells by all methods from tibble results
utils.select_unassigned <- function(assign_results){
  if("label" %in% names(assign_results)){
    assign_results1 <- select(assign_results,-label)
  }
  labeled_idx <- reduce(map(assign_results1,~which(.!='unassigned')),intersect)
  assign_results[-labeled_idx,]
}


