source('R/config.R')
library(tidyverse)
library(stringr)
###load a dataset
utils.load_datasets <- function(file_names){
  data_paths <- paste(data_home,file_names,sep='')
  print(str_glue('data path is {data_paths}'))
  if(length(data_paths)==1){return(read_rds(data_paths))}
  sces <- list()
  for(i in seq_along(data_paths)){
    sces[i] <- read_rds(data_paths[[i]])
  }
  utils.combine_SCEdatasets(sces)
}


###combine multiple singlecellexperiment into one
utils.combine_SCEdatasets <- function(sces){
    colDatas <- purrr::reduce((map(sces,~colData(.))),rbind)
    intersected_genes <- purrr::reduce(map(sces,~rownames(rowData(.))),intersect)
    metaDatas <- purrr::reduce((map(sces,~metadata(.))),bind_rows)
    allcounts <- purrr::reduce(map(sces,~counts(.)[intersected_genes,]),function(x,y){
                                                                            utils.col2rowNames(merge(x,y,by=0,all=FALSE),1)
                                                                          })
    sce<-SingleCellExperiment(list(counts=as.matrix(allcounts)),
                              colData=colDatas,
                              metadata=metaDatas)
    sce <- addPerCellQC(sce)
    rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE))
    sce
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
utils.select_assigned <- function(results){
    assign_results <- results$assign_results
    if("label" %in% names(assign_results)){
      assign_results1 <- select(assign_results,-label)
    }
    labeled_idx <- purrr::reduce(map(assign_results1,~which(.!='unassigned')),intersect)
    map(results,~.[labeled_idx,])
}

###select unassigned cells by all methods from tibble results
utils.select_unassigned <- function(results){
  assign_results <- results$assign_results
  if("label" %in% names(assign_results)){
    assign_results1 <- select(assign_results,-label)
  }
  labeled_idx <- purrr::reduce(map(assign_results1,~which(.!='unassigned')),intersect)
  map(results,~.[-labeled_idx,])
}

###convert column into rownames
utils.col2rowNames <- function(data,col){
  stopifnot(is(data,"data.frame"))
  rownames(data) <- data[,col]
  data[,col] <- NULL
  data
}

###remove batch effects
utils.remove_batch_effects <- function(data1, data2){
  
}
