source('R/config.R')

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

###select different sequencing depth cells
utils.seqDepthSelector <- function(data, depth_quantile){
  
}