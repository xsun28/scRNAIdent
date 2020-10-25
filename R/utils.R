source('R/config.R')
source('R/dataset_config.R')
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
  sces
}


###combine multiple singlecellexperiment into one
utils.combine_SCEdatasets <- function(sces,if_combined=TRUE){
    require(purrr)
    require(SingleCellExperiment)
    require(scater)
    intersected_genes <- purrr::reduce(map(sces,~rownames(rowData(.))),intersect)
    if(if_combined){
      colDatas_cols <- if('batch' %in% colnames(colData(sces[[1]]))) c('batch','label','sampleId') else c('label','sampleId')
      colDatas <- purrr::reduce((map(sces,~colData(.)[,colDatas_cols])),rbind)
     
      metaDatas <- purrr::reduce((map(sces,~metadata(.))),bind_rows)
      allcounts <- purrr::reduce(map(sces,~as.matrix(counts(.)[intersected_genes,])),function(x,y){
                                                                              utils.col2rowNames(merge(x,y,by=0,all=FALSE),1)
                                                                            })
      sce <- SingleCellExperiment(list(counts=as.matrix(allcounts)),
                                colData=colDatas,
                                metadata=metaDatas)
      sce <- addPerCellQC(sce)
      rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE))
      return(sce)
    }else{
      return(map(sces,~.[intersected_genes,]))
    }
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
utils.filter <- function(data,filter_gene=TRUE, filter_cells=TRUE, filter_cell_type=TRUE){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  study <- metadata(data)$study[[1]]
  properties <- dataset.properties[[study]]
  threshold <- properties$sample_threshold
  if(!is_null(properties$cell_types)){
    print(str_glue("only keep cell types with marker genes: {properties$cell_types}"))
    data <- data[,colData(data)$label %in% properties$cell_types] ###if only keeps cell types with marker gene files
  }
  if(filter_cell_type){
    cell_type_num <- table(colData(data)$label)
    cell_types_kept <- names(cell_type_num[cell_type_num>threshold])
    print(str_glue("only keep cell types:{cell_types_kept} which have more cells than {threshold}"))
    data <- data[,colData(data)$label %in% cell_types_kept]
  }
  if(filter_gene & filter_cells){
    return(data[rowData(data)$count>0,colData(data)$detected>0])
  }else if (filter_cells){
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
  cluster_results <- results$cluster_results
  if("label" %in% names(assign_results)){
    assign_results <- dplyr::select(assign_results,-label)
    cluster_results <- dplyr::select(cluster_results,-label)
  }
  labeled_idx_assign <- purrr::reduce(map(assign_results,~which(.!='unassigned')),intersect)
  labeled_idx_cluster <- purrr::reduce(map(cluster_results,~which(!is.na(.))),intersect)
  labeled_idx <- intersect(labeled_idx_assign, labeled_idx_cluster)
  map(results,~.[labeled_idx,])
}

###select unassigned cells by all methods from tibble results
utils.select_unassigned <- function(results){
  assign_results <- results$assign_results
  cluster_results <- results$cluster_results
  if("label" %in% names(assign_results)){
    assign_results <- dplyr::select(assign_results,-label)
    cluster_results <- dplyr::select(cluster_results,-label)
  }
  labeled_idx_assign <- purrr::reduce(map(assign_results,~which(.!='unassigned')),intersect)
  labeled_idx_cluster <- purrr::reduce(map(cluster_results,~which(!is.na(.))),intersect)
  labeled_idx <- intersect(labeled_idx_assign, labeled_idx_cluster)
  map(results,~.[-labeled_idx,])
}

###convert column into rownames
utils.col2rowNames <- function(data,col){
  stopifnot(is(data,"data.frame"))
  rownames(data) <- data[,col]
  data[,col] <- NULL
  data
}

###remove batch effects using MNN
utils.remove_batch_effects <- function(batches){
  require(batchelor)
  require(scater)
  require(scran)
  batches <- multiBatchNorm(batches)
  decs <- map(batches,modelGeneVar)
  combined.decs <- do.call('combineVar',decs)
  chosen.hvgs <- combined.decs$bio > 0
  f.out <- fastMNN(batches,subset.row=chosen.hvgs)
  sces_batch_free <- map(batches,~f.out[,colnames(.)])
  add_col_row_data <- function(x,y){
    colData(x)[,c('label','sampleId')] <- colData(y)[,c('label','sampleId')]
    rowData(x)[,'count'] <- rowData(y)[,'count']
    counts(x) <- exp(assays(x)$reconstructed)
    metadata(x) <- metadata(y)
    x <- addPerCellQC(x)
  }
  map2(sces_batch_free,batches,add_col_row_data)
}

####append one singlecellexperiment onto another
utils.append_sce <- function(sce1,sce2){
  colDatas <- rbind(colData(sce1),colData(sce2))
  metaDatas <- bind_rows(metadata(sce1),metadata(sce2))
  allcounts <- utils.col2rowNames(merge(as.matrix(counts(sce1)),as.matrix(counts(sce2)),by=0),1)  

  sce <- SingleCellExperiment(list(counts=as.matrix(allcounts)),
                              colData=colDatas,
                              metadata=metaDatas)
  rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE))
  return(sce)
}