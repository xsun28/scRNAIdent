source('R/marker_gene.R')
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
utils.combine_SCEdatasets <- function(sces,if_combined=TRUE,colDatas_cols=NULL){
    require(purrr)
    require(SingleCellExperiment)
    require(scater)
    intersected_genes <- purrr::reduce(purrr::map(sces,~rownames(rowData(.))),intersect)
    if(if_combined){
      if(purrr::is_null(colDatas_cols)){
        colDatas_cols <- if('batch' %in% colnames(colData(sces[[1]]))) c('batch','label','sampleId') else c('label','sampleId')
        }
      colDatas <- purrr::reduce((purrr::map(sces,~colData(.)[,colDatas_cols])),rbind)
     
      metaDatas <- purrr::reduce((purrr::map(sces,~metadata(.))),bind_rows)
      allcounts <- purrr::reduce(purrr::map(sces,~as.matrix(counts(.)[intersected_genes,])),function(x,y){
                                                                              utils.col2rowNames(merge(x,y,by=0,all=FALSE),1)
                                                                            })
      sce <- SingleCellExperiment(list(counts=as.matrix(allcounts)),
                                colData=colDatas,
                                metadata=metaDatas)
      sce <- addPerCellQC(sce)
      rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE),geneName=intersected_genes)
      return(sce)
    }else{
      return(purrr::map(sces,~.[intersected_genes,]))
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
  rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE),geneName=genes)
  sce
}

utils.get_dataset_paths <- function(data_home,dataset_names){
  require(tidyverse)
  paths <- purrr::map_chr(dataset_names,~paste(data_home,.,sep=''))
  paths
}

###filter out cells with all 0s reads and genes with all 0s reads
utils.filter <- function(data,filter_gene=TRUE, filter_cells=TRUE, filter_cell_type=TRUE){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  study_name <- metadata(data)$study_name[[1]]
  properties <- dataset.properties[[study_name]]
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
utils.sampler <- function(data,sample_num,types,column="label"){
  sample_idx <- purrr::map(types,~ which(colData(data)[[column]]==.)) %>% 
    purrr::map(~sample(.,min(sample_num,length(.)),replace=FALSE)) 
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
  purrr::map(assign_results,~f(.))
}

###select assigned cells by all methods from tibble results

utils.select_assigned <- function(results){
  assign_results <- results$assign_results
  cluster_results <- results$cluster_results
  labels <- unique(assign_results$label)
  if("label" %in% names(assign_results)){
    assign_results <- dplyr::select(assign_results,-label)
    cluster_results <- dplyr::select(cluster_results,-label)
  }
  labeled_idx_assign <- purrr::reduce(purrr::map(assign_results,~which(. %in% labels)),intersect)
  labeled_idx_cluster <- purrr::reduce(purrr::map(cluster_results,~which(!is.na(.))),intersect)
  labeled_idx <- intersect(labeled_idx_assign, labeled_idx_cluster)
  purrr::map(results,~.[labeled_idx,])
}

###select unassigned cells by all methods from tibble results
utils.select_unassigned <- function(results){
  assign_results <- results$assign_results
  cluster_results <- results$cluster_results
  labels <- unique(assign_results$label)
  if("label" %in% names(assign_results)){
    assign_results <- dplyr::select(assign_results,-label)
    cluster_results <- dplyr::select(cluster_results,-label)
  }

  labeled_idx_assign <- purrr::reduce(purrr::map(assign_results,~which(. %in% labels)),intersect)
  labeled_idx_cluster <- purrr::reduce(purrr::map(cluster_results,~which(!is.na(.))),intersect)
  labeled_idx <- intersect(labeled_idx_assign, labeled_idx_cluster)
  if(length(labeled_idx)==0){
    return(results)
  }
  purrr::map(results,~.[-labeled_idx,])
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
  batches <- do.call(multiBatchNorm,as.list(batches))
  decs <- purrr::map(batches,modelGeneVar)
  combined.decs <- do.call('combineVar',decs)
  chosen.hvgs <- combined.decs$bio > 0
  args <- list()
  args$subset.row <- chosen.hvgs
  f.out <- do.call(fastMNN,c(as.list(batches),args))
  sces_batch_free <- purrr::map(batches,~f.out[,colnames(.)])
  add_col_row_data <- function(x,y){
    colData(x)[,c('label','sampleId')] <- colData(y)[,c('label','sampleId')]
    rowData(x)[,'count'] <- rowData(y)[,'count']
    rowData(x)[,'geneName'] <- rowData(y)[,'geneName']
    counts(x) <- exp(assays(x)$reconstructed)
    metadata(x) <- metadata(y)
    x <- addPerCellQC(x)
  }
  purrr::map2(sces_batch_free,batches,add_col_row_data)
}

####append one singlecellexperiment onto another
utils.append_sce <- function(sce1,sce2){
  common_cols <- intersect(colnames(colData(sce1)),colnames(colData(sce2)))
  colDatas <- rbind(colData(sce1)[common_cols],colData(sce2)[common_cols])
  metaDatas <- bind_rows(metadata(sce1),metadata(sce2))
  allcounts <- utils.col2rowNames(merge(as.matrix(counts(sce1)),as.matrix(counts(sce2)),by=0,all=F),1)  

  allcounts <- SingleCellExperiment(list(counts=as.matrix(allcounts)),
                              colData=colDatas,
                              metadata=metaDatas)
  rowData(allcounts) <- tibble(count=nexprs(allcounts,byrow=TRUE),geneName=rownames(allcounts))
  return(allcounts)
}

utils.get_methods <- function(data){
  stopifnot(is(data,"data.frame"))
  data$methods <- unlist(purrr::map(rownames(data),~str_split(.,"\\.\\.\\.")[[1]][[1]]))
  data
}

utils.check_marker_genes <- function(data,marker_gene_file,generated_marker_gene_file,marker_gene_method){
  
  if(purrr::is_null(marker_gene_file)){
    if(!file.exists(str_glue("{marker_home}/{generated_marker_gene_file}.RDS"))){
      print(str_glue("Generating {marker_home}/{generated_marker_gene_file}.RDS"))
      markers_mat <- markergene.generate_marker_genes(marker_gene_method,data,str_glue("{marker_home}/{generated_marker_gene_file}"))
    }
    else{
      print(str_glue("Loading {marker_home}/{generated_marker_gene_file}.RDS"))
      markers_mat <- read_rds(str_glue("{marker_home}/{generated_marker_gene_file}.RDS"))
    }
    matchidx <- match(rownames(markers_mat), rownames(data))
    markers_mat <- markers_mat[!is.na(matchidx),]
  }else{
    markers<- read.csv(str_glue("{marker_home}/{marker_gene_file}"))
    markers_mat <- matrix(0, nrow(markers), length(unique(markers$CellType)))
    for(i in 1:nrow(markers)) {
      idx0 <- match(markers$CellType[i], unique(markers$CellType))
      markers_mat[i, idx0] <- 1
    }
    rownames(markers_mat) <- markers$Marker
    colnames(markers_mat) <- unique(markers$CellType)
    gene_names <- if(length(rowData(data)$geneName)>0) rowData(data)$geneName else rownames(data) ##PBMC data gene name is not same as rownames
    matchidx <- match(markers[,1], gene_names)
    markers_mat <- markers_mat[!is.na(matchidx),]
  }
  list(markers_mat=markers_mat,matchidx=matchidx)
}

utils.update_batch_effects_free_config <- function(experiment){
  experiments.methods[[experiment]]$cluster <<- experiments.methods[[experiment]]$cluster_batch_free
  print(str_glue("batch effects free cluster methods are: {experiments.methods[[experiment]]$cluster}"))
  for(m in experiments.methods[[experiment]]$cluster){
    if(exists(str_glue("methods.config.{m}.batch_free"))){
      assign(str_glue("methods.config.{m}"), get(str_glue("methods.config.{m}.batch_free")),envir = .GlobalEnv)
    }
  }
  experiments.methods[[experiment]]$assign <<- experiments.methods[[experiment]]$assign_batch_free
  print(str_glue("batch effects free assign methods are: {experiments.methods[[experiment]]$assign}"))
  for(m in experiments.methods[[experiment]]$assign){
    if(exists(str_glue("methods.config.{m}.batch_free"))){
      assign(str_glue("methods.config.{m}"),get(str_glue("methods.config.{m}.batch_free")),envir = .GlobalEnv)
    }
  }
}

utils.try_catch_method_error <- function(code, silent=FALSE){
  tryCatch(code, error=function(c) {
    msg <- conditionMessage(c)
    if (!silent) warning(c)
    invisible(structure(msg, class = "try-error"))
  })
}


###### create cell hierarchy for wNMI
utils.createCellTypeHierarchy <- function(data,labels){
  require(Wind)
  stopifnot(is(data,"SingleCellExperiment"))
  data <- as.matrix(counts(data))
  createRef(data,labels)
}


####### create weights for wRI
utils.createCellTypeWeights <- function(data,labels){
  require(Wind)
  stopifnot(is(data,"SingleCellExperiment"))
  data <- as.matrix(counts(data))
  createWeights(data,labels)
}



#####convert mouse genes to human genes
utils.convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  names(genesV2) <- c("genes","human_gene")
  genesV2
}

#####convert from ensemblId to gene symbols 
utils.convert2GeneSymbols <- function(ensemblIDs){
  require("org.Hs.eg.db") 
  symbols <- mapIds(org.Hs.eg.db, keys = ensemblIDs, keytype = "ENSEMBL", column="SYMBOL")
  symbols
}

#####convert from gene symbols to ensemblId 
utils.convert2EnsemblIDs <- function(symbols){
  require("org.Hs.eg.db") 
  ensembls <- mapIds(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", column="ENSEMBL")
  ensembls
}

