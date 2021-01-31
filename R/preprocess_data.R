#####convert datasets into singleCellExperiment R objects for later use

source("R/config.R")
source("R/utils.R")
source("R/dataset_config.R")
preprocess_dataset <- function(dataset=c("PBMC","pancreas","midbrain","ADASD","lung_cancer")) {
  switch(dataset,
         PBMC = preprocess_PBMC(dataset),
         pancreas = preprocess_pancreas(dataset),
         midbrain = preprocess_midbrain(dataset),
         ADASD = preprocess_ADASD(dataset),
         lung_cancer = preprocess_lung_cancer(dataset),
         stop("Unknown datasets")
  )  
}

##convert PBMC dataset into singcellexperiment R objects and save
preprocess_PBMC <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
  ###load PBMC_AllCells_withLabels.RData
  load(dataset_paths[[1]])
  #data <- eval(as.name(x)) ##rename PBMC dataset to variable data
  data <- c()
  labels <- vector("list", length(AllCounts))
  names(labels) <- names(AllCounts)
  for(nm in names(AllCounts)){
    data <- cbind(data, as(AllCounts[[nm]],'sparseMatrix'))
    labels[[nm]] <- rep(nm,dim(AllCounts[[nm]])[2])
  }
  labels <- flatten_chr(labels)
  rm(AllCounts)
  gc(TRUE)
  data <- utils.convert_to_SingleCellExperiment(data,rownames(data),colnames(data),tibble(label=labels),list(study='PBMC', study_name="All_PBMC"))
  rowData(data)$EnsembleId <- purrr::map_chr(rownames(data),~str_split(.,"\\t")[[1]][[1]])
  rowData(data)$geneName <- purrr::map_chr(rownames(data), ~str_split(.,"\\t")[[1]][[2]])
  new_dataset_path <- utils.get_dataset_paths(data_home,datasets[[dataset]][[1]])
  write_rds(data,new_dataset_path)
  
  ##load GSE96583_batch1_3_samples.RData
  dataset_names <- load(dataset_paths[[2]])
  sces <- map(map(dataset_names,~eval(as.name(.))), ~{colData(.)$sampleId <- colData(.)$ind
                                                    colData(.)$label <- colData(.)$cell.type
                                                    return(.)})
  sces <- utils.combine_SCEdatasets(sces,if_combined=T)
  metadata(sces) <- list(study='PBMC', study_name="GSE96583_batch1_3_samples")
  new_dataset_path <- utils.get_dataset_paths(data_home,datasets[[dataset]][[2]])
  write_rds(sces,new_dataset_path)
  
  ###load GSE96583_8_Stim_Pats.RData
  dataset_names <- load(dataset_paths[[3]])
  sces <- map(eval(as.name(dataset_names[[1]])), ~{colData(.)$sampleId <- colData(.)$ind
                                                   colData(.)$label <- colData(.)$cell
                                                   return(.)})
  sces <- utils.combine_SCEdatasets(sces,if_combined=T)
  metadata(sces) <- list(study='PBMC', study_name="GSE96583_8_Stim_Pats")
  new_dataset_path <- utils.get_dataset_paths(data_home,datasets[[dataset]][[3]])
  write_rds(sces,new_dataset_path)
  
  ###load GSE96583_8_Ctrl_Pats.RData
  dataset_names <- load(dataset_paths[[4]])
  sces <- map(eval(as.name(dataset_names[[1]])), ~{colData(.)$sampleId <- colData(.)$ind
                                                   colData(.)$label <- colData(.)$cell
                                                   return(.)})
  sces <- utils.combine_SCEdatasets(sces,if_combined=T)
  metadata(sces) <- list(study='PBMC', study_name="GSE96583_8_Ctrl_Pats")
  new_dataset_path <- utils.get_dataset_paths(data_home,datasets[[dataset]][[4]])
  write_rds(sces,new_dataset_path)
}

##convert Pancreas dataset into singlecellexperiment R objects and save
preprocess_pancreas <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
  new_dataset_path <- utils.get_dataset_paths(data_home,datasets[[dataset]])
  for(i in seq_along(dataset_paths)){
    study_name <- unlist(str_split(raw_datasets[[dataset]][[i]],"\\."))[[1]]
    print(str_glue('Preprocessing {study_name} dataset'))
    x <- load(dataset_paths[[i]])
    cnt <- eval(as.name(x[[1]])) %>%
      as('sparseMatrix')
    sampleId <- eval(as.name(x[[2]]))
    labels <- eval(as.name(x[[3]]))
    sampleId <- if(length(sampleId)!=length(labels)) rep(NA,length(labels)) else sampleId
    genes <- purrr::map_chr(rownames(cnt),~str_split(.,"__")[[1]][1])

    data <- utils.convert_to_SingleCellExperiment(cnt,genes,colnames(cnt),tibble(label=labels,sampleId=sampleId),
                                            list(study='pancreas',study_name=study_name))
    rowData(data)$geneName <- rownames(data)
    write_rds(data,new_dataset_path[[i]])
  }
}

###preprocess AD and autism datasets into singlecellexperiment object and save
preprocess_ADASD <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
  x <- load(dataset_paths[[1]])
  AD <- x[1:2]
  AD_cnt <- eval(as.name(AD[[1]])) %>%
    as('sparseMatrix')
  rm(list=AD[[1]])
  gc(T)
  AD_labels <- eval(as.name(AD[[2]]))
  AD_data <- utils.convert_to_SingleCellExperiment(AD_cnt,rownames(AD_cnt),colnames(AD_cnt),tibble(label=AD_labels),
                                                   list(study='ADASD',study_name="AD"))
  rowData(AD_data)$geneName <- rownames(AD_data)
  write_rds(AD_data,paste(data_home,"ADASD_AD.RDS"))
  
  autism <- x[3:4]
  autism_cnt <- eval(as.name(autism[[1]])) %>%
    as('sparseMatrix')
  rm(list=autism[[1]])
  gc(T)
  autism_labels <- eval(as.name(autism[[2]]))
  autism_data <- utils.convert_to_SingleCellExperiment(autism_cnt,rownames(autism_cnt),colnames(autism_cnt),tibble(label=autism_labels),
                                                       list(study='ADASD',study_name="autism"))
  rowData(autism_data)$geneName <- rownames(autism_data)
  write_rds(autism_data,paste(data_home,"ADASD_autism.RDS"))
}


##convert PBMC dataset into singcellexperiment R objects and save
preprocess_midbrain <- function(dataset){
  
}


##remove batch effects of datasets and saves
preprocess.remove_batch_effects <- function(sces,file_names){
  sces <- utils.remove_batch_effects(sces) 
  for(i in seq_along(sces)){
    path <- utils.get_dataset_paths(data_home,file_names[[i]])
    write_rds(sces[[i]],path)
  }
}

###intersect singlecellexperiments objects with different genes and save
preprocess.intersect_sces <- function(sces,file_names){
  paths <- utils.get_dataset_paths(data_home,file_names)
  if(sum(unlist(purrr::map(paths,file.exists)))>=1) {
    print("intersected datasets already exist, skipping intersection")  
    return()
  }
  sces <- utils.combine_SCEdatasets(sces,if_combined=FALSE)
  purrr::walk2(sces,paths,write_rds)
}