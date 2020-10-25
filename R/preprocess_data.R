#####convert datasets into singleCellExperiment R objects for later use

source("R/config.R")
source("R/utils.R")
source("R/dataset_config.R")
preprocess_dataset <- function(dataset=c("PBMC","pancreas")) {
  switch(dataset,
         PBMC = preprocess_PBMC(dataset),
         pancreas = preprocess_pancreas(dataset),
         stop("Unknown datasets")
  )  
}

##convert PBMC dataset into singcellexperiment R objects and save
preprocess_PBMC <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
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
  data <- utils.convert_to_SingleCellExperiment(data,rownames(data),colnames(data),tibble(label=labels),list(study='PBMC'))
  rowData(data)$EnsembleId <- map_chr(rownames(data),~str_split(.,"\\t")[[1]][[1]])
  rowData(data)$geneName <- map_chr(rownames(data), ~str_split(.,"\\t")[[1]][[2]])
  new_dataset_path <- utils.get_dataset_paths(data_home,datasets[[dataset]])
  write_rds(data,new_dataset_path)
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
    genes <- map_chr(rownames(cnt),~str_split(.,"__")[[1]][1])

    data <- utils.convert_to_SingleCellExperiment(cnt,genes,colnames(cnt),tibble(label=labels,sampleId=sampleId),
                                            list(study='pancreas',study_name=study_name))
    write_rds(data,new_dataset_path[[i]])
  }
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