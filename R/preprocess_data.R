#####convert datasets into singleCellExperiment R objects for later use

source("R/config.R")
source("R/utils.R")
source("R/dataset_config.R")

preprocesss.datasets <- list(PBMC=list("PBMC_AllCells_withLabels.RDS","GSE96583_batch{i}_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS"),
                             pancreas=list(muraro="Muraro_pancreas_clean.RDS",
                                           seger="Segerstolpe_pancreas_clean.RDS",
                                           xin="Xin_pancreas_clean.RDS"
                             ),
                             hcl=list("human_cell_landscape.RDS"),
                             mca=list("mouse_cell_landscape.RDS"),
                             ADASD=list(AD="ADASD_AD.RDS",autism="ADASD_autism.RDS"),
                             midbrain=list(human="midbrain_human.RDS",mouse="midbrain_mouse.RDS"),
                             cellbench=list(tenx="cellbench_10x.RDS",CELseq2="cellbench_CELseq2.RDS",Dropseq="cellbench_Dropseq.RDS")
)



preprocess_dataset <- function(dataset=c("PBMC","pancreas","midbrain","ADASD","cellbench","hcl")) {
  switch(dataset,
         PBMC = preprocess_PBMC(dataset),
         pancreas = preprocess_pancreas(dataset),
         midbrain = preprocess_midbrain(dataset),
         ADASD = preprocess_ADASD(dataset),
         cellbench = preprocess_cellbench(dataset),
         hcl = preprocess_human_cell_landscape(dataset),
         mca = preprocess_mouse_cell_atlas(dataset),
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
  data <- utils.convert_to_SingleCellExperiment(data,rownames(data),colnames(data),tibble(label=labels),list(study='PBMC', study_name="PBMC_AllCells_withLabels"))
  rowData(data)$EnsembleId <- purrr::map_chr(rownames(data),~str_split(.,"\\t")[[1]][[1]])
  rowData(data)$geneName <- purrr::map_chr(rownames(data), ~str_split(.,"\\t")[[1]][[2]])
  rownames(data) <- rowData(data)$geneName
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]][[1]])
  write_rds(data,new_dataset_path)
  
  ##load GSE96583_batch1_3_samples.RData
  map_path <- str_glue("{type_home}/{dataset.properties$GSE96583_batch1_samples$cell_type_map}")
  type_map <- read_csv(map_path)
  dataset_names <- load(dataset_paths[[2]])
  sces <- purrr::map(purrr::map(dataset_names,~eval(as.name(.))), ~{
    colData(.)$sampleId <- rownames(colData(.))
    colData(.)$label <- utils.convertCellTypes(colData(.)$cell.type,type_map)
    counts(.) <- as(counts(.),'sparseMatrix')
    return(.)})
  # sces <- purrr::map(sces,addPerCellQC)
  colData_cols <- colnames(colData(sces[[1]]))
  # sces <- utils.combine_SCEdatasets(sces,if_combined=T,colData_cols)
  for(i in seq_along(sces)){
    sce <- sces[[i]]
    sce <- sce[,!is.na(sce$label)]
    sce <- addPerCellQC(sce)
    sce <- sce[,which(!quickPerCellQC(sce)$discard)]
    metadata(sce) <- list(study='PBMC', study_name=str_glue("GSE96583_batch{i}_samples"))
    rowData(sce)$count <- nexprs(sce,byrow=TRUE)
    rowData(sce)$geneName <- rownames(sce)
    rowData(sce)$EnsembleId <- utils.convert2EnsemblIDs(rownames(sce))
    colData(sce)$unique_id <- 1:ncol(sce)
    new_dataset_path <- utils.get_dataset_paths(data_home,str_glue(preprocesss.datasets[[dataset]][[2]]))
    write_rds(sce,new_dataset_path)
  }
  
  
  ###load GSE96583_8_Stim_Pats.RData
  map_path <- str_glue("{type_home}/{dataset.properties$GSE96583_8_Stim_Pats$cell_type_map}")
  type_map <- read_csv(map_path)
  dataset_names <- load(dataset_paths[[3]])
  sces <- purrr::map(eval(as.name(dataset_names[[1]])), ~{colData(.)$sampleId <- rownames(colData(.))
  colData(.)$label <- utils.convertCellTypes(colData(.)$cell,type_map)
  counts(.) <- as(counts(.),'sparseMatrix')
  return(.)})
  
  colData_cols <- colnames(colData(sces[[1]]))
  sces <- utils.combine_SCEdatasets(sces,if_combined=T,colData_cols)
  sces <- sces[,!is.na(sces$label)]
  sces <- sces[,which(!quickPerCellQC(sces)$discard)]
  metadata(sces) <- list(study='PBMC', study_name="GSE96583_8_Stim_Pats")
  rowData(sces)$geneName <- rownames(sces)
  rowData(sces)$EnsembleId <- utils.convert2EnsemblIDs(rownames(sces))
  colData(sces)$unique_id <- 1:ncol(sces)
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]][[3]])
  write_rds(sces,new_dataset_path)
  
  ###load GSE96583_8_Ctrl_Pats.RData
  map_path <- str_glue("{type_home}/{dataset.properties$GSE96583_8_Ctrl_Pats$cell_type_map}")
  type_map <- read_csv(map_path)
  dataset_names <- load(dataset_paths[[4]])
  sces <- purrr::map(eval(as.name(dataset_names[[1]])), ~{colData(.)$sampleId <- rownames(colData(.))
  colData(.)$label <- utils.convertCellTypes(colData(.)$cell,type_map)
  counts(.) <- as(counts(.),'sparseMatrix')
  return(.)})
  
  colData_cols <- colnames(colData(sces[[1]]))
  sces <- utils.combine_SCEdatasets(sces,if_combined=T,colData_cols)
  sces <- sces[,!is.na(sces$label)]
  sces <- sces[,which(!quickPerCellQC(sces)$discard)]
  metadata(sces) <- list(study='PBMC', study_name="GSE96583_8_Ctrl_Pats")
  rowData(sces)$geneName <- rownames(sces)
  rowData(sces)$EnsembleId <- utils.convert2EnsemblIDs(rownames(sces))
  colData(sces)$unique_id <- 1:ncol(sces)
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]][[4]])
  write_rds(sces,new_dataset_path)
}

##convert Pancreas dataset into singlecellexperiment R objects and save
preprocess_pancreas <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]])
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
    
    data <- utils.convert_to_SingleCellExperiment(cnt,genes,colnames(cnt),tibble(label=labels,sampleId=sampleId,ind=sampleId),
                                                  list(study='pancreas',study_name=study_name))
    rowData(data)$geneName <- rownames(data)
    rowData(data)$EnsembleId <- utils.convert2EnsemblIDs(rownames(data))
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
                                                   list(study='ADASD',study_name="ADASD_AD"))
  rowData(AD_data)$geneName <- rownames(AD_data)
  rowData(AD_data)$EnsembleId <- utils.convert2EnsemblIDs(rownames(AD_data))
  write_rds(AD_data,paste(data_home,"ADASD_AD.RDS",sep = ""))
  
  autism <- x[3:4]
  autism_cnt <- eval(as.name(autism[[1]])) %>%
    as('sparseMatrix')
  rm(list=autism[[1]])
  gc(T)
  autism_labels <- eval(as.name(autism[[2]]))
  autism_data <- utils.convert_to_SingleCellExperiment(autism_cnt,rownames(autism_cnt),colnames(autism_cnt),tibble(label=autism_labels),
                                                       list(study='ADASD',study_name="ADASD_autism"))
  rowData(autism_data)$EnsembleId <- purrr::map_chr(rownames(autism_data),~str_split(.,"\\t")[[1]][[1]])
  rowData(autism_data)$geneName <- purrr::map_chr(rownames(autism_data), ~str_split(.,"\\t")[[1]][[2]])
  rownames(autism_data) <- rowData(autism_data)$geneName
  write_rds(autism_data,paste(data_home,"ADASD_autism.RDS",sep = ""))
}


##convert PBMC dataset into singcellexperiment R objects and save
preprocess_midbrain <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
  dataset_names <- load(dataset_paths[[1]])
  sces <- purrr::map(purrr::map(dataset_names,~eval(as.name(.))), ~{colData(.)$sampleId <- rownames(colData(.))
  colData(.)$label <- colData(.)$cps
  colData(.)$species <- purrr::map_chr(colData(.)$celltypes,~str_sub(.,1,1))
  counts(.) <- as(counts(.),'sparseMatrix')
  return(.)})
  
  sces <- purrr::map(sces,addPerCellQC)
  converted_genes <- utils.convertMouseGeneList(rowData(sces[[1]])$genes)
  converted_genes <- merge(rowData(sces[[1]]),converted_genes,by="genes",all.x=T)
  converted_genes$EnsembleId <- utils.convert2EnsemblIDs(converted_genes$human_gene)
  converted_genes$geneName <- converted_genes$human_gene
  new_mouse_sces <- sces[[1]][converted_genes$genes,]
  new_mouse_sces <- new_mouse_sces[,which(!quickPerCellQC(new_mouse_sces)$discard)]
  rowData(new_mouse_sces) <- converted_genes
  rowData(new_mouse_sces)$count <- nexprs(new_mouse_sces,byrow=TRUE)
  colData(new_mouse_sces)$unique_id <- 1:ncol(new_mouse_sces)
  colData(new_mouse_sces)$label <- map_chr(new_mouse_sces$label, ~{new <- dataset.properties[["midbrain_mouse"]]$cell_type_map[[.]]
  if(is_null(new)) return(.) else return(new)})
  rownames(new_mouse_sces) <- converted_genes$human_gene
  new_mouse_sces <- new_mouse_sces[-which(purrr::map_lgl(rownames(new_mouse_sces),is.na)),]
  
  sces[[2]] <- sces[[2]][,which(!quickPerCellQC(sces[[2]])$discard)]
  rowData(sces[[2]])$human_gene <- rowData(sces[[2]])$genes
  rowData(sces[[2]])$EnsembleId <- utils.convert2EnsemblIDs(rowData(sces[[2]])$human_gene)
  rowData(sces[[2]])$geneName <- rowData(sces[[2]])$human_gene
  rowData(sces[[2]])$count <- nexprs(sces[[2]],byrow=TRUE)
  colData(sces[[2]])$unique_id <- 1:ncol(sces[[2]])
  colData(sces[[2]])$label <- map_chr(sces[[2]]$label, ~{new <- dataset.properties[["midbrain_human"]]$cell_type_map[[.]]
  if(is_null(new)) return(.) else return(new)})
  
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]])
  metadata(new_mouse_sces) <- list(study="midbrain",study_name="midbrain_mouse")
  metadata(sces[[2]]) <- list(study="midbrain",study_name="midbrain_human")
  write_rds(new_mouse_sces,new_dataset_path[[2]])
  write_rds(sces[[2]],new_dataset_path[[1]])
}

#####convert lung cancer cells  into singcellexperiment R objects and save
preprocess_cellbench <- function(dataset){
  require(Matrix)
  require(tidyverse)
  dataset_paths <- utils.get_dataset_paths(raw_data_home,raw_datasets[[dataset]])
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]])
  dataset_names <- load(dataset_paths[[1]])
  for(i in seq_along(dataset_names)){
    study_name <- unlist(str_split(preprocesss.datasets$cellbench[[i]],"[\\.]"))[[1]]
    protocol <- unlist(str_split(preprocesss.datasets[[dataset]][[i]],"[\\._]"))[[2]]
    print(str_glue('Preprocessing {protocol} protocol'))
    sce <- eval(as.name(dataset_names[[i]]))
    colData(sce)$sampleId <- rownames(colData(sce))
    colData(sce)$label <- colData(sce)$cell_line_demuxlet
    colData(sce)$unique_id <- 1:ncol(sce)
    sce <- addPerCellQC(sce)
    sce <- sce[,which(!quickPerCellQC(sce)$discard)]
    counts(sce) <- as(counts(sce),'sparseMatrix')
    metadata(sce) <- list(study='cellbench', study_name=study_name,protocol=protocol)
    rowData(sce)$EnsembleId <- rownames(sce)
    rowData(sce)$geneName <- utils.convert2GeneSymbols(rownames(sce))
    rowData(sce)$count <- nexprs(sce,byrow=TRUE)
    write_rds(sce,new_dataset_path[[i]])
  }
}

#####human cell landscape
preprocess_human_cell_landscape <- function(dataset){
  require(tidyverse)
  hcl <- as.data.frame(read_tsv(str_glue("{raw_data_home}/human_cell_landscape.tsv")))
  rownames(hcl) <- hcl[,1]
  hcl <- hcl[,-1]
  hcl_meta <- (read_tsv(str_glue("{raw_data_home}/human_cell_landscape_meta.tsv")))
  exp_trans_hcl <- exp(hcl)-1
  data <- utils.convert_to_SingleCellExperiment(as.matrix(exp_trans_hcl),rownames(exp_trans_hcl),colnames(exp_trans_hcl),tibble(label=hcl_meta$cell_type),list(study='hcl', study_name="human_cell_landscape"))
  assay(data)$logcounts <- as.matrix(hcl)
  rowData(data)$EnsembleId <- utils.convert2EnsemblIDs(rownames(data))
  rowData(data)$geneName <- rownames(data)
  rowData(data)$count <- nexprs(data,byrow=TRUE)
  colData(data)$unique_id <- 1:ncol(data)
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]][[1]])
  write_rds(data,new_dataset_path)
  gc(T)
}

#####mouse cell atlas
preprocess_mouse_cell_atlas <- function(dataset){
  require(tidyverse)
  mca <- as.data.frame(read_tsv(str_glue("{raw_data_home}/mouse_cell_atlas.tsv")))
  rownames(mca) <- mca[,1]
  mca <- mca[,-1]
  mca_meta <- (read_tsv(str_glue("{raw_data_home}/mouse_cell_atlas_meta.tsv")))
  data <- utils.convert_to_SingleCellExperiment(as.matrix(mca),rownames(mca),colnames(mca),tibble(label=mca_meta$celltype),list(study='mca', study_name="mouse_cell_atlas"))
  rowData(data)$genes <- rownames(data)
  converted_genes <- utils.convertMouseGeneList(rowData(data)$genes)
  converted_genes <- merge(rowData(data),converted_genes,by="genes",all.x=T)
  converted_genes$EnsembleId <- utils.convert2EnsemblIDs(converted_genes$human_gene)
  converted_genes$geneName <- converted_genes$human_gene
  
  new_mouse_data <- data[converted_genes$genes,]
  new_mouse_data <- new_mouse_data[,which(!quickPerCellQC(new_mouse_data)$discard)]
  rowData(new_mouse_data) <- converted_genes
  
  new_dataset_path <- utils.get_dataset_paths(data_home,preprocesss.datasets[[dataset]][[1]])
  write_rds(data,new_dataset_path)
  gc(T)
}

##remove batch effects of datasets and saves

preprocess.remove_batch_effects <- function(sces, method="MNN"){
  sces <- utils.remove_batch_effects(sces, method) 
  return(sces)
}

# preprocess.remove_batch_effects <- function(sces,file_names){
#   sces <- utils.remove_batch_effects(sces) 
#   for(i in seq_along(sces)){
#     path <- utils.get_dataset_paths(data_home,file_names[[i]])
#     write_rds(sces[[i]],path)
#   }
# }

###intersect singlecellexperiments objects with different genes and save
preprocess.intersect_sces <- function(sces,file_names){
  paths <- utils.get_dataset_paths(data_home,file_names)
  # if(sum(unlist(purrr::map(paths,file.exists)))>=1) {
  #   print("intersected datasets already exist, skipping intersection")  
  #   return()
  # }
  sces <- utils.combine_SCEdatasets(sces,if_combined=FALSE)
  purrr::walk2(sces,paths,write_rds)
}