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
    sce <- read_rds(data_paths[[i]])
    sces[i] <- sce[!duplicated(rownames(sce)),!duplicated(colnames(sce))]
  }
  sces
}

utils.intersect_on_genes <- function(data1, data2){
  stopifnot(is(data1,"SingleCellExperiment"))
  stopifnot(is(data2,"SingleCellExperiment"))
  common_genes <- intersect(rownames(data1),rownames(data2))
  list(data1=data1[common_genes,],data2=data2[common_genes,])
}

###combine multiple singlecellexperiment into one
utils.combine_SCEdatasets <- function(sces,if_combined=TRUE,colDatas_cols=NULL){
    require(purrr)
    require(SingleCellExperiment)
    require(scater)
    intersected_genes <- purrr::reduce(purrr::map(sces,~rownames(rowData(.))),intersect)
    if(if_combined){
      if(purrr::is_null(colDatas_cols)){
        colDatas_cols <- if('batch' %in% colnames(colData(sces[[1]]))) c('batch','label') else 'label'
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
  require(tidyverse)
  colnames(data_matrix) <- cellIds
  rownames(data_matrix) <- genes
  coldata$unique_id <- 1:nrow(coldata)
  sce <- SingleCellExperiment(list(counts=data_matrix),colData=coldata,metadata=meta)
  sce <- addPerCellQC(sce)
  sce <- sce[,which(!quickPerCellQC(sce)$discard)]
  rowData(sce) <- tibble(count=nexprs(sce,byrow=TRUE),geneName=genes)
  sce
}

utils.get_dataset_paths <- function(data_home,dataset_names){
  require(tidyverse)
  paths <- purrr::map_chr(dataset_names,~paste(data_home,.,sep=''))
  paths
}

###filter out cells with all 0s reads and genes with all 0s reads
utils.filter <- function(data,filter_gene=TRUE, filter_cells=TRUE, filter_cell_type=TRUE,kept_cell_types=NULL){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  study_name <- metadata(data)$study_name[[1]]
  properties <- dataset.properties[[study_name]]
  threshold <- properties$sample_threshold
  if(!purrr::is_null(kept_cell_types)){
    print(str_glue("only keep cell types with marker genes: {kept_cell_types}"))
    data <- data[,colData(data)$label %in% kept_cell_types] ###if only keeps cell types with marker gene files
  }else if(!is_null(properties$cell_types)){
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
utils.sampler <- function(data,sample_pctg=NULL,types,column="label",sample_seed=NULL){
  if(!purrr::is_null(sample_seed))
    set.seed(sample_seed) #####ensures sample the same index by setting seed
  sample_idx <- purrr::map(types,~ which(colData(data)[[column]]==.)) %>% 
    purrr::map(~sample(.,ceiling(length(.)*sample_pctg),replace=FALSE)) 
  return(data[,unlist(sample_idx)])
}

###select different sequencing depth cells, right means take the deeper sequencing data
utils.seqDepthSelector <- function(data, low_quantile,high_quantile){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  sampled <- dplyr::group_by(dplyr::mutate(as.data.frame(colData(data)),id=colnames(data)),label) %>% filter(percent_rank(sum)>=low_quantile,percent_rank(sum)<=high_quantile)
  return(data[,sampled$id])
}

###label unassigned cells in the tibble results
utils.label_unassigned <- function(cell_types,test_results,assign=T){
  for(col in colnames(test_results)){
    if(col=="label") next
    if(sum(purrr::map_lgl(test_results[[col]],is.na)) == length(test_results[[col]])) next
    if(sum(purrr::map_lgl(test_results[[col]],purrr::is_null)) == length(test_results[[col]])) next
    if(assign){
      test_results[which(!test_results[[col]]%in% cell_types),col] <- "unassigned"
    }else{
      test_results[which(purrr::map_lgl(test_results[,col],~{!(is_numeric(.)&&(.>0)&&(!is.na(.)))})),col] <- -1L
    }
  }
  test_results
}

utils.assigned_index <- function(results){
  assign_results <- results$assign_results
  cluster_results <- results$cluster_results
  labels <- unique(assign_results$label)
  if("label" %in% names(assign_results)){
    assign_results <- dplyr::select(assign_results,-label)
    cluster_results <- dplyr::select(cluster_results,-label)
  }
  for(col in colnames(assign_results)){
    if(
      purrr::is_null(assign_results[[col]]) || 
      sum(purrr::map_lgl(assign_results[[col]],purrr::is_null))==length(assign_results[[col]]) ||
      sum(purrr::map_lgl(assign_results[[col]],is.na))==length(assign_results[[col]])
    )
      assign_results <- dplyr::select(assign_results,-col)
  }
  for(col in colnames(cluster_results)){
    if(
      purrr::is_null(cluster_results[[col]]) || 
      sum(purrr::map_lgl(cluster_results[[col]],purrr::is_null))==length(cluster_results[[col]]) ||
      sum(purrr::map_lgl(cluster_results[[col]],is.na))==length(cluster_results[[col]])
    )
      cluster_results <- dplyr::select(cluster_results,-col)
  }
  
  labeled_idx_assign <- purrr::reduce(purrr::map(assign_results,~which(. %in% labels)),intersect)
  labeled_idx_cluster <- purrr::reduce(purrr::map(cluster_results,~which(.!=-1L)),intersect)
  labeled_idx <- intersect(labeled_idx_assign, labeled_idx_cluster)
  labeled_idx
}
###select assigned cells by all methods from tibble results

utils.select_assigned <- function(results){
  labeled_idx <- utils.assigned_index(results)
  purrr::map(results,~.[labeled_idx,])
}

###select unassigned cells by all methods from tibble results
utils.select_unassigned <- function(results){
  labeled_idx <- utils.assigned_index(results)
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

######remove batch effects
utils.remove_batch_effects <- function(batches, method){
  switch(method,
         MNN = utils.batch_effects.MNN(batches),
         seurat = utils.batch_effects.seurat(batches),
         harmony = utils.batch_effects.harmony(batches),
         combat_seq = utils.batch_effects.combat_seq(batches),
         stop("Unkown experiments")
  )
}


###remove batch effects using MNN
utils.batch_effects.MNN <- function(batches){
  require(batchelor)
  require(scater)
  require(scran)
  print("using MNN to remove batch effects")
  
  batches <- do.call(multiBatchNorm,as.list(batches))
  decs <- purrr::map(batches,modelGeneVar)
  combined.decs <- do.call('combineVar',decs)
  chosen.hvgs <- combined.decs$bio > 0
  args <- list()
  args$subset.row <- chosen.hvgs
  f.out <- do.call(fastMNN,c(as.list(batches),args))
  sces_batch_free <- purrr::map(batches,~f.out[,colnames(.)])
  add_col_row_data <- function(x,y){
    colData(x)[,'label'] <- colData(y)[,'label']
    rowData(x)[,'count'] <- rowData(y)[,'count']
    rowData(x)[,'geneName'] <- rowData(y)[,'geneName']
    counts(x) <- exp(assays(x)$reconstructed)
    metadata(x) <- metadata(y)
    x <- addPerCellQC(x)
  }
  purrr::map2(sces_batch_free,batches,add_col_row_data)
}

###remove batch effects using harmony

utils.batch_effects.harmony <- function(batches){
  require(harmony)
  require(MUDAN)
  print("using harmony to remove batch effects")
  
  batch1 <- counts(batches[[1]])
  batch2 <- counts(batches[[2]])
  cd <- cbind(batch1, batch2)
  meta <- c(rep('batch1', ncol(batch1)), rep('batch2', ncol(batch2)))
  names(meta) <- c(colnames(batch1), colnames(batch2))
  meta <- factor(meta)
  mat <- MUDAN::normalizeCounts(cd, verbose=FALSE) 
  ## variance normalize, identify overdispersed genes
  matnorm.info <- MUDAN::normalizeVariance(mat, details=TRUE, verbose=FALSE) 
  ## log transform
  matnorm <- log10(matnorm.info$mat+1) 
  ## 30 PCs on overdispersed genes
  harmonized <- HarmonyMatrix(matnorm, meta, do_pca = FALSE, verbose = FALSE)
  # harmonized[harmonized<0] = 0
  corrected_batch_1 <- t(harmonized[names(meta[meta=="batch1"]),])
  corrected_batch_1 <- utils.convert_to_SingleCellExperiment(corrected_batch_1,rownames(batch1),colnames(batch1),colData(batches[[1]]),metadata(batches[[1]]))
  corrected_batch_2 <- t(harmonized[names(meta[meta=="batch2"]),])
  corrected_batch_2 <- utils.convert_to_SingleCellExperiment(corrected_batch_2,rownames(batch2),colnames(batch2),colData(batches[[2]]),metadata(batches[[2]]))
  return(list(corrected_batch_1,corrected_batch_2))
}

##########
utils.batch_effects.combat_seq <- function(batches){
  require(sva)
  print("using combat-seq to remove batch effects")
  batch1 <- as.matrix(counts(batches[[1]]))
  batch2 <- as.matrix(counts(batches[[2]]))
  cd <- cbind(batch1, batch2)
  meta <- c(rep(1, ncol(batch1)), rep(2, ncol(batch2)))
  names(meta) <- c(colnames(batch1), colnames(batch2))
  meta <- factor(meta)
  adjusted <- ComBat_seq(cd, batch=meta, group=NULL)
  corrected_batch_1 <- adjusted[,names(meta[meta==1])]
  corrected_batch_1 <- utils.convert_to_SingleCellExperiment(corrected_batch_1,rownames(batch1),colnames(batch1),colData(batches[[1]]),metadata(batches[[1]]))
  corrected_batch_2 <- adjusted[,names(meta[meta==2])]
  corrected_batch_2 <- utils.convert_to_SingleCellExperiment(corrected_batch_2,rownames(batch2),colnames(batch2),colData(batches[[2]]),metadata(batches[[2]]))
  return(list(corrected_batch_1,corrected_batch_2))
}


####append one singlecellexperiment onto another
utils.append_sce <- function(sce1,sce2){
  require(scater)
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
    # if(!file.exists(str_glue("{marker_home}/{generated_marker_gene_file}.RDS"))){
    print(str_glue("Generating {marker_home}/{generated_marker_gene_file}.RDS"))
    markers_mat <- markergene.generate_marker_genes(marker_gene_method,data,str_glue("{marker_home}/{generated_marker_gene_file}"))
    i <- 0;
    while(is.na(rownames(markers_mat))){
      if(i >= 10) break;
      i = i+1
      print(str_glue("in seurat generating marker gene, markers mat is zero length, adjusting logfc_threshold
            from {markergene.config.seurat$logfc_threshold} to {markergene.config.seurat$logfc_threshold/2}"))
      markergene.config.seurat$logfc_threshold <<- markergene.config.seurat$logfc_threshold/2
      markers_mat <- markergene.generate_marker_genes(marker_gene_method,data,str_glue("{marker_home}/{generated_marker_gene_file}"))
    }  
    # }
    # else{
    #   print(str_glue("Loading {marker_home}/{generated_marker_gene_file}.RDS"))
    #   markers_mat <- read_rds(str_glue("{marker_home}/{generated_marker_gene_file}.RDS"))
    # }
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
    gene_names <- rownames(data) 
    matchidx <- match(markers[,1], gene_names)
    markers_mat <- markers_mat[!is.na(matchidx),]
  }
  list(markers_mat=markers_mat,matchidx=matchidx)
}

utils.update_batch_effects_free_config <- function(experiment){
  # experiments.methods[[experiment]]$cluster <<- experiments.methods[[experiment]]$cluster_batch_free
  print(str_glue("batch effects free cluster methods are: {experiments.methods[[experiment]]$cluster}"))
  for(m in experiments.methods[[experiment]]$cluster){
    if(exists(str_glue("methods.config.{m}.batch_free"))){
      assign(str_glue("methods.config.{m}"), get(str_glue("methods.config.{m}.batch_free")),envir = .GlobalEnv)
    }
  }
  
  # experiments.methods[[experiment]]$assign <<- experiments.methods[[experiment]]$assign_batch_free
  print(str_glue("batch effects free assign methods are: {experiments.methods[[experiment]]$assign}"))
  for(m in experiments.methods[[experiment]]$assign){
    if(exists(str_glue("methods.config.{m}.batch_free"))){
      assign(str_glue("methods.config.{m}"),get(str_glue("methods.config.{m}.batch_free")),envir = .GlobalEnv)
    }
  }
  
  # experiments.methods[[experiment]]$marker_gene_assign <<- experiments.methods[[experiment]]$marker_gene_assign_batch_free
  # print(str_glue("batch effects free marker gene methods are: {experiments.methods[[experiment]]$marker_gene_assign}"))
  # for(m in experiments.methods[[experiment]]$marker_gene_assign){
  #   if(exists(str_glue("methods.config.{m}.batch_free"))){
  #     assign(str_glue("methods.config.{m}"),get(str_glue("methods.config.{m}.batch_free")),envir = .GlobalEnv)
  #   }
  # }
}

utils.try_catch_method_error <- function(code, silent=FALSE){
  tryCatch(code, error=function(c) {
    msg <- conditionMessage(c)
    if (!silent) warning(c)
    print(str_glue("error occured {msg}"))
    error(logger, str_glue("error occured {msg}"))
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


utils.convertCellTypes <- function(from_type, type_map){
  for(i in seq_along(rownames(type_map))){
    from_type[which(from_type==type_map$from[[i]])]=type_map$to[[i]]
  }
  from_type
}

utils.assign_test_cluster_label_from_training <- function(train_dataset_label,train_dataset_num,test_dataset_num){
  require(modeest)
  unique_cluster_num <- union(unique(train_dataset_num),unique(test_dataset_num))
  train_label_num <- tibble(label=train_dataset_label,num=train_dataset_num)
  test_label_num <- tibble(num=test_dataset_num) 
  num_to_label <- train_label_num %>% 
                  dplyr::group_by(num) %>%
                  summarize(label=mfv(label)[[1]])
  
  test_label_num <- dplyr::left_join(test_label_num,num_to_label,by="num")
  test_label_num[is.na(test_label_num$num),'label'] <- "unknown"
  test_label_num[is.na(test_label_num$num),'num'] <- -1
  test_label_num
}

# utils.remove_duplicated_sampleid <- function(...){
#   datasets <- list(...)
#   stopifnot(0==sum(purrr::map_lgl(datasets,~{!is(.,"SingleCellExperiment")})))
#   sampleids <- purrr::reduce(purrr::map(datasets,~{rownames(colData(.))}),union)
#   for(i in 1:length(datasets)){
#     sampleid <- rownames(colData(datasets[[i]]))
#     duplicated_sampleid <- purrr:map_lgl(sampleid,~{.%in%sampleids})
#   }
# }

utils.clean_marker_files <- function(){
  marker_files <- c(dir(path=marker_home,  pattern=".RDS",full.names = T),dir(path=marker_home,  pattern="Garnett*",full.names = T))
  file.remove(marker_files)
}

utils.get_logger <- function(level="DEBUG", file){
  require(log4r)
  if(!dir.exists(log_home)) dir.create(log_home)
  log_file <- str_glue("{log_home}/{file}")
  if(file.exists(log_file)) file.remove(log_file)
  file_logger <- logger(level, appenders = file_appender(log_file))
  file_logger
}

utils.remove_files <- function(file_names){
  paths <- utils.get_dataset_paths(data_home,file_names)
  file.remove(paths)
}


utils.manhattan_dist <- function(a,b){
  dist <- abs(a-b)
  dist <- sum(dist)
  return(dist)
}


utils.calculate_sampling_pctg <- function(dataset, target_num, exp_config, if_train){
  total_num <- length(colnames(dataset))
  if(exp_config$use_inter_dataset){
    print("use inter dataset")
    target_sampling_pctg <- min(1,target_num/total_num)
  }else{
    print("use intra dataset")
    if(exp_config$cv){
      print("using CV")
      target_sampling_pctg <- min(1,target_num/total_num)
    }else{
      print("using same train test dataset")
      if(if_train){
        target_sampling_pctg <- min(0.8,target_num/total_num)
      }else{
        target_sampling_pctg <- min(1,target_num/total_num)
      }
    }
  }
  target_sampling_pctg
}


utils.calc_entropy <- function(probs){
  entropy <- sum(map_dbl(probs,~{-log(.)*.}))
  entropy
}

utils.get_train_test_types <- function(train_dataset,test_dataset){
  train_dataset_name <- dataset.name.map1[[train_dataset]]
  test_dataset_name <- dataset.name.map1[[test_dataset]]
  if(!purrr::is_null(dataset.properties[[train_dataset_name]]$cell_types)){
    train_type <- as.list(dataset.properties[[train_dataset_name]]$cell_types)
  }else{
    train_data <- utils.load_datasets(train_dataset) %>% utils.filter(filter_gene=F, filter_cells=TRUE, filter_cell_type=TRUE)

    train_type <- unique(colData(train_data)$label) 
  }
  
  if(!purrr::is_null(dataset.properties[[test_dataset_name]]$cell_types)){
    test_type <- as.list(dataset.properties[[test_dataset_name]]$cell_types)
  }else{
    test_data <- utils.load_datasets(test_dataset) %>% utils.filter(filter_gene=F, filter_cells=TRUE, filter_cell_type=TRUE)
    test_type <- unique(colData(test_data)$label) 
  }
  return(list(train_type=train_type,test_type=test_type))
}

utils.cor.stable <- function (x, y, method="pearson", ...) {
  omit1 <- which(apply(x, 2, sd) == 0)
  omit2 <- which(apply(y, 2, sd) == 0)
  if (length(omit1) > 0 && length(omit2) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[-omit1,-omit2] = cor(x[,-omit1], y[,-omit2], method=method, ...)
  } else if (length(omit1) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[-omit1,] = cor(x[,-omit1], y, method=method, ...)
  } else if (length(omit2) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[,-omit2] = cor(x, y[,-omit2], method=method, ...)
  } else {
    r = cor(x, y, method=method, ...)
  }
  
  return(r)
}

utils.medianMatrix <- function(mat,groups) {
  require(matrixStats)
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if(sum(A)==1) {
      mat[,A]
    } else {
      rowMedians(mat[,A],na.rm=T)
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}


utils.quantileMatrix <- function(mat,groups,q=0.5) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if (nrow(mat)==1) {
      quantile(mat[A],na.rm=T,probs=q)
    } else {
      if(sum(A)==1) {
        mat[,A]
      } else {
        rowQuantiles(mat[,A],na.rm=T,probs=q)
      }
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}
