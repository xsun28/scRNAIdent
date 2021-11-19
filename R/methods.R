source("R/methods_config.R")
run_assign_methods <- function(method,train_data, test_data,exp_config){
  # source("R/methods_config.R")
  switch(method,
         scmap_cluster = assign.scmap_cluster(train_data, test_data, exp_config),
         scmap_cell = assign.scmap_cell(train_data, test_data, exp_config),
         chetah = assign.chetah(train_data, test_data),
         cellassign = assign.cellassign(train_data,test_data,exp_config),
         garnett = assign.garnett(train_data,test_data,exp_config),
         singlecellnet = assign.singlecellnet(train_data,test_data, exp_config),
         singleR = assign.singleR(train_data,test_data, exp_config),
         seurat_mapping = assign.seurat(train_data,test_data, exp_config),
         stop("No such assigning method")
  )
  
}

run_cluster_methods <- function(method,data,exp_config){
  switch(method,
         sc3 = cluster.sc3(data,exp_config),
         liger = cluster.liger(data,exp_config),
         seurat_clustering = cluster.seurat(data,exp_config),
         cidr = cluster.cidr(data,exp_config),
         tscan = cluster.tscan(data,exp_config),
         cidr = cluster.cidr(data,exp_config),
         monocle3 = cluster.monocle3(data,exp_config),
         pcaReduce = cluster.pcaReduce(data,exp_config),
         raceID3 = cluster.raceID3(data,exp_config),
         same_clustering = cluster.same_clustering(data,exp_config),
         sharp = cluster.sharp(data,exp_config),
         stop("No such cluster method")
  )
}

######assigning using suerat
assign.seurat <- function(train_data, test_data, exp_config) {
  require(Seurat)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  train_study_name <- metadata(train_data)$study_name
  test_study_name <- metadata(test_data)$study_name
  train_data_prop <- dataset.properties[[train_study_name]]
  test_data_prop <- dataset.properties[[test_study_name]]
  
  train_cnts <- counts(train_data)
  test_cnts <- counts(test_data)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.seurat.batch_free")) methods.config.seurat else methods.config.seurat.batch_free
  ## PCA dimention redcution
  
  seuset <- CreateSeuratObject(train_cnts, project='supervised')
  if(!train_data_prop$logged)
    seuset <- NormalizeData(object = seuset)
  seuset <- FindVariableFeatures(object = seuset,nfeatures = m_config[['nfeatures']])
  seuset <- ScaleData(object = seuset)
  seuset <- RunPCA(object = seuset)
  seuset$celltype <- train_data$label
  seuset1 <- CreateSeuratObject(test_cnts, project='supervised')
  if(!test_data_prop$logged)
    seuset1 <- NormalizeData(seuset1, verbose = FALSE)
  seuset1 <- FindVariableFeatures(seuset1, selection.method = "vst", nfeatures =  m_config[['nfeatures']],
                                  verbose = FALSE)
  anchors <- FindTransferAnchors(reference = seuset, query = seuset1, dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = seuset$celltype,
                              dims = 1:30)
  query <- AddMetaData(seuset1, metadata = predictions)
  
  unname(query$predicted.id)
}



####assigning using scmap cluster
assign.scmap_cluster <- function(train_data, test_data, exp_config){
  require(scmap)
  require(scater)
  train_study_name <- metadata(train_data)$study_name
  test_study_name <- metadata(test_data)$study_name
  train_data_prop <- dataset.properties[[train_study_name]]
  test_data_prop <- dataset.properties[[test_study_name]]
  if(experiment %in% c("cell_number","sequencing_depth"))
    saved_trained_path <- str_glue("{pretrained_home}/{experiment}_scmap_cluster_trained.RData")
  if(experiment %in% c("cell_number","sequencing_depth") && exp_config$trained){
    stopifnot(is(test_data,"SingleCellExperiment"))
    print(str_glue("scmap_cluster already trained, load from {saved_trained_path}"))
    load(saved_trained_path)
  }else{
    stopifnot(is(train_data,"SingleCellExperiment"))
    stopifnot(is(test_data,"SingleCellExperiment"))
    counts(train_data) <- as.matrix(counts(train_data))
    if(!train_data_prop$logged)
      train_data <- logNormCounts(train_data)
    rowData(train_data)$feature_symbol <- rownames(train_data)
    train_data <- selectFeatures(train_data,n_features=methods.config.scmap[['nfeatures']]) %>%
      indexCluster(cluster_col = "label")
    index_list <- list(metadata(train_data)$scmap_cluster_index)
    if(experiment %in% c("cell_number","sequencing_depth")){
      print(str_glue("saving {experiment} trained scmap_cluster to {saved_trained_path}"))
      save(index_list,file=saved_trained_path)
    }
  }
  
  counts(test_data) <- as.matrix(counts(test_data))
  if(!test_data_prop$logged)
    test_data <- logNormCounts(test_data)
  rowData(test_data)$feature_symbol <- rownames(test_data)
  results <- scmapCluster(projection=test_data, index_list=index_list)
  return(results$combined_labs)
}

####assigning using scmap cell
assign.scmap_cell <- function(train_data, test_data, exp_config){
  require(scmap)
  require(scater)
  train_study_name <- metadata(train_data)$study_name
  test_study_name <- metadata(test_data)$study_name
  train_data_prop <- dataset.properties[[train_study_name]]
  test_data_prop <- dataset.properties[[test_study_name]]
  if(experiment %in% c("cell_number","sequencing_depth"))
    saved_trained_path <- str_glue("{pretrained_home}/{experiment}_scmap_cell_trained.RData")
  if(experiment %in% c("cell_number","sequencing_depth") && exp_config$trained){
    stopifnot(is(test_data,"SingleCellExperiment"))
    print(str_glue("scmap_cell already trained, load from {saved_trained_path}"))
    load(saved_trained_path)
  }else{
    stopifnot(is(train_data,"SingleCellExperiment"))
    stopifnot(is(test_data,"SingleCellExperiment"))
    counts(train_data) <- as.matrix(counts(train_data))
    if(!train_data_prop$logged)
      train_data <- logNormCounts(train_data)
    rowData(train_data)$feature_symbol <- rownames(train_data)
    train_data <- selectFeatures(train_data,n_features=methods.config.scmap[['nfeatures']]) %>%
      indexCell()
    index_list <- list(metadata(train_data)$scmap_cell_index)
    train_labels <- list(as.character(colData(train_data)$label))
    if(experiment %in% c("cell_number","sequencing_depth")){
      print(str_glue("saving {experiment} trained scmap_cell to {saved_trained_path}"))
      save(index_list,train_labels,file=saved_trained_path)
    }
  }
  
  counts(test_data) <- as.matrix(counts(test_data))
  if(purrr::is_null(methods.config.scmap[['seed']]))
    set.seed(1)
  else
    set.seed(methods.config.scmap[['seed']])
  if(!test_data_prop$logged)
    test_data <- logNormCounts(test_data)
  rowData(test_data)$feature_symbol <- rownames(test_data)
  cell_results <- scmapCell(projection=test_data, index_list = index_list)
  clusters_results <- scmapCell2Cluster(
    cell_results,train_labels,threshold = methods.config.scmap[['threshold']]
  )
  clusters_results$combined_labs
}

###assigning using chetah
assign.chetah <- function(train_data, test_data){
  require(SingleCellExperiment)
  require(CHETAH)
  source("R/CHETAH.R")
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  colData(train_data)$celltypes <- colData(train_data)$label
  
  #run classifier
  test_data <- CHETAHclassifier(input = test_data, ref_cells = train_data)
  
  ## Extract celltypes:
  #CHETAH generates generic names for the intermediate types: “Unassigned” for cells that are classified to the very first node, 
  #and “Node1”, “Node2”, etc for the additional nodes
  unname(test_data$celltype_CHETAH)
}
# getAnywhere("CHETAHclassifier")
# rlang::env_unlock(env = asNamespace('CHETAH'))
# rlang::env_binding_unlock(env = asNamespace('CHETAH'))
# assign('CHETAHclassifier', assign.CHETAHclassifier, envir = asNamespace('CHETAH'))
# rlang::env_unlock(env = as.environment('package:CHETAH'))
# rlang::env_binding_unlock(env = as.environment('package:CHETAH'))
# assign('CHETAHclassifier', assign.CHETAHclassifier, envir = as.environment('package:CHETAH'))
# rlang::env_binding_lock(env = asNamespace('CHETAH'))
# rlang::env_lock(asNamespace('CHETAH'))



#######assigning using garnet
assign.garnett <- function(train_data,test_data,exp_config){
  require(garnett)
  require(monocle)
  require(SingleCellExperiment)
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  print("start method garnett")
  if(experiment %in% c("cell_number","sequencing_depth")){
    saved_trained_path <- str_glue("{pretrained_home}/{experiment}_garnett_trained.RData")
    if(exp_config$trained){
      stopifnot(is(test_data,"SingleCellExperiment"))
      print(str_glue("garnett already trained, load from {saved_trained_path}"))
      load(saved_trained_path)
    }
  }
  
  stopifnot(is(test_data,"SingleCellExperiment"))
  stopifnot(is(train_data,"SingleCellExperiment"))
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.garnett.batch_free")) methods.config.garnett else methods.config.garnett.batch_free
  study <- metadata(train_data)$study[[1]]
  study_name <- metadata(train_data)$study_name[[1]]
  gene_name_type <- dataset.properties[[study_name]]$gene_name_type
  pretrained_classifier_name <- m_config[[study_name]]$pretrained_classifier
  
  if(is_null(pretrained_classifier_name)){
    if(experiment %in% c("celltype_structure")){
      marker_file_path <- str_glue("Garnett_{study_name}_marker_{gene_name_type}_{exp_config$level}.txt")
    }else{
      marker_file_path <- m_config[[study_name]]$marker_file_path
    }
    if(purrr::is_null(marker_file_path)||!file.exists(str_glue("{marker_home}/{marker_file_path}"))){
      print("Garnett marker file not exist, generating marker file...")
      if(purrr::is_null(marker_file_path))
        marker_file_path <- str_glue("Garnett_{study_name}_marker_{gene_name_type}.txt")
      marker_gene_file <- exp_config$marker_gene_file
      if(experiment %in% c("celltype_structure")){
        generated_marker_gene_file <- str_glue("{study_name}_markergene_{m_config[[study_name]]$marker_gene_method}_{gene_name_type}_{exp_config$level}")
      }else{
        generated_marker_gene_file <- str_glue("{study_name}_markergene_{m_config[[study_name]]$marker_gene_method}_{gene_name_type}")
      }
      marker_gene_method <- m_config[[study_name]]$marker_gene_method
      check_results <- utils.check_marker_genes(train_data,marker_gene_file, generated_marker_gene_file,marker_gene_method) 
      markers_mat <- check_results$markers_mat
      matchidx <- check_results$matchidx
      # if(!file.exists(str_glue("{marker_home}/{marker_file_path}"))){
      print(str_glue("{marker_file_path} generating...."))
      assign.garnett.generate_marker_file(markers_mat,marker_file_path,gene_name_type)
      # }
    }else{
      print(str_glue("{marker_file_path} exists,skipping generating..."))
    }
    
    set.seed(260)
    db <- org.Hs.eg.db
    if("species" %in% colnames(colData(train_data))){
      if(colData(train_data)$species[[1]]=="h"){
        db <- org.Hs.eg.db
      }else if(colData(train_data)$species[[1]]=="m"){
        db <- org.Hs.eg.db
      }else{
        stop("unkown train species")
      }
    }
    num_unknown <- if(purrr::is_null(m_config[[study_name]]$num_unkown)) 100 else m_config[[study_name]]$num_unknown
    min_observations <- if(purrr::is_null(m_config[[study_name]]$min_observations)) 8 else m_config[[study_name]]$min_observations
    if(!(experiment %in% c("cell_number","sequencing_depth") && exp_config$trained)){
      cds_train <- assign.garnett.process_data(train_data, gene_name_type)
      classifier <- train_cell_classifier(cds = cds_train,
                                          marker_file = str_glue("{marker_home}/{marker_file_path}"),
                                          db=db,
                                          min_observations=min_observations,
                                          cds_gene_id_type = gene_name_type,
                                          num_unknown = num_unknown,
                                          marker_file_gene_id_type = gene_name_type)
      write_rds(classifier,str_glue("{pretrained_home}/pretrained_{study_name}_{gene_name_type}_classifier.rds"))
      if((experiment %in% c("cell_number","sequencing_depth") && !exp_config$trained)){
        print(str_glue("saving {experiment} trained garnett classifier to {saved_trained_path}"))
        save(classifier,file=saved_trained_path)
      }
    }
  }else{
    classifier_path <- str_glue("{pretrained_home}/{pretrained_classifier_name}")
    classifier <- readRDS(classifier_path)
  }
  cds_test <- assign.garnett.process_data(test_data, gene_name_type)
  
  cds_test <- classify_cells(cds_test, classifier,
                             db = db,
                             cluster_extend = TRUE,
                             cds_gene_id_type = gene_name_type)
  utils.clean_marker_files()
  cds_test$cell_type
}


assign.garnett.process_data <- function(data, gene_name_type){
  mat <- counts(data)
  if(gene_name_type=="SYMBOL"){
    geneData <- as.tibble(rowData(data))
    geneData$fullGeneName <- rownames(data) 
    geneData <- dplyr::group_by(geneData,fullGeneName) %>%
      dplyr::filter(count==max(count))
    mat <- counts(data[geneData$fullGeneName])
    rownames(mat) <- geneData$fullGeneName
    # geneData$gene_short_name <- geneData$geneName
    geneData <- as.data.frame(geneData)
    rownames(geneData) <- geneData$fullGeneName
  }else if(gene_name_type=="ENSEMBL"){
    geneData <- rowData(data)
    mat <- counts(data)
    rownames(mat) <- geneData$EnsembleId
    # geneData$gene_short_name <- geneData$geneName
    geneData <- as.data.frame(geneData)
    rownames(geneData) <- geneData$EnsembleId
  }else{
    stop(str_glue("Unknown gene name type: {gene_name_type}"))
  }
  
  
  cellData <- as.data.frame(colData(data))
  # create a new CDS object
  phenoData <- new("AnnotatedDataFrame", data = cellData)
  featureData <- new("AnnotatedDataFrame", data = geneData)
  cds <- newCellDataSet(as(mat, "dgCMatrix"),
                        phenoData = phenoData,
                        featureData = featureData)
  
  # generate size factors for normalization later
  cds <- estimateSizeFactors(cds)
}

assign.garnett.generate_marker_file <- function(markers_mat,marker_file_path,gene_name_type){
  if(length(markers_mat)==0) throw("markers_mat failed to generate in Garnett")
  if(gene_name_type == "SYMBOL"){
    gene_names <- rownames(markers_mat) %>%
      purrr::map_chr(~{names <- str_split(.,"\\t")[[1]]
      return(names[which(!str_detect(names,"^ENSG0+"))])
      })
  }else{
    gene_names <- rownames(markers_mat) %>%
      purrr::map_chr(~{names <- str_split(.,"\\t")[[1]]
      return(names[which(str_detect(names,"^ENSG0+"))])
      })
  }
  rownames(markers_mat) <- gene_names
  cell_types <- colnames(markers_mat)
  marker_file_path <- str_glue("{marker_home}/{marker_file_path}")
  print(str_glue("Generating Garnett marker file: {marker_file_path}"))
  for(t in cell_types){
    genes <- names(markers_mat[which(markers_mat[,t] > 0),t])
    if(length(genes)>0){
      
      line <- str_glue(">{t}")
      write(line,file=marker_file_path, append = TRUE)
      
      line <- str_c(genes,collapse=", ")
      
      line <- str_c("expressed: ",line)
      write(line,file=marker_file_path,append = TRUE)
    }
  }
}




###assigning using cellassign
assign.cellassign <- function(train_data,test_data,exp_config){
  require(tensorflow)
  require(cellassign)
  require(scran)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.cellassign.batch_free")) methods.config.cellassign else methods.config.cellassign.batch_free
  marker_gene_file <- exp_config$marker_gene_file
  study_name <- metadata(train_data)$study_name[[1]]
  gene_name_type <- dataset.properties[[study_name]]$gene_name_type
  
  if(experiment %in% c("celltype_structure")){
    generated_marker_gene_file <- str_glue("{study_name}_markergene_{m_config$marker_gene_method}_{gene_name_type}_{exp_config$level}")
  }
  else{
    generated_marker_gene_file <- str_glue("{study_name}_markergene_{m_config$marker_gene_method}_{gene_name_type}")
  }
  marker_gene_method <- m_config$marker_gene_method
  check_results <- utils.check_marker_genes(train_data,marker_gene_file, generated_marker_gene_file,marker_gene_method) 
  markers_mat <- check_results$markers_mat
  if(length(markers_mat)==0) throw("markers_mat failed to generate in cellassign")
  matchidx <- check_results$matchidx
  counts(test_data) <- as.matrix(counts(test_data))
  test_data$celltypes <- colData(test_data)$label
  s <- computeSumFactors(test_data) %>% 
    sizeFactors()
  
  fit <- cellassign(exprs_obj = test_data[na.omit(matchidx),], 
                    marker_gene_info = markers_mat, 
                    s = s, 
                    learning_rate = m_config$learning_rate, 
                    shrinkage = m_config$shrinkage,
                    verbose = FALSE)
  utils.clean_marker_files()
  fit$cell_type
}

######assigning using singlecellnet
assign.singlecellnet <- function(train_data, test_data, exp_config){
  require(singleCellNet)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  #' extract sampTab and expDat sce object into regular S3 objects
  #' @param sce_object
  #' @param exp_type
  #' @param list
  #' @export
  extractSCE <- function(sce_object, exp_type = "counts"){
    #extract metadata
    sampTab = as.data.frame(colData(sce_object, internal = TRUE))
    sampTab$sample_name = rownames(sampTab)
    
    #extract expression matrix
    if(exp_type == "counts"){
      expDat = counts(sce_object)
    }
    
    if(exp_type == "normcounts"){
      expDat = normcounts(sce_object)
    }
    
    if(exp_type == "logcounts"){
      expDat = logcounts(sce_object)
    }
    
    return(list(sampTab = sampTab, expDat = expDat))
  }
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.singlecellnet.batch_free")) methods.config.singlecellnet else methods.config.singlecellnet.batch_free
  set.seed(100) #can be any random seed number
  train_scefile <- extractSCE(train_data, exp_type = "counts") 
  train_metadata <- train_scefile$sampTab
  train_expdata <- train_scefile$expDat
  test_scefile <- extractSCE(test_data, exp_type = "counts") 
  test_metadata <- test_scefile$sampTab
  test_expdata <- test_scefile$expDat
  
  if(m_config$cross_species){
    common_gene_file <- str_glue("{data_home}{m_config$common_gene_file}")
    oTab <- utils_loadObject(fname = "../data/human_mouse_genes_Jul_24_2018.rda")
    aa <- csRenameOrth(expQuery = test_expdata, expTrain = train_expdata, orthTable = oTab)
    test_expdata <- aa[['expQuery']]
    train_expdata <- aa[['expTrain']]
  }else{
    commonGenes<-intersect(rownames(train_expdata), rownames(test_expdata))
    train_expdata <- train_expdata[commonGenes, ]
    test_expdata <- test_expdata[commonGenes, ]
  }
  ncells <- if(purrr::is_null(m_config$ncells)) 100 else m_config$ncells
  nTopGenes <- if(purrr::is_null(m_config$nTopGenes)) 10 else m_config$nTopGenes
  nRand <- if(purrr::is_null(m_config$nRand)) 70 else m_config$nRand
  nTrees <- if(purrr::is_null(m_config$nTrees)) 1000 else m_config$nTrees
  nTopGenePairs <- if(purrr::is_null(m_config$nTopGenePairs)) 25 else m_config$nTopGenePairs
  
  class_info<-scn_train(stTrain = train_metadata, expTrain = train_expdata, nTopGenes = nTopGenes, nRand = nRand, nTrees = nTrees, nTopGenePairs = nTopGenePairs, 
                        dLevel = "label", colName_samp = "sample_name")
  #predict
  pred_results <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_expdata, nrand = 0)
  pred_labels <- assign_cate(classRes = pred_results, sampTab = test_metadata, cThresh = 0.5)
  pred_labels$category
}

#### assigning using singleR
assign.singleR <- function(train_data, test_data, exp_config=NULL){
  require(SingleR)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.singleR.batch_free")) methods.config.singleR else methods.config.singleR.batch_free
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  train_data <- logNormCounts(train_data) 
  test_data <- logNormCounts(test_data) 
  preds <- SingleR(test=test_data, ref=train_data, labels=train_data$label, de.method="wilcox")
  preds$labels
}



### clustering using Seurat
cluster.seurat <- function(data,exp_config) {
  require(Seurat)
  stopifnot(is(data,"SingleCellExperiment"))
  cnts <- counts(data)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.seurat.batch_free")) methods.config.seurat else methods.config.seurat.batch_free
  ## PCA dimention redcution
  seuset <- CreateSeuratObject(cnts, project='unsupervised')
  seuset <- NormalizeData(object = seuset)
  seuset <- FindVariableFeatures(object = seuset,nfeatures = m_config[['nfeatures']])
  seuset <- ScaleData(object = seuset)
  seuset <- RunPCA(object = seuset)
  seuset <- FindNeighbors(object = seuset,dims = 1:m_config[['pc_dims']])
  seuset <- FindClusters(object = seuset, resolution = m_config[['resolution']])
  as.integer(unname(seuset$seurat_clusters))
}

### TSCAN
cluster.tscan <- function(data,exp_config) {
  require(TSCAN)
  stopifnot(is(data,"SingleCellExperiment"))
  cnts <- counts(data)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.tscan.batch_free")) methods.config.tscan else methods.config.tscan.batch_free
  minexpr_value = m_config[['minexpr_value']]
  minexpr_percent = m_config[['minexpr_percent']] 
  cvcutoff=m_config[['cvcutoff']]
  dataset_prop <- dataset.properties[[metadata(data)$study_name]]
  procdata <- preprocess(as.matrix(cnts),takelog = !dataset_prop$logged,minexpr_value = minexpr_value, minexpr_percent =minexpr_percent, cvcutoff=cvcutoff)
  i = 0;
  while(dim(procdata)[1]<=10){
    if(i > 10) return(NULL)
    i = i + 1
    print(str_glue({"minexpr_value={minexpr_value},minexpr_percent = {minexpr_percent},cvcutoff = {cvcutoff}
      preprocess data in tscan results in 0 dim data, adjusting paramters... "}))
    minexpr_value = minexpr_value/2.0
    minexpr_percent = minexpr_percent/2.0
    cvcutoff = cvcutoff/2.0
    procdata <- preprocess(as.matrix(cnts),takelog = F,minexpr_value = minexpr_value, minexpr_percent =minexpr_percent, cvcutoff=cvcutoff)
  }
  
  if(exp_config$known_cluster_num){
    K = length(unique(colData(data)$label))
  }else{
    K = m_config[['k']]
  }
  
  if(!purrr::is_null(K) )
    ret <- tryCatch(exprmclust(procdata, clusternum=K), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured {msg}"))
      error(logger, str_glue("error occured in tscan {msg}"))
      structure(msg, class = "try-error")
    })
  else
    ret <- tryCatch(exprmclust(procdata), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured {msg}"))
      error(logger, str_glue("error occured in tscan {msg}"))
      structure(msg, class = "try-error")
    })
  i = 0
  K = if(length(K)>1) K[1] else  K/2
  while(inherits(ret,"try-error")){ ###most try 10 iterations
    if(i > 10 || K<=1 || exp_config$known_cluster_num) return(NULL)
    i = i+1
    K = ceiling(K/2)
    print(str_glue("in tscan:{ret}"))
    print(str_glue("new k: {K}"))
    ret <- tryCatch(exprmclust(procdata, clusternum=K), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured {msg}"))
      error(logger, str_glue("error occured in tscan {msg}"))
      structure(msg, class = "try-error")
    })
  }
  
  as.integer(unname(ret$clusterid))
}

### SC3
cluster.sc3 <- function(data,exp_config) {
  require(SC3)
  stopifnot(is(data,"SingleCellExperiment"))
  counts(data) <- as.matrix(counts(data))
  data <- logNormCounts(data)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.sc3.batch_free")) methods.config.sc3 else methods.config.sc3.batch_free
  gene_filter <- m_config$gene_filter
  print(str_glue("gene filter for sc3 is {gene_filter}"))
  rowData(data)$feature_symbol <- rownames(data)
  data <- data[!duplicated(rowData(data)$feature_symbol),]
  data <- sc3_prepare(data,gene_filter = gene_filter)
  if(exp_config$known_cluster_num){
    k <- length(unique(colData(data)$label))
  }else{
    k <- m_config$k
  }
  if(is_null(k)){
    ret <- tryCatch(sc3_estimate_k(data), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured {msg}"))
      error(logger, str_glue("error occured in sc3 {msg}"))
      structure(msg, class = "try-error")
    })
    if(!inherits(ret,"try-error")){
      k <- metadata(ret)$sc3$k_estimation
    }else{
      k <- 8
    }
  }
  data <- sc3_calc_dists(data)
  data <- sc3_calc_transfs(data)
  data <- sc3_kmeans(data, ks = k)
  data <- sc3_calc_consens(data)
  col_data <- colData(data)
  colTb = as.vector(col_data[ , grep("sc3_", colnames(col_data))])
  as.integer(colTb)
}

#####liger
cluster.liger <- function(data,exp_config){
  require(rliger)
  num.cores <- detectCores()-2
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.liger.batch_free")) methods.config.liger else methods.config.liger.batch_free
  stopifnot(is(data,"SingleCellExperiment"))
  data_list <- list(counts(data))
  names(data_list) <- metadata(data)$study[[1]]
  liger_data <- createLiger(data_list)
  liger_data <- rliger::normalize(liger_data)
  liger_data <- selectGenes(liger_data, var.thresh = c(0.3, 0.875), do.plot = F)
  liger_data <- scaleNotCenter(liger_data)
  if(m_config$suggestK==TRUE){
    print("suggesting K")
    ret <- tryCatch(suggestK(liger_data, num.cores = num.cores, gen.new = T, plot.log2 = F,nrep = 5), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured {msg}"))
      error(logger, str_glue("error occured in liger {msg}"))
      structure(msg, class = "try-error")
    })
    if(!inherits(ret,"try-error")){
      k.suggest <- ret
    }else{
      info(logger,"in liger can not suggest K")
      k.suggest <- if(purrr::is_null(m_config$k.suggest))  25 else m_config$k.suggest
    }
  }else{
    k.suggest <- if(purrr::is_null(m_config$k.suggest))  25 else m_config$k.suggest
  }
  thresh <- if(purrr::is_null(m_config$thresh)) 5e-5 else m_config$thresh
  lambda <- if(purrr::is_null(m_config$lambda)) 5 else m_config$lambda
  resolution <- if(purrr::is_null(m_config$resolution)) 1.0 else m_config$resolution
  
  ret <- tryCatch(optimizeALS(liger_data, k=k.suggest, thresh = thresh, lambda=lambda,nrep = 3), error=function(c) {
    msg <- conditionMessage(c)
    print(str_glue("error occured {msg}"))
    error(logger, str_glue("error occured in liger {msg}"))
    structure(msg, class = "try-error")
  })
  i = 0
  while(inherits(ret,"try-error")){ ###most try 10 iterations
    if(i > 10) return(NULL)
    i = i+1
    print(str_glue("in liger:{ret}"))
    k.suggest = k.suggest -5
    print(str_glue("new k: {k.suggest}"))
    ret <- tryCatch(optimizeALS(liger_data, k=k.suggest, thresh = thresh, lambda=lambda,nrep = 3), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured in liger: {msg}"))
      error(logger, str_glue("error occured in liger {msg}"))
      structure(msg, class = "try-error")
    })
  }
  liger_data <- ret
  liger_data <- quantileAlignSNF(liger_data,resolution = resolution) #SNF clustering and quantile alignment
  names(liger_data@clusters) <- as.integer.Array( names(liger_data@clusters))
  str_names <- 1:dim(data)[2]
  as.integer(unname(liger_data@clusters[str_names]))
}


### CIDR
cluster.cidr <- function(data,exp_config) {
  require(cidr)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.cidr.batch_free")) methods.config.cidr else methods.config.cidr.batch_free
  stopifnot(is(data,"SingleCellExperiment"))
  
  sData <- scDataConstructor(as.matrix(counts(data)))
  min2 <- if(purrr::is_null(m_config$min2)) 8 else m_config$min2
  
  ret <- tryCatch(determineDropoutCandidates(sData,min2=min2), error=function(c) {
    msg <- conditionMessage(c)
    print(str_glue("error occured in cidr: {msg}"))
    error(logger, str_glue("error occured in cidr {msg}"))
    structure(msg, class = "try-error")
  })
  
  i = 0
  while(inherits(ret,"try-error")){
    if(i > 10||min2 <= 0) return(NULL)
    i = i+1
    print(str_glue("in cidr:{ret}"))
    min2 = min2 -2
    print(str_glue("new min2: {min2}"))
    ret <- tryCatch(determineDropoutCandidates(sData,min2=min2), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured in cidr: {msg}"))
      error(logger, str_glue("error occured in cidr {msg}"))
      structure(msg, class = "try-error")
    })
  }
  
  sData <- wThreshold(ret)
  sData <- scDissim(sData)
  sData <- scPCA(sData, plotPC=FALSE)
  sData <- nPC(sData)
  if(exp_config$known_cluster_num){
    nCluster <- length(unique(colData(data)$label))
  }else{
    nCluster <- m_config$nCluster
  }
  if(nCluster <= 1) return(NULL)
  n <- m_config$n
  if(purrr::is_null(nCluster))
    ret <- tryCatch(scCluster(sData, n=n), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured in cidr: {msg}"))
      error(logger, str_glue("error occured in cidr {msg}"))
      structure(msg, class = "try-error")
    })
  else
    ret <- tryCatch(scCluster(sData, nCluster=nCluster), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured in cidr: {msg}"))
      error(logger, str_glue("error occured in cidr {msg}"))
      structure(msg, class = "try-error")
    })
  
  i = 0
  K = 10
  while(inherits(ret,"try-error")){
    if(i > 10||K <= 1||exp_config$known_cluster_num) return(NULL)
    i = i+1
    print(str_glue("in cidr:{ret}"))
    K = K -2
    print(str_glue("new k: {K}"))
    ret <- tryCatch(scCluster(sData, nCluster=K), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured in cidr: {msg}"))
      error(logger, str_glue("error occured in cidr {msg}"))
      structure(msg, class = "try-error")
    })
  }
  
  
  return(ret@clusters)
}

### Monocle
cluster.monocle3 <- function(data,exp_config) {
  stopifnot(is(data,"SingleCellExperiment"))
  require(monocle3)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.monocle3.batch_free")) methods.config.monocle3 else methods.config.monocle3.batch_free
  num_dim <- if(purrr::is_null(m_config$num_dim)) 100 else m_config$num_dim
  resolution <- if(purrr::is_null(m_config$resolution)) 1e-5 else m_config$resolution
  cds <- new_cell_data_set(as.matrix(counts(data)),
                           cell_metadata = colData(data),
                           gene_metadata = rowData(data))
  
  cds <- preprocess_cds(cds, num_dim = num_dim)
  cds <- reduce_dimension(cds,reduction_method = 'UMAP')
  cds <- cluster_cells(cds, resolution=resolution)
  as.integer(monocle3::clusters(cds))
  
}

###pcaReduce
cluster.pcaReduce <- function(data,exp_config) {
  stopifnot(is(data,"SingleCellExperiment"))
  require(pcaReduce)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.pcaReduce.batch_free")) methods.config.pcaReduce else methods.config.pcaReduce.batch_free
  counts(data) <- as.matrix(counts(data))
  data <- logNormCounts(data)
  Input <- t(logcounts(data))
  if(exp_config$known_cluster_num){
    q <- length(unique(colData(data)$label))-1
    K <- 1  ###get results from Kth level in the hiearchy clustering
  }else{
    q <- if(purrr::is_null(m_config$q)) 25 else m_config$q
    K <- if(purrr::is_null(m_config$K)) 1  else m_config$K
  }
  Output_S <- PCAreduce(Input, nbt=1, q=q, method='S')

  Output_S[[1]][,K]
}


cluster.raceID3 <- function(data, exp_config){
  stopifnot(is(data,"SingleCellExperiment"))
  require(RaceID)
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.raceID3.batch_free")) methods.config.raceID3 else methods.config.raceID3.batch_free
  mintotal <- if(purrr::is_null(m_config$mintotal)) 1000 else m_config$mintotal
  
  if(exp_config$known_cluster_num){
    cln <- length(unique(colData(data)$label))
  }else{
    cln <- m_config$k
  }
  
  sat <- if(purrr::is_null(cln)) T else F
  sc <- SCseq(as.matrix(counts(data)))
  sc <- filterdata(sc,mintotal=mintotal)
  sc <- compdist(sc,metric="pearson")
  dist_m <- sc@distances
  for(i in seq_along(dist_m)){
    if(is.na(dist_m[i])||is.null(dist_m[i])){
      dist_m[i] = 0
    }
  }
  sc@distances <- dist_m
  ret <- tryCatch(clustexp(sc,cln=cln,sat=sat), error=function(c) {
    msg <- conditionMessage(c)
    print(str_glue("error occured in raceID3: {msg}"))
    error(logger, str_glue("error occured in raceID3 {msg}"))
    structure(msg, class = "try-error")
  })
  
  i = 0
  while(inherits(ret,"try-error")){
    if(i > 10||mintotal <= 0) return(NULL)
    i = i+1
    print(str_glue("in raceID3:{ret}"))
    mintotal = mintotal/2
    print(str_glue("new mintotal: {mintotal}"))
    sc <- SCseq(as.matrix(counts(data)))
    sc <- filterdata(sc,mintotal=mintotal)
    sc <- compdist(sc,metric="pearson")
    ret <- tryCatch(clustexp(sc,cln=cln,sat=sat), error=function(c) {
      msg <- conditionMessage(c)
      print(str_glue("error occured in raceID3: {msg}"))
      error(logger, str_glue("error occured in raceID3 {msg}"))
      structure(msg, class = "try-error")
    })
  }
  if(exp_config$known_cluster_num){
    return(ret@cluster$kpart)
  }else{
    sc <- findoutliers(ret)
    return(as.integer(sc@cpart))
  }
}

cluster.same_clustering <- function(data, exp_config){
  require(SAMEclustering)
  stopifnot(is(data,"SingleCellExperiment"))
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.same_clustering.batch_free")) methods.config.same_clustering else methods.config.same_clustering.batch_free
  percent_dropout <- if(purrr::is_null(m_config$percent_dropout)) 0 else m_config$percent_dropout
  resolution <- if(purrr::is_null(m_config$resolution)) 0.7 else m_config$resolution
  dimensions <- if(purrr::is_null(m_config$dimensions)) 2 else m_config$dimensions
  perplexity <- if(purrr::is_null(m_config$perplexity)) 30 else m_config$perplexity
  mt.cutoff <- if(purrr::is_null(m_config$mt.cutoff)) 0.8 else m_config$mt.cutoff
  cluster.result <- cluster.same_clustering.individual_clustering(inputTags = counts(data), mt_filter = TRUE,
                                          percent_dropout = percent_dropout, SC3 = T, CIDR = T, nPC.cidr = NULL, Seurat =T, nGene_filter = FALSE,
                                          nPC.seurat = NULL, resolution = resolution, tSNE = T, dimensions = dimensions, perplexity = perplexity, SIMLR = F, diverse = F, 
                                          save.results = FALSE,mt.cutoff=mt.cutoff, SEED = 123)
  cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = 123)
  as.integer(cluster.ensemble$BICcluster)
}



cluster.same_clustering.individual_clustering <- function(inputTags, mt_filter = TRUE, mt.pattern = "^MT-", mt.cutoff = 0.1, percent_dropout = 10,
                                                          SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL,
                                                          Seurat = TRUE, nGene_filter = TRUE, low.genes = 200, high.genes = 8000, nPC.seurat = NULL, resolution = 0.7, 
                                                          tSNE = TRUE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                                                          SIMLR = TRUE, diverse = TRUE, save.results = FALSE, SEED = 1){
  
  cluster_number <- NULL
  cluster_results <- NULL
  inputTags = as.matrix(inputTags)
  
  # Filter out cells that have mitochondrial genes percentage over 5%
  if (mt_filter == TRUE){
    mito.genes <- grep(pattern = mt.pattern, x = rownames(x = inputTags), value = TRUE)
    percent.mito <- Matrix::colSums(inputTags[mito.genes, ])/Matrix::colSums(inputTags)
    inputTags <- inputTags[,which(percent.mito <= mt.cutoff)]
  }
  
  ##### SC3
  if(SC3 == TRUE){
    message("Performing SC3 clustering...")
    
    sc3OUTPUT <- cluster.same_clustering.sc3_SAME(inputTags = inputTags, gene_filter = gene_filter, svm_num_cells = svm_num_cells, save.results = save.results, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))
  }
  
  
  ##### CIDR
  if(CIDR == TRUE){
    message("Performing CIDR clustering...")
    
    cidrOUTPUT <- cidr_SAME(inputTags = inputTags, percent_dropout = percent_dropout, nPC.cidr = nPC.cidr, save.results = save.results, SEED = SEED)
    
    if(is.null(nPC.cidr)) {
      nPC.cidr <- cidrOUTPUT@nPC
    }
    
    cluster_results <- rbind(cluster_results, matrix(c(cidrOUTPUT@clusters), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number,  cidrOUTPUT@nCluster)
  }
  
  
  ##### Seurat
  if (Seurat == TRUE){
    message("Performing Seurat clustering...")
    
    if(is.null(nPC.seurat)) {
      nPC.seurat <- nPC.cidr
    }
    
    seurat_output <- seurat_SAME(inputTags = inputTags, nGene_filter = nGene_filter, low.genes = low.genes, high.genes = high.genes, 
                                 nPC.seurat = nPC.seurat, resolution = resolution, save.results = save.results, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(seurat_output), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(!is.na(seurat_output)))
  }
  
  
  ##### tSNE+kmeans
  if(tSNE == TRUE){
    message("Performing tSNE + k-means clustering...")
    
    ### Dimensionality reduction by Rtsne
    if(length(inputTags[1,]) < tsne_min_cells) {
      perplexity = tsne_min_perplexity
    }
    
    tsne_kmeansOUTPUT <- tSNE_kmeans_SAME(inputTags = inputTags, percent_dropout = percent_dropout, dimensions = dimensions, perplexity = perplexity, 
                                          k.min = 2, k.max = max(cluster_number), var_genes = var_genes, save.results = save.results, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(as.numeric(tsne_kmeansOUTPUT$cluster)))
  }
  
  ##### SIMLR
  if(SIMLR == TRUE){
    message("Performing SIMLR clustering...")
    
    simlrOUTPUT <- SIMLR_SAME(inputTags = inputTags, percent_dropout = percent_dropout, k.min = 2, k.max = max(cluster_number), save.results = save.results, SEED = SEED)
    cluster_results <- rbind(cluster_results, simlrOUTPUT$y$cluster)
  }
  
  ##### Individual method selection
  if (dim(cluster_results)[1] == 5 && diverse == TRUE){
    message("Selecting clusteirng methods for ensemble...")
    
    rownames(cluster_results) <- c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")
    
    ARI=matrix(0,5,5)
    rownames(ARI) <- c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")
    colnames(ARI) <- c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")
    
    for(i in c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")){
      for(j in c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")){
        ARI[i,j] <- adjustedRandIndex(unlist(cluster_results[i,]), unlist(cluster_results[j,]))
      }
    }
    m1 <- which.min(apply(ARI,1,var))
    cluster_results <- cluster_results[-m1,]
  }
  
  return(cluster_results)
}


cluster.same_clustering.sc3_SAME <- function(inputTags, gene_filter, svm_num_cells, save.results, SEED){
  exp_cell_exprs <- NULL
  sc3OUTPUT <- NULL
  n_cores <- detectCores()-2
  # cell expression
  ### For count data, it would be normalized by the total cound number and then log2 transformed
  exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
  normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
  logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
  
  rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
  exp_cell_exprs <- exp_cell_exprs[!duplicated(rowData(exp_cell_exprs)$feature_symbol), ]
  
  ### Estimating optimal number of clustering
  exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
  optimal_K <- metadata(exp_cell_exprs)$sc3$k_estimation
  
  ### Clustering by SC3 at the optimal K
  if (ncol(inputTags) < svm_num_cells){
    #print(optimal_K)
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter, n_cores = n_cores, rand_seed = SEED)
  } else if (ncol(inputTags) >= svm_num_cells){
    ### Runing SVM
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter,
                          svm_max = svm_num_cells, svm_num_cells = svm_num_cells, n_cores = n_cores, rand_seed = SEED)
    exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
  }
  
  if(save.results == TRUE){
    save(exp_cell_exprs, file = "sc3OUTPUT.Rdata")
  }
  
  ### Exporting SC3 results
  p_Data <- colData(exp_cell_exprs)
  col_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
  sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
  return(sc3OUTPUT)
}





cluster.sharp <- function(data,exp_config){
  require(SHARP)
  stopifnot(is(data,"SingleCellExperiment"))
  m_config <- if(purrr::is_null(exp_config$batch_free) || !exp_config$batch_free || !exists("methods.config.sharp.batch_free")) methods.config.sharp else methods.config.sharp.batch_free
  rN = 2103;
  partition.ncells = if(purrr::is_null(m_config$partition.ncells)) 2000 else m_config$partition.ncells
  ensize.K = m_config$ensize.K
  print(str_glue("partition.ncells = {partition.ncells }"))
  logflag = if(purrr::is_null(m_config$logflag)) T else m_config$logflag
  if(exp_config$known_cluster_num){
    N.cluster <- length(unique(colData(data)$label))
  }else{
    N.cluster <- m_config$N.cluster
  }
  ret <- tryCatch(SHARP(counts(data), prep = TRUE, rN.seed=rN, 
                        forview = FALSE, N.cluster = N.cluster,
                        logflag=logflag, partition.ncells = partition.ncells,ensize.K = ensize.K), 
                  error=function(c) {
    msg <- conditionMessage(c)
    print(str_glue("error occured in sharp: {msg}"))
    error(logger, str_glue("error occured in sharp {msg}"))
    structure(msg, class = "try-error")
  })
  
  i = 1
  while(inherits(ret,"try-error")){
    if(i > 5 || ensize.K <= 1) return(NULL)
    i = i+1
    print(str_glue("in sharp:{ret}"))
    ensize.K = ceiling(ensize.K/2)
    print(str_glue("new ensize.K: {ensize.K}"))
    ret <- tryCatch(SHARP(counts(data), prep = TRUE, rN.seed=rN, 
                          forview = FALSE, N.cluster = N.cluster,
                          logflag=logflag, partition.ncells = partition.ncells,ensize.K = ensize.K), 
                    error=function(c) {
                      msg <- conditionMessage(c)
                      print(str_glue("error occured in sharp: {msg}"))
                      error(logger, str_glue("error occured in sharp {msg}"))
                      structure(msg, class = "try-error")
                    })
  }
  
  
  
  as.integer(ret$pred_clusters)
}


