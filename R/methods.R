source("R/methods_config.R")
source("R/config.R")
source("R/marker_gene.R")
run_assign_methods <- function(method,train_data, test_data,exp_config){
  switch(method,
         scmap_cluster = assign.scmap_cluster(train_data, test_data),
         scmap_cell = assign.scmap_cell(train_data, test_data),
         chetah = assign.chetah(train_data, test_data),
         cellassign = assign.cellassign(train_data,exp_config),
         stop("No such assigning method")
  )
  
}

run_cluster_methods <- function(method,data){
  switch(method,
         sc3 = cluster.sc3(data),
         seurat = cluster.seurat(data),
         cidr = cluster.cidr(data),
         tscan = cluster.tscan(data),
         stop("No such cluster method")
  )
}


####assigning using scmap cluster
assign.scmap_cluster <- function(train_data, test_data){
  require(scmap)
  require(scater)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  train_data <- logNormCounts(train_data)
  rowData(train_data)$feature_symbol <- rownames(train_data)
  train_data <- selectFeatures(train_data,n_features=methods.config.scmap[['nfeatures']]) %>%
    indexCluster(cluster_col = "label")
  test_data <- logNormCounts(test_data)
  rowData(test_data)$feature_symbol <- rownames(test_data)
  results <- scmapCluster(projection=test_data, index_list=list(metadata(train_data)$scmap_cluster_index))
  results$combined_labs
}

####assigning using scmap cell
assign.scmap_cell <- function(train_data, test_data){
  require(scmap)
  require(scater)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  train_data <- logNormCounts(train_data)
  rowData(train_data)$feature_symbol <- rownames(train_data)
  if(purrr::is_null(methods.config.scmap[['seed']]))
    set.seed(1)
  else
    set.seed(methods.config.scmap[['seed']])
  train_data <- selectFeatures(train_data,n_features=methods.config.scmap[['nfeatures']]) %>%
    indexCell()
  test_data <- logNormCounts(test_data)
  rowData(test_data)$feature_symbol <- rownames(test_data)
  cell_results <- scmapCell(projection=test_data, index_list = list(metadata(train_data)$scmap_cell_index))
  clusters_results <- scmapCell2Cluster(
    cell_results, 
    list(
      as.character(colData(train_data)$label)
    ),threshold = methods.config.scmap[['threshold']]
  )
  clusters_results$combined_labs
}

###assigning using chetah
assign.chetah <- function(train_data, test_data){
    require(SingleCellExperiment)
    require(CHETAH)
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

###assigning using cellassign
assign.cellassign <- function(data,exp_config){
  require(tensorflow)
  require(cellassign)
  require(scran)
  m_config <- methods.config.cellassign
  marker_gene_file <- exp_config$marker_gene_file
  study <- metadata(data)$study[[1]]
  if(purrr::is_null(marker_gene_file)){
    if(!file.exists(str_glue("{marker_home}/{study}_markergene_{m_config$marker_gene_method}.RDS")))
      markers_mat <- generate_marker_genes(m_config$marker_gene_method,data,str_glue("{marker_home}/{study}_markergene"))
    else
      markers_mat <- read_rds(str_glue("{marker_home}/{study}_markergene_{m_config$marker_gene_method}.RDS"))
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

  counts(data) <- as.matrix(counts(data))
  data$celltypes <- colData(data)$label
  s <- computeSumFactors(data) %>% 
    sizeFactors()
  
  fit <- cellassign(exprs_obj = data[na.omit(matchidx),], 
                    marker_gene_info = markers_mat, 
                    s = s, 
                    learning_rate = m_config$learning_rate, 
                    shrinkage = m_config$shrinkage,
                    verbose = FALSE)
  fit$cell_type
}

### clustering using Seurat
cluster.seurat <- function(data) {
  require(Seurat)
  stopifnot(is(data,"SingleCellExperiment"))
  cnts <- counts(data)
  m_config <- methods.config.seurat
  ## PCA dimention redcution
  seuset <- CreateSeuratObject(cnts, project='simple_accuracy')
  seuset <- NormalizeData(object = seuset)
  seuset <- FindVariableFeatures(object = seuset,nfeatures = m_config[['nfeatures']])
  seuset <- ScaleData(object = seuset)
  seuset <- RunPCA(object = seuset)
  seuset <- FindNeighbors(object = seuset,dims = 1:m_config[['pc_dims']])
  seuset <- FindClusters(object = seuset, resolution = m_config[['resolution']])
  unname(seuset$seurat_clusters)
}

### TSCAN
cluster.tscan <- function(data) {
  require(TSCAN)
  stopifnot(is(data,"SingleCellExperiment"))
  cnts <- counts(data)
  m_config <- methods.config.tscan
  procdata <- preprocess(as.matrix(cnts),cvcutoff=m_config[['cvcutoff']])
  if(!purrr::is_null(m_config[['k']]) )
    lpsmclust <- exprmclust(procdata, clusternum=m_config[['k']])
  else
    lpsmclust <- exprmclust(procdata)
  unname(lpsmclust$clusterid)
}

### SC3
cluster.sc3 <- function(data) {
  require(SC3)
  stopifnot(is(data,"SingleCellExperiment"))
  counts(data) <- as.matrix(counts(data))
  data <- logNormCounts(data)
  m_config <- methods.config.sc3
  k <- m_config$k
  rowData(data)$feature_symbol <- rownames(data)
  data <- data[!duplicated(rowData(data)$feature_symbol),]
  data <- sc3_prepare(data)
  if(purrr::is_null(k)){
    data <- sc3_estimate_k(data)## estimate number of clusters
    k <- metadata(sce)$sc3$k_estimation
  }
  data <- sc3_calc_dists(data)
  data <- sc3_calc_transfs(data)
  data <- sc3_kmeans(data, ks = k)
  data <- sc3_calc_consens(data)
  colTb = colData(data)[,1]
}











### Monocle
cluster.monocle <- function(counts, K) {
  require(monocle)
  dat <- newCellDataSet(counts)
  dat <- estimateSizeFactors(dat)
  dat <- estimateDispersions(dat)
  
  ## pick marker genes for cell clustering
  disp_table <- dispersionTable(dat)
  ix <- order(disp_table[,"mean_expression"], decreasing=TRUE)
  unsup_clustering_genes <- disp_table[ix[50:1000], "gene_id"]
  ## unsup_clustering_genes <- subset(disp_table, mean_expression >= 1)
  dat <- setOrderingFilter(dat, unsup_clustering_genes)
  ## the follwoing step can be slow. Need to keep marker genes number low.
  dat <- reduceDimension(dat, reduction_method="tSNE")
  
  ## clustering
  if( !missing(K) )
    dat <- clusterCells(dat, num_clusters = K)
  else
    dat <- clusterCells(dat)
  
  pData(dat)$Cluster
}

### CIDR
cluster.cidr <- function(counts, K) {
  require(cidr)
  sData <- scDataConstructor(counts)
  sData <- determineDropoutCandidates(sData)
  sData <- wThreshold(sData)
  sData <- scDissim(sData)
  sData <- scPCA(sData, plotPC=FALSE)
  sData <- nPC(sData)
  if(missing(K))
    sData <- scCluster(sData)
  else
    sData <- scCluster(sData, nCluster=K)
  
  return(sData@clusters)
}


