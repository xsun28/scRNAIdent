source("R/methods_config.R")

run_assign_methods <- function(method,train_data, test_data){
  switch(method,
         scmap = assign.scmap(train_data, test_data),
         chetah = assign.chetah(train_data, test_data),
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


####assigning using scmap
assign.scmap <- function(train_data, test_data){
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
assign.cellassign <- function(train_data, test_data){
  
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
  if(!is_null(m_config[['k']]) )
    lpsmclust <- exprmclust(procdata, clusternum=m_config[['k']])
  else
    lpsmclust <- exprmclust(procdata)
  unname(lpsmclust$clusterid)
}

### SC3
cluster.sc3 <- function(data, K) {
  require(SC3)
  
}









### clustering using SC3
cluster.sc3 <- function(counts, K) {
  require(SC3)
  sce = SingleCellExperiment(
    assays = list(
      counts = as.matrix(counts),
      logcounts = log2(as.matrix(counts) + 1)
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
  sce = sc3_prepare(sce)
  if( missing(K) ) { ## estimate number of clusters
    sce = sc3_estimate_k(sce)
    K = metadata(sce)$sc3$k_estimation
  }
  
  sce = sc3_calc_dists(sce)
  sce = sc3_calc_transfs(sce)
  sce = sc3_kmeans(sce, ks = K)
  sce = sc3_calc_consens(sce)
  colTb = colData(sce)[,1]
  return(colTb)
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


