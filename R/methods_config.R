methods.config.scmap <- list(nfeatures=500,threshold=0.5,seed=1)
methods.config.seurat <- list(nfeatures=2000,pc_dims=10,resolution=0.5)
methods.config.tscan <- list(cvcutoff=0.01,k=8)
methods.config.sc3 <- list(nfeatures=500,k=8)
methods.config.cellassign <- list(learning_rate=1e-2,shrinkage=TRUE,marker_gene_method='seurat')
methods.config.liger <- list(k.suggest=25,lambda=NULL,resolution=NULL,thresh=NULL)
methods.config.singlecellnet <- list(cross_species=FALSE,common_gene_file=NULL,ncells=50,nRand=70,nTrees=1000,nTopGenes=10,nTopGenePairs=25)
methods.config.garnett <- list(PBMC=list(marker_file_path="Garnett_PBMC_marker_ENSEMBL.txt",gene_name_type="ENSEMBL",pretrained_classifier=NULL,marker_gene_method='seurat'),
                               pancreas=list(marker_file_path="Garnett_pancreas_marker_SYMBOL.txt",gene_name_type="SYMBOL",pretrained_classifier=NULL,marker_gene_method='seurat'))
