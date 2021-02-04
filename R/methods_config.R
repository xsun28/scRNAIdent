methods.config.scmap <- list(nfeatures=500,threshold=0.5,seed=1)
methods.config.seurat <- list(nfeatures=2000,pc_dims=10,resolution=0.5)
methods.config.tscan <- list(cvcutoff=0.01,k=8)
methods.config.sc3 <- list(nfeatures=500,k=8,gene_filter = T)
methods.config.sc3.batch_free <- list(nfeatures=500,k=8,gene_filter = FALSE)
methods.config.cellassign <- list(learning_rate=1e-2,shrinkage=TRUE,marker_gene_method='seurat')
methods.config.liger <- list(suggestK=F,k.suggest=25,lambda=NULL,resolution=NULL,thresh=NULL)
methods.config.singlecellnet <- list(cross_species=FALSE,common_gene_file=NULL,ncells=50,nRand=70,nTrees=1000,nTopGenes=10,nTopGenePairs=25)

methods.config.garnett <- list(PBMC_AllCells_withLabels=list(marker_file_path=str_glue("Garnett_PBMC_AllCells_withLabels_marker_{dataset.properties$PBMC_AllCells_withLabels$gene_name_type}.txt"),
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               GSE96583_batch1_3_samples=list(marker_file_path=str_glue("Garnett_GSE96583_batch1_3_samples_marker_{dataset.properties$GSE96583_batch1_3_samples$gene_name_type}.txt"),
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               GSE96583_8_Stim_Pats=list(marker_file_path=str_glue("Garnett_GSE96583_8_Stim_Pats_marker_{dataset.properties$GSE96583_8_Stim_Pats$gene_name_type}.txt"),
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               GSE96583_8_Ctrl_Pats=list(marker_file_path=str_glue("Garnett_GSE96583_8_Ctrl_Pats_marker_{dataset.properties$GSE96583_8_Ctrl_Pats$gene_name_type}.txt"),
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               ADASD_AD=list(marker_file_path="Garnett_AD_marker_{dataset.properties$ADASD_AD$gene_name_type}.txt",
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               ADASD_autism=list(marker_file_path="Garnett_autism_marker_{dataset.properties$ADASD_autism$gene_name_type}.txt",
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               Muraro_pancreas=list(marker_file_path="Garnett_pancreas_marker_{dataset.properties$Muraro_pancreas$gene_name_type}.txt",
                                                              pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               Segerstolpe_pancreas=list(marker_file_path="Garnett_pancreas_marker_{dataset.properties$Segerstolpe_pancreas$gene_name_type}.txt",
                                                    pretrained_classifier=NULL,marker_gene_method='seurat'),
                               
                               Xin_pancreas=list(marker_file_path="Garnett_pancreas_marker_{dataset.properties$Xin_pancreas$gene_name_type}.txt",
                                                    pretrained_classifier=NULL,marker_gene_method='seurat')
                               )
