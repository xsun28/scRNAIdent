experiments.data <- list(simple_accuracy="PBMC_AllCells_withLabels.RDS", 
                                 cell_number="PBMC_AllCells_withLabels.RDS", 
                                 sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                                 celltype_structure="GSE96583_8_Stim_Pats.RDS",
                                 batch_effects_no_free=list(muraro="Muraro_pancreas_clean.RDS",seger="Segerstolpe_pancreas_clean.RDS",xin="Xin_pancreas_clean.RDS"),
                                 # batch_effects_no_free=list("PBMC_AllCells_withLabels.RDS","GSE96583_8_Ctrl_Pats.RDS"),
                                 inter_diseases = list("GSE96583_8_Ctrl_Pats.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_batch1_3_samples.RDS"),
                                 celltype_complexity = list(),
                                 inter_species = list(),
                                 random_noise = list(),
                                 inter_protocol = list("cellbench_10x","cellbench_CELseq2","cellbench_Dropseq")
                                 )

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy="PBMC_AllCells_withLabels.RDS", 
                     cell_number="PBMC_AllCells_withLabels.RDS", 
                     sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                     celltype_structure="GSE96583_8_Stim_Pats.RDS"),
                     # inter_diseases="GSE96583_8_Ctrl_Pats.RDS"),
  test_dataset=list(simple_accuracy="PBMC_AllCells_withLabels.RDS", 
                    cell_number="PBMC_AllCells_withLabels.RDS", 
                    sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                    celltype_structure="GSE96583_8_Stim_Pats.RDS")
                    # inter_diseases="GSE96583_8_Stim_Pats.RDS")
  )
  
##for batch effects removed, scmap and singlecellnet doesn't work
experiments.methods <- list(
  simple_accuracy=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")), 
  cell_number=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  sequencing_depth=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  celltype_structure=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  batch_effects=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign"),
                     cluster_batch_free=c('seurat',"sc3",'tscan','liger'), assign_batch_free=c('chetah','garnett'), marker_gene_assign_batch_free=c("cellassign")
                     ),
  inter_diseases = list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  celltype_complexity = list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  inter_species = list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  random_noise = list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")), 
  inter_protocol = list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign"))
)

experiments.parameters <- list(
  simple_accuracy=list(train_sample_pctg=0.1,train_sample_num=300,test_sample_pctg=0.1,test_sample_num=300, cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'),
                       marker_gene_file=NULL),
  cell_number=list(sample_pctg=c(0.05,0.1,0.15,0.2),sample_num=c(100,200,300,400),train_sample_pctg=0.1,train_sample_num=400, test_sample_pctg=0.1,test_sample_num=300, train_sampling=F, test_sampling=T,
                   cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL),
  sequencing_depth=list(quantile=list(low=0.2,high=0.8),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                        marker_gene_file=NULL),
  celltype_structure=list(train_sample_pctg=0.1,test_sample_pctg=0.1,train_sample_num=300,test_sample_num=300,cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                          structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  batch_effects=list(train_sample_pctg=0.1,test_sample_pctg=0.1,train_sample_num=500,test_sample_num=500,cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                     marker_gene_file=NULL),
  inter_diseases = list(train_sample_pctg=0.1,test_sample_pctg=0.1,train_sample_num=500,test_sample_num=500,cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL),
  celltype_complexity = list(),
  inter_species = list(),
  random_noise = list(),
  inter_protocol = list()
)
