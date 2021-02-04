experiments.cluster.data <- list(simple_accuracy="GSE96583_8_Stim_Pats.RDS", 
                                 cell_number="GSE96583_8_Stim_Pats.RDS", 
                                 sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                                 celltype_structure="GSE96583_8_Stim_Pats.RDS",
                                 batch_effects_no_free=list(muraro="Muraro_pancreas_clean.RDS",seger="Segerstolpe_pancreas_clean.RDS",xin="Xin_pancreas_clean.RDS")
                                 # batch_effects_no_free=list("PBMC_AllCells_withLabels.RDS","GSE96583_8_Ctrl_Pats.RDS")
                                 )

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy="GSE96583_8_Stim_Pats.RDS", 
                     cell_number="GSE96583_8_Stim_Pats.RDS", 
                     sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                     celltype_structure="GSE96583_8_Stim_Pats.RDS"),
  
  test_dataset=list(simple_accuracy="GSE96583_8_Stim_Pats.RDS", 
                    cell_number="GSE96583_8_Stim_Pats.RDS", 
                    sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                    celltype_structure="GSE96583_8_Stim_Pats.RDS")
  )
  
##for batch effects removed, scmap and singlecellnet doesn't work
experiments.methods <- list(
  simple_accuracy=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")), 
  cell_number=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  sequencing_depth=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  celltype_structure=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  batch_effects=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign"),
                     cluster_batch_free=c('seurat',"sc3",'tscan','liger'), assign_batch_free=c('chetah','garnett'), marker_gene_assign_batch_free=c("cellassign")
                     )
)

experiments.parameters <- list(
  simple_accuracy=list(sample_num=300, cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'),
                       marker_gene_file=NULL),
  cell_number=list(sample_num=c(100,200,400),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                   marker_gene_file=NULL),
  sequencing_depth=list(quantile=list(low=0.2,high=0.8),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                        marker_gene_file=NULL),
  celltype_structure=list(sample_num=300,cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                          structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  batch_effects=list(sample_num=NA,cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                     marker_gene_file=NULL)
)
