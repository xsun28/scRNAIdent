source('R/config.R')
source('R/dataset_config.R')
experiments.data <- list(simple_accuracy="PBMC", cell_number="PBMC", sequencing_depth="PBMC",
                         celltype_structure="PBMC",batch_effects="pancreas")
  
experiments.cluster.data <- list(simple_accuracy=datasets[[experiments.data$simple_accuracy]], 
                                 cell_number=datasets[[experiments.data$cell_number]], 
                                 sequencing_depth=datasets[[experiments.data$sequencing_depth]],
                                 celltype_structure=datasets[[experiments.data$celltype_structure]],
                                 batch_effects_no_free=datasets[[experiments.data$batch_effects]])

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy=datasets[[experiments.data$simple_accuracy]], 
                     cell_number=datasets[[experiments.data$cell_number]], 
                     sequencing_depth=datasets[[experiments.data$sequencing_depth]],
                     celltype_structure=datasets[[experiments.data$celltype_structure]]),
  
  test_dataset=list(simple_accuracy=datasets[[experiments.data$simple_accuracy]], 
                    cell_number=datasets[[experiments.data$cell_number]], 
                    sequencing_depth=datasets[[experiments.data$sequencing_depth]],
                    celltype_structure=datasets[[experiments.data$celltype_structure]])
  )
  
##for batch effects removed, scmap and singlecellnet doesn't work
experiments.methods <- list(
  simple_accuracy=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")), 
  cell_number=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  sequencing_depth=list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign")),
  celltype_structure=list(cluster=c('sc3','seurat','tscan','liger'),assign=c('scmap','chetah','garnet','cellassign','singlecellnet'),marker_gene_assign=c("cellassign")),
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
  celltype_structure=list(),
  batch_effects=list(sample_num=NA,cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                     marker_gene_file=NULL)
)