source('R/config.R')

experiments.cluster.data <- list(simple_accuracy=datasets$PBMC, cell_number=datasets$PBMC, sequencing_depth=datasets$PBMC,
                         cell_types=datasets$PBMC,batch_effects_no_free=datasets$pancreas, 
                         batch_effects_free=batch_effects_free_datasets$pancreas)

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy=datasets$PBMC, cell_number=datasets$PBMC, sequencing_depth=datasets$PBMC,
                     cell_types=datasets$PBMC),
  test_dataset=list(simple_accuracy=datasets$PBMC, cell_number=datasets$PBMC, sequencing_depth=datasets$PBMC,
                    cell_types=datasets$PBMC)
  )
  
experiments.methods <- list(
  simple_accuracy=list(cluster=c('seurat','tscan'),assign=c('scmap','chetah')), 
  cell_number=list(cluster=c('seurat','tscan'),assign=c('scmap','chetah')),
  sequencing_depth=list(cluster=c('seurat','tscan'),assign=c('scmap','chetah')),
  cell_types=list(cluster=c('sc3','seurat','tscan'),assign=c('scmap','chetah','garnet','cellassign')),
  batch_effects=list(cluster=c('seurat','tscan'),assign=c('scmap','chetah'))
)

experiments.parameters <- list(
  simple_accuracy=list(sample_num=700, cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI')),
  cell_number=list(sample_num=c(100,200,400,700),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI')),
  sequencing_depth=list(quantile=list(low=0.2,high=0.8),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI')),
  cell_types=list(),
  batch_effects=list(sample_num=NA,cv=FALSE,metrics=c('ARI','AMI','FMI'))
)