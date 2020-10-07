source('R/config.R')

experiments.cluster.data <- list(simple_accuracy=datasets$PBMC, cell_number=datasets$PBMC, sequencing_depth=datasets$PBMC,
                         cell_types=datasets$PBMC, batch_effects=datasets$pancreas)

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy=datasets$PBMC, cell_number=datasets$PBMC, sequencing_depth=datasets$PBMC,
                     cell_types=datasets$PBMC, batch_effects=datasets$pancreas),
  test_dataset=list(simple_accuracy=datasets$PBMC, cell_number=datasets$PBMC, sequencing_depth=datasets$PBMC,
                    cell_types=datasets$PBMC, batch_effects=datasets$pancreas)
  )
  
experiments.methods <- list(
  simple_accuracy=list(cluster=c('seurat','tscan'),assign=c('scmap','chetah')), 
  cell_number=list(cluster=c('sc3','seurat','tscan','cidr'),assign=c('scmap','chetah','garnet','cellassign')),
  sequencing_depth=list(cluster=c('sc3','suerat','tscan','cidr'),assign=c('scmap','chetah','garnet','cellassign')),
  cell_types=list(cluster=c('sc3','suerat','tscan','cidr'),assign=c('scmap','chetah','garnet','cellassign')),
  batch_effects=list(cluster=c('sc3','suerat','tscan','cidr'),assign=c('scmap','chetah','garnet','cellassign'))
)

experiments.parameters <- list(
  simple_accuracy=list(sample_num=700, cv=TRUE, cv_fold=5, metrics=c('ARI','AMI','FMI')),
  cell_number=list(sample_num=c(100,200,400,700),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI')),
  sequencing_depth=list(),
  cell_types=list(),
  batch_effects=list()
)