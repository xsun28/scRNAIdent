experiment <- "batch_effects"

experiments.data <- list(simple_accuracy="GSE96583_batch1_3_samples.RDS", 
                                 cell_number="midbrain_mouse.RDS", 
                                 sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                                 celltype_structure="GSE96583_8_Stim_Pats.RDS",
                                 batch_effects=list(muraro="Muraro_pancreas_clean.RDS",seger="Segerstolpe_pancreas_clean.RDS",xin="Xin_pancreas_clean.RDS"),
                                 # batch_effects_no_free=list("PBMC_AllCells_withLabels.RDS","GSE96583_8_Ctrl_Pats.RDS"),
                                 inter_diseases = list("GSE96583_8_Ctrl_Pats.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_batch1_3_samples.RDS"),
                                 celltype_complexity = list(),
                                 inter_species = list(),
                                 random_noise = list(),
                                 inter_protocol = list("cellbench_10x","cellbench_CELseq2","cellbench_Dropseq")
                                 )

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy="GSE96583_batch1_3_samples.RDS", 
                     cell_number="midbrain_mouse.RDS", 
                     sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                     celltype_structure="GSE96583_8_Stim_Pats.RDS"),
                     # inter_diseases="GSE96583_8_Ctrl_Pats.RDS"),
  test_dataset=list(simple_accuracy="GSE96583_batch1_3_samples.RDS", 
                    cell_number="midbrain_mouse.RDS", 
                    sequencing_depth="GSE96583_8_Stim_Pats.RDS",
                    celltype_structure="GSE96583_8_Stim_Pats.RDS")
                    # inter_diseases="GSE96583_8_Stim_Pats.RDS")
  )
  
##for batch effects removed, scmap and singlecellnet doesn't work
experiments.methods.base_config <- list(cluster=c('seurat','tscan','sc3','liger'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett'),marker_gene_assign=c("cellassign"))
experiments.methods <- list(
  simple_accuracy=experiments.methods.base_config, 
  cell_number=experiments.methods.base_config,
  sequencing_depth=experiments.methods.base_config,
  celltype_structure=experiments.methods.base_config,
  batch_effects=experiments.methods.base_config, 
                       # list(cluster_batch_free=c('seurat',"sc3",'tscan','liger'), assign_batch_free=c('chetah','garnett'), marker_gene_assign_batch_free=c("cellassign"))),
  inter_diseases = experiments.methods.base_config,
  celltype_complexity = experiments.methods.base_config,
  inter_species = experiments.methods.base_config,
  random_noise = experiments.methods.base_config, 
  inter_protocol = experiments.methods.base_config
)



######simple accuracy config
experiments.parameters.simple_accuracy.base_config <- list(train_sample_pctg=0.8,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL, 
                                                           cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.PBMC <- list(train_sample_pctg=0.2,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL, 
                                                                cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.pancreas <- experiments.parameters.simple_accuracy.base_config

experiments.parameters.simple_accuracy.base_config.ADASD <- list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL, 
                                                                 cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.midbrain <- experiments.parameters.simple_accuracy.base_config

experiments.parameters.simple_accuracy.base_config.cellbench <- experiments.parameters.simple_accuracy.base_config


experiments.parameters.simple_accuracy <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL, 
                                  cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.simple_accuracy.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.simple_accuracy.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.simple_accuracy.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.simple_accuracy.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.simple_accuracy.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.simple_accuracy.base_config.cellbench,
  cellbench_Dropseq.RDS = experiments.parameters.simple_accuracy.base_config.cellbench
)


#####cell number config
experiments.parameters.cell_number.base_config <- list(sample_pctg=c(0.4,0.6,0.8,1.0),sample_num=NULL,train_sample_pctg=0.7,
                                                       train_sample_num=NULL, test_sample_pctg=NULL,test_sample_num=NULL,
                                                       train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                       marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.PBMC <- list(sample_pctg=c(0.05,0.1,0.2,0.3),sample_num=NULL,train_sample_pctg=0.1,
                                                            train_sample_num=NULL, test_sample_pctg=0.1,test_sample_num=NULL,
                                                            train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                            marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.pancreas <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number.base_config.ADASD <- list(sample_pctg=c(0.02,0.03,0.04,0.05),sample_num=NULL,train_sample_pctg=0.06,
                                                             train_sample_num=NULL, test_sample_pctg=0.1,test_sample_num=NULL,
                                                             train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                             marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.midbrain <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number.base_config.cellbench <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number <- list(
  PBMC_AllCells_withLabels.RDS = list(sample_pctg=c(0.03,0.05,0.08,0.1),sample_num=NULL,train_sample_pctg=0.1,
                                      train_sample_num=NULL, test_sample_pctg=0.1,test_sample_num=NULL,
                                      train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                      marker_gene_file=NULL,trained=F),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.cell_number.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.cell_number.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.cell_number.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.cell_number.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.cell_number.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.cell_number.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.cell_number.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.cell_number.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.cell_number.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.cell_number.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.cell_number.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.cell_number.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.cell_number.base_config.cellbench
)

####sequencing depth config
experiments.parameters.sequencing_depth.base_config <- list(quantile=list(low=0.3,high=0.7),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                                                            marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.PBMC <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth.base_config.pancreas <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth.base_config.ADASD <- list(quantile=list(low=0.04,high=0.96),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                                                                  marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.midbrain <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth.base_config.cellbench <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth <- list(
  PBMC_AllCells_withLabels.RDS = list(quantile=list(low=0.05,high=0.95),cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                                      marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.sequencing_depth.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.sequencing_depth.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.sequencing_depth.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.sequencing_depth.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.sequencing_depth.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.sequencing_depth.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.sequencing_depth.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.sequencing_depth.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.sequencing_depth.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.sequencing_depth.base_config.cellbench
)

####celltype_structure config
experiments.parameters.celltype_structure.base_config <- list(train_sample_pctg=0.8,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL,
                                                            cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                                            structure_file="PBMC_celltypes.csv",marker_gene_file=NULL)

experiments.parameters.celltype_structure.base_config.PBMC <- list(train_sample_pctg=0.2,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL,
                                                                   cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                                                   structure_file="PBMC_celltypes.csv",marker_gene_file=NULL)

experiments.parameters.celltype_structure.base_config.pancreas <- experiments.parameters.celltype_structure.base_config

experiments.parameters.celltype_structure.base_config.ADASD <- list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                                                    cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                                                    structure_file="PBMC_celltypes.csv",marker_gene_file=NULL)

experiments.parameters.celltype_structure.base_config.midbrain <- experiments.parameters.celltype_structure.base_config

experiments.parameters.celltype_structure.base_config.cellbench <- experiments.parameters.celltype_structure.base_config

experiments.parameters.celltype_structure <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                      cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                      structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.celltype_structure.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.celltype_structure.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.celltype_structure.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.celltype_structure.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.celltype_structure.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.celltype_structure.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.celltype_structure.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.celltype_structure.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.celltype_structure.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.celltype_structure.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.celltype_structure.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.celltype_structure.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.celltype_structure.base_config.cellbench
)

####batch effects config

experiments.parameters.batch_effects.base_config  <- list(train_sample_pctg=0.8,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL,
                                                          cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                          marker_gene_file=NULL,fixed_train=T,fixed_test=T)


experiments.parameters.batch_effects.base_config.PBMC <- list(train_sample_pctg=0.2,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL,
                                                              cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                              marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.pancreas <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects.base_config.ADASD <- list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                                               cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                               marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.midbrain <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects.base_config.cellbench <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                      cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                      marker_gene_file=NULL,fixed_train=T,fixed_test=T),
  
  GSE96583_batch1_3_samples.RDS = experiments.parameters.batch_effects.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.batch_effects.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.batch_effects.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.batch_effects.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.batch_effects.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.batch_effects.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.batch_effects.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.batch_effects.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.batch_effects.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.batch_effects.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.batch_effects.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.batch_effects.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.batch_effects.base_config.cellbench
)

####inter diseases config
experiments.parameters.inter_diseases.base_config = list(train_sample_pctg=0.8,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL,
                                                         cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)


experiments.parameters.inter_diseases.base_config.PBMC <- list(train_sample_pctg=0.2,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL, 
                                                                   cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.pancreas <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.inter_diseases.base_config.ADASD <- list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                                                    cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.midbrain <- experiments.parameters.inter_diseases.base_config

experiments.parameters.inter_diseases.base_config.cellbench <- experiments.parameters.inter_diseases.base_config

experiments.parameters.inter_diseases <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                      cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.inter_diseases.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.inter_diseases.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.inter_diseases.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.inter_diseases.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.inter_diseases.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.inter_diseases.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.inter_diseases.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.inter_diseases.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.inter_diseases.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.inter_diseases.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.inter_diseases.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.inter_diseases.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.inter_diseases.base_config.cellbench
)

##########
experiments.parameters <- list(
  simple_accuracy=list(train_sample_pctg=experiments.parameters.simple_accuracy[[experiments.assign.data$train_dataset$simple_accuracy]]$train_sample_pctg,
                       train_sample_num=experiments.parameters.simple_accuracy[[experiments.assign.data$train_dataset$simple_accuracy]]$train_sample_num,
                       test_sample_pctg=experiments.parameters.simple_accuracy[[experiments.assign.data$test_dataset$simple_accuracy]]$test_sample_pctg,
                       test_sample_num=experiments.parameters.simple_accuracy[[experiments.assign.data$test_dataset$simple_accuracy]]$test_sample_num, 
                       cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'),
                       marker_gene_file=NULL),
  
  cell_number=list(sample_pctg=experiments.parameters.cell_number[[experiments.data$cell_number]]$sample_pctg,
                   sample_num=experiments.parameters.cell_number[[experiments.data$cell_number]]$sample_num,
                   train_sample_pctg=experiments.parameters.cell_number[[experiments.assign.data$train_dataset$cell_number]]$train_sample_pctg,
                   train_sample_num=experiments.parameters.cell_number[[experiments.assign.data$train_dataset$cell_number]]$train_sample_num,
                   test_sample_pctg=experiments.parameters.cell_number[[experiments.assign.data$test_dataset$cell_number]]$test_sample_pctg,
                   test_sample_num=experiments.parameters.cell_number[[experiments.assign.data$test_dataset$cell_number]]$test_sample_num,
                   train_sampling=experiments.parameters.cell_number[[experiments.assign.data$train_dataset$cell_number]]$train_sampling, 
                   test_sampling=experiments.parameters.cell_number[[experiments.assign.data$test_dataset$cell_number]]$test_sampling, 
                   cv=F,cv_fold=F, metrics=c('ARI','AMI','FMI'), 
                   marker_gene_file=NULL,trained=F),
  
  sequencing_depth=list(quantile=experiments.parameters.sequencing_depth[[experiments.data$sequencing_depth]]$quantile,
                        cv=TRUE,cv_fold=5,metrics=c('ARI','AMI','FMI'),
                        marker_gene_file=NULL),
  
  celltype_structure=list(train_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_pctg,
                          train_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_num,
                          test_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_pctg,
                          test_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_num,
                          cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                          structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  
  batch_effects=list(train_sample_pctg=NULL,train_sample_num=NULL,
                     test_sample_pctg=NULL,test_sample_num=NULL,
                     cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                     marker_gene_file=NULL,fixed_train=T,fixed_test=T),
  
  inter_diseases = list(train_sample_pctg=NULL,train_sample_num=NULL,
                        test_sample_pctg=NULL,test_sample_num=NULL,
                        cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL),
  
  celltype_complexity = list(),
  inter_species = list(),
  random_noise = list(),
  inter_protocol = list()
)
