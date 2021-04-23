experiment <- "cell_number"

# experiments.data <- list(simple_accuracy="PBMC_AllCells_withLabels.RDS", 
#                                  cell_number="ADASD_autism.RDS",
#                                  celltype_number = "midbrain_human.RDS",
#                                  sequencing_depth="ADASD_AD.RDS",####only use PBMC and ADASD datasets
#                                  celltype_structure="GSE96583_8_Stim_Pats.RDS",
#                                  # batch_effects=list(muraro="Muraro_pancreas_clean.RDS",seger="Segerstolpe_pancreas_clean.RDS"),
#                                  batch_effects=list("PBMC_AllCells_withLabels.RDS","GSE96583_8_Ctrl_Pats.RDS"),
#                                  inter_diseases = list("GSE96583_8_Ctrl_Pats.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_batch1_3_samples.RDS"),
#                                  # inter_diseases = list("ADASD_AD.RDS","ADASD_autism.RDS"),
#                                  celltype_complexity = list(),
#                                  inter_species = list(),
#                                  random_noise = list(),
#                                  inter_protocol = list("cellbench_10x","cellbench_CELseq2","cellbench_Dropseq")
#                                  )

experiments.assign.data <- list(
  train_dataset=list(simple_accuracy="PBMC_AllCells_withLabels.RDS",
                     cell_number="GSE96583_8_Ctrl_Pats.RData",
                     celltype_number = "midbrain_human.RDS",
                     sequencing_depth="ADASD_AD.RDS",
                     celltype_structure="GSE96583_8_Stim_Pats.RDS"),
                     # inter_diseases="GSE96583_8_Ctrl_Pats.RDS"),
  test_dataset=list(simple_accuracy="PBMC_AllCells_withLabels.RDS",
                    cell_number="PBMC_AllCells_withLabels.RDS",
                    celltype_number = "midbrain_human.RDS",
                    sequencing_depth="ADASD_AD.RDS",
                    celltype_structure="GSE96583_8_Stim_Pats.RDS")
                    # inter_diseases="GSE96583_8_Stim_Pats.RDS")
  )
  
##for batch effects removed, scmap and singlecellnet doesn't work
experiments.methods.base_config <- list(cluster=c('seurat','tscan','sc3','liger','cidr','monocle3','pcaReduce'),assign=c('scmap_cluster','scmap_cell','chetah','singlecellnet','garnett','singleR'),marker_gene_assign=c("cellassign"))
experiments.methods <- list(
  simple_accuracy=experiments.methods.base_config, 
  cell_number=experiments.methods.base_config,
  celltype_number = experiments.methods.base_config,
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
experiments.parameters.simple_accuracy.base_config <- list(train_sample_pctg=0.7,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL, 
                                                           cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.PBMC <- list(train_sample_pctg=0.08,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL, 
                                                                cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.pancreas <- list(train_sample_pctg=0.5,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL, 
                                                                    cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.ADASD <- list(train_sample_pctg=0.015,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL, 
                                                                 cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.midbrain <- experiments.parameters.simple_accuracy.base_config

experiments.parameters.simple_accuracy.base_config.cellbench <- experiments.parameters.simple_accuracy.base_config


experiments.parameters.simple_accuracy <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.02,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL, 
                                  cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.simple_accuracy.base_config.ADASD,
  ADASD_autism.RDS = list(train_sample_pctg=0.01,train_sample_num=NULL,test_sample_pctg=0.005,test_sample_num=NULL, 
                          cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL),
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

experiments.parameters.cell_number.base_config.PBMC <- list(sample_pctg=c(0.1,0.2,0.3,0.4),sample_num=NULL,train_sample_pctg=0.08,
                                                            train_sample_num=NULL, test_sample_pctg=0.1,test_sample_num=NULL,
                                                            train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                            marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.pancreas <- list(sample_pctg=c(0.3,0.5,0.6,0.8),sample_num=NULL,train_sample_pctg=0.5,
                                                                train_sample_num=NULL, test_sample_pctg=NULL,test_sample_num=NULL,
                                                                train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                                marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.ADASD <- list(sample_pctg=c(0.01,0.015,0.03,0.04),sample_num=NULL,train_sample_pctg=0.015,
                                                             train_sample_num=NULL, test_sample_pctg=0.1,test_sample_num=NULL,
                                                             train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                             marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.midbrain <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number.base_config.cellbench <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number <- list(
  PBMC_AllCells_withLabels.RDS = list(sample_pctg=c(0.03,0.05,0.08,0.1),sample_num=NULL,train_sample_pctg=0.02,
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
  ADASD_autism.RDS = list(sample_pctg=c(0.008,0.01,0.015,0.03),sample_num=NULL,train_sample_pctg=0.01,
                          train_sample_num=NULL, test_sample_pctg=0.1,test_sample_num=NULL,
                          train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                          marker_gene_file=NULL,trained=F),
  midbrain_human.RDS = experiments.parameters.cell_number.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.cell_number.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.cell_number.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.cell_number.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.cell_number.base_config.cellbench
)

#####celltype number config
experiments.parameters.celltype_number.base_config <- list(type_pctg=c(0.2,0.5,0.8,1),train_sample_pctg=0.7,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL, 
                                                           cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.celltype_number.base_config.PBMC <- list(type_pctg=c(0.2,0.5,0.8,1),train_sample_pctg=0.7,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL, 
                                                                cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.celltype_number.base_config.pancreas <- list(type_pctg=c(0.2,0.5,0.8,1),train_sample_pctg=0.7,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL, 
                                                                    cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.celltype_number.base_config.ADASD <- list(type_pctg=c(0.2,0.5,0.8,1),train_sample_pctg=0.7,train_sample_num=NULL,test_sample_pctg=1,test_sample_num=NULL, 
                                                                 cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.celltype_number.base_config.midbrain <- experiments.parameters.celltype_number.base_config

experiments.parameters.celltype_number.base_config.cellbench <- experiments.parameters.celltype_number.base_config

experiments.parameters.celltype_number <- list(
  PBMC_AllCells_withLabels.RDS = experiments.parameters.celltype_number.base_config.PBMC,
  GSE96583_batch1_3_samples.RDS = experiments.parameters.celltype_number.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.celltype_number.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.celltype_number.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.celltype_number.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.celltype_number.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.celltype_number.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.celltype_number.base_config.ADASD,
  ADASD_autism.RDS = experiments.parameters.celltype_number.base_config.ADASD,
  midbrain_human.RDS = experiments.parameters.celltype_number.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.celltype_number.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.celltype_number.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.celltype_number.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.celltype_number.base_config.cellbench
)


####sequencing depth config
experiments.parameters.sequencing_depth.base_config <- list(train_quantiles=list(c(low=0.0, high=0.1),c(low=0.3, high=0.4),c(low=0.6, high=0.7),c(low=0.9, high=1)),
                                                            test_quantiles=list(c(low=0.0, high=0.1),c(low=0.3, high=0.4),c(low=0.6, high=0.7),c(low=0.9, high=1)),
                                                            cv=F, cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,train_sample_pctg=0.1, test_sample_pctg=0.1,
                                                            marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.PBMC <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth.base_config.pancreas <- list(train_quantiles=list(c(low=0.0, high=0.25),c(low=0.25, high=0.5),c(low=0.5, high=0.75),c(low=0.75, high=1)),
                                                                     test_quantiles=list(c(low=0.0, high=0.25),c(low=0.25, high=0.5),c(low=0.5, high=0.75),c(low=0.75, high=1)),
                                                                     cv=F,
                                                                     cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,train_sample_pctg=0.4, test_sample_pctg=0.4,
                                                                     marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.ADASD <- list(train_quantiles=list(c(low=0., high=0.02),c(low=0.24, high=0.26),c(low=0.74, high=0.76),c(low=0.98, high=1)),
                                                                  test_quantiles=list(c(low=0., high=0.02),c(low=0.24, high=0.26),c(low=0.74, high=0.76),c(low=0.98, high=1)),
                                                                  cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,train_sample_pctg=0.02, test_sample_pctg=0.02,
                                                                  marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.midbrain <- list(train_quantiles=list(c(low=0.0, high=0.25),c(low=0.25, high=0.5),c(low=0.5, high=0.75),c(low=0.75, high=1)),
                                                                     test_quantiles=list(c(low=0.0, high=0.25),c(low=0.25, high=0.5),c(low=0.5, high=0.75),c(low=0.75, high=1)),
                                                                     cv=F,
                                                                     cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,train_sample_pctg=0.7, test_sample_pctg=0.7,
                                                                     marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.cellbench <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth <- list(
  PBMC_AllCells_withLabels.RDS = experiments.parameters.sequencing_depth.base_config.ADASD,
  GSE96583_batch1_3_samples.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  Muraro_pancreas_clean.RDS = list(train_quantiles=list(c(low=0.0, high=0.25),c(low=0.25, high=0.5),c(low=0.5, high=0.75),c(low=0.75, high=1)),
                                   test_quantiles=list(c(low=0.0, high=0.25),c(low=0.25, high=0.5),c(low=0.5, high=0.75),c(low=0.75, high=1)),
                                   cv=F,
                                   cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,train_sample_pctg=0.6, test_sample_pctg=0.6,
                                   marker_gene_file=NULL),
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


experiments.parameters.batch_effects.base_config.PBMC <- list(train_sample_pctg=0.07,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL,
                                                              cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                              marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.pancreas <- list(train_sample_pctg=0.5,train_sample_num=NULL,test_sample_pctg=0.5,test_sample_num=NULL,
                                                                  cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                                  marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.ADASD <- list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.05,test_sample_num=NULL,
                                                               cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                               marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.midbrain <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects.base_config.cellbench <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.05,train_sample_num=NULL,test_sample_pctg=0.017,test_sample_num=NULL,
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


experiments.parameters.inter_diseases.base_config.PBMC <- list(train_sample_pctg=0.08,train_sample_num=NULL,test_sample_pctg=0.08,test_sample_num=NULL, 
                                                               cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.pancreas <- list(train_sample_pctg=0.7,train_sample_num=NULL,test_sample_pctg=0.2,test_sample_num=NULL, 
                                                                   cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.ADASD <- list(train_sample_pctg=0.016,train_sample_num=NULL,test_sample_pctg=0.016,test_sample_num=NULL,
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
  ADASD_AD.RDS =  experiments.parameters.inter_diseases.base_config.ADASD,
  ADASD_autism.RDS = list(train_sample_pctg=0.011,train_sample_num=NULL,test_sample_pctg=0.011,test_sample_num=NULL,
                          cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL),
  midbrain_human.RDS = experiments.parameters.inter_diseases.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.inter_diseases.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.inter_diseases.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.inter_diseases.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.inter_diseases.base_config.cellbench
)

##########
experiments.parameters <- list(
  simple_accuracy=list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'),batch_free=F,
                       marker_gene_file=NULL,use_intra_dataset=T,intra_dataset=dataset.datasets,
                       use_inter_dataset=F,inter_dataset=NULL),
  
  cell_number=list( cv=F,cv_fold=F, metrics=c('ARI','AMI','FMI'), batch_free=F,
                   marker_gene_file=NULL,trained=F,use_intra_dataset=F,intra_dataset=list(),
                   use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                          dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC5,
                                                          dataset.interdatasets$PBMC6,dataset.interdatasets$PBMC8,
                                                          dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                          dataset.interdatasets$pancreas5,dataset.interdatasets$ADASD2,
                                                          dataset.interdatasets$midbrain2
                                                          )),

  
  sequencing_depth=list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
                        marker_gene_file=NULL,batch_free=F,use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                               dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC5,
                                                               dataset.interdatasets$PBMC6,dataset.interdatasets$PBMC8,
                                                               dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                               dataset.interdatasets$pancreas5,dataset.interdatasets$ADASD2,
                                                               dataset.interdatasets$midbrain2)),
  
  # celltype_structure=list(train_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_pctg,
  #                         train_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_num,
  #                         test_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_pctg,
  #                         test_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_num,
  #                         cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),batch_free=F,
  #                         structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  
  batch_effects=list(train_sample_pctg=NULL,train_sample_num=NULL,
                     test_sample_pctg=NULL,test_sample_num=NULL,
                     cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                     marker_gene_file=NULL,fixed_train=T,fixed_test=T,batch_correct_method="MNN",use_intra_dataset=T,intra_dataset=list(),
                     use_inter_dataset=F,inter_dataset=list()),
  
  inter_diseases = list(train_sample_pctg=NULL,train_sample_num=NULL,
                        test_sample_pctg=NULL,test_sample_num=NULL,batch_free=F,
                        cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL,use_intra_dataset=T,intra_dataset=list(),
                        use_inter_dataset=F,inter_dataset=list()),
  
  celltype_complexity = list(),
  inter_species = list(),
  random_noise = list(),
  inter_protocol = list()
)


experiments.config.update <- function(experiment, train_dataset, test_dataset=NULL, exp_config){
  switch(experiment,
         simple_accuracy=experiments.config.update.simple_accuracy(train_dataset, test_dataset, exp_config),
         cell_number = experiments.config.update.cell_number(train_dataset, test_dataset, exp_config),
         celltype_number = experiments.config.update.celltype_number(train_dataset, test_dataset, exp_config),
         sequencing_depth = experiments.config.update.sequencing_depth(train_dataset, test_dataset, exp_config),
         celltype_structure = experiments.config.update.celltype_structure(train_dataset, test_dataset, exp_config),
         batch_effects = experiments.config.update.batch_effects(train_dataset, test_dataset, exp_config),
         inter_diseases = experiments.config.update.inter_diseases(train_dataset, test_dataset, exp_config),
         celltype_complexity = experiments.config.update.celltype_complexity(train_dataset, test_dataset, exp_config),
         inter_species = experiments.config.update.inter_species(train_dataset, test_dataset, exp_config),
         random_noise = experiments.config.update.random_noise(train_dataset, test_dataset, exp_config),
         inter_protocol = experiments.config.update.inter_protocol(train_dataset, test_dataset, exp_config),
         stop("Unkown experiments")
         )
}


experiments.config.update.simple_accuracy <-function(train_dataset, test_dataset=NULL, exp_config){
  # train_test_sample_pctgs <- utils.calculate_sampling_pctg(train_dataset, test_dataset,exp_config)
  # train_sample_pctg <- train_test_pctgs[[1]]
  # test_sample_pctg <- train_test_pctgs[[2]]
  exp_config$train_sample_pctg <- experiments.parameters.simple_accuracy[[train_dataset]]$train_sample_pctg
  exp_config$train_sample_num <- experiments.parameters.simple_accuracy[[train_dataset]]$train_sample_num
  exp_config$test_sample_pctg <- experiments.parameters.simple_accuracy[[test_dataset]]$test_sample_pctg
  exp_config$test_sample_num <- experiments.parameters.simple_accuracy[[test_dataset]]$test_sample_num
  exp_config
}

experiments.config.update.cell_number <-function(train_dataset, test_dataset=NULL, exp_config){
  if(purrr::is_null(test_dataset)){
    print("cell number is intra_dataset, using train dataset properties...")
    exp_config$sample_pctg <- experiments.parameters.cell_number[[train_dataset]]$sample_pctg
    exp_config$sample_num <- experiments.parameters.cell_number[[train_dataset]]$sample_num
  }
  if(exp_config$train_sampling&&!exp_config$test_sampling){
    exp_config$fixed_test <- T
  }else if(!exp_config$train_sampling&&exp_config$test_sampling){
    exp_config$fixed_train <- T
  }
  
  exp_config$train_sample_pctg <- experiments.parameters.cell_number[[train_dataset]]$train_sample_pctg
  exp_config$train_sample_num <- experiments.parameters.cell_number[[train_dataset]]$train_sample_num
  exp_config$test_sample_pctg <- experiments.parameters.cell_number[[test_dataset]]$test_sample_pctg
  exp_config$test_sample_num <- experiments.parameters.cell_number[[test_dataset]]$test_sample_num
  exp_config$train_sampling <- experiments.parameters.cell_number[[train_dataset]]$train_sampling
  exp_config$test_sampling <- experiments.parameters.cell_number[[test_dataset]]$test_sampling
  exp_config
}


experiments.config.update.sequencing_depth <-function(train_dataset, test_dataset=NULL,exp_config){
  exp_config$train_quantiles <- experiments.parameters.sequencing_depth[[train_dataset]]$train_quantiles
  exp_config$train_sample_pctg <- experiments.parameters.sequencing_depth[[train_dataset]]$train_sample_pctg
  exp_config$test_quantiles <- experiments.parameters.sequencing_depth[[test_dataset]]$test_quantiles
  exp_config$test_sample_pctg <- experiments.parameters.sequencing_depth[[test_dataset]]$test_sample_pctg
  exp_config
}


experiments.config.update.batch_effects <- function(train_dataset, test_dataset=NULL,exp_config){
  exp_config$train_sample_num <- experiments.parameters.batch_effects[[train_dataset]]$train_sample_num
  exp_config$train_sample_pctg <- experiments.parameters.batch_effects[[train_dataset]]$train_sample_pctg
  exp_config$test_sample_num <- experiments.parameters.batch_effects[[test_dataset]]$test_sample_num
  exp_config$test_sample_pctg <- experiments.parameters.batch_effects[[test_dataset]]$test_sample_pctg
  exp_config
}


experiments.config.update.inter_diseases <- function(train_dataset, test_dataset=NULL, exp_config){
  experiments.parameters$inter_diseases$train_sample_num <<- experiments.parameters.inter_diseases[[train_dataset]]$train_sample_num
  experiments.parameters$inter_diseases$train_sample_pctg <<- experiments.parameters.inter_diseases[[train_dataset]]$train_sample_pctg
  experiments.parameters$inter_diseases$test_sample_num <<- experiments.parameters.inter_diseases[[test_dataset]]$test_sample_num
  experiments.parameters$inter_diseases$test_sample_pctg <<- experiments.parameters.inter_diseases[[test_dataset]]$test_sample_pctg
}

experiments.config.update.train_test_datasets <- function(experiment, train_dataset, test_dataset){
  experiments.assign.data$train_dataset[[experiment]] <<- train_dataset
  print(str_glue("assign train data is {experiments.assign.data$train_dataset[[experiment]]}"))
  experiments.assign.data$test_dataset[[experiment]] <<- test_dataset
  print(str_glue("assign test data is {experiments.assign.data$test_dataset[[experiment]]}"))
}

experiments.config.check_config <- function(exp_config){
  stopifnot(!(purrr::is_null(exp_config$use_intra_dataset)&purrr::is_null(exp_config$use_inter_dataset)))
  stopifnot(xor(exp_config$use_intra_dataset,exp_config$use_inter_dataset))
  if(exp_config$use_intra_dataset){
    stopifnot(!purrr::is_null(exp_config$intra_dataset))
  }else{
    stopifnot(!purrr::is_null(exp_config$inter_dataset))
  }
  if(!purrr::is_null(exp_config$cv)){
    if(exp_config$cv)
    {stopifnot(!purrr::is_null(exp_config$intra_dataset))}
  }
}