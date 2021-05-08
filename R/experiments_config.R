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
  inter_protocol = experiments.methods.base_config,
  imbalance_impacts = experiments.methods.base_config
)



######simple accuracy config
experiments.parameters.simple_accuracy.base_config <- list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.PBMC <- list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.pancreas <- list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.ADASD <- list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL)

experiments.parameters.simple_accuracy.base_config.midbrain <- experiments.parameters.simple_accuracy.base_config

experiments.parameters.simple_accuracy.base_config.cellbench <- experiments.parameters.simple_accuracy.base_config


experiments.parameters.simple_accuracy <- list(
  PBMC_AllCells_withLabels.RDS = list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.simple_accuracy.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.simple_accuracy.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.simple_accuracy.base_config.ADASD,
  ADASD_autism.RDS = list(train_sample_pctg=0.01,test_sample_pctg=0.005,
                          cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI'), marker_gene_file=NULL),
  midbrain_human.RDS = experiments.parameters.simple_accuracy.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.simple_accuracy.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.simple_accuracy.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.simple_accuracy.base_config.cellbench,
  cellbench_Dropseq.RDS = experiments.parameters.simple_accuracy.base_config.cellbench
)


#####cell number config
experiments.parameters.cell_number.base_config <- list(train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                       marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.PBMC <- list(train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                            marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.pancreas <- list(train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                                marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.ADASD <- list(train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                                             marker_gene_file=NULL,trained=F)

experiments.parameters.cell_number.base_config.midbrain <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number.base_config.cellbench <- experiments.parameters.cell_number.base_config

experiments.parameters.cell_number <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                                      marker_gene_file=NULL,trained=F),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.cell_number.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.cell_number.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.cell_number.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.cell_number.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.cell_number.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.cell_number.base_config.pancreas,
  ADASD_AD.RDS = experiments.parameters.cell_number.base_config.ADASD,
  ADASD_autism.RDS = list(train_sampling=F, test_sampling=T, cv=F,cv_fold=NULL,metrics=c('ARI','AMI','FMI'), 
                          marker_gene_file=NULL,trained=F),
  midbrain_human.RDS = experiments.parameters.cell_number.base_config.midbrain,
  midbrain_mouse.RDS = experiments.parameters.cell_number.base_config.midbrain,
  cellbench_10x.RDS = experiments.parameters.cell_number.base_config.cellbench,
  cellbench_CELseq2.RDS = experiments.parameters.cell_number.base_config.cellbench,
  cellbench_Dropseq=experiments.parameters.cell_number.base_config.cellbench
)

####sequencing depth config
experiments.parameters.sequencing_depth.base_config <- list(cv=F, cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
                                                            marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.PBMC <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth.base_config.pancreas <- list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
                                                                     marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.ADASD <- list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
                                                                  marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.midbrain <- list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
                                                                     marker_gene_file=NULL)

experiments.parameters.sequencing_depth.base_config.cellbench <- experiments.parameters.sequencing_depth.base_config

experiments.parameters.sequencing_depth <- list(
  PBMC_AllCells_withLabels.RDS = experiments.parameters.sequencing_depth.base_config.ADASD,
  GSE96583_batch1_3_samples.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.sequencing_depth.base_config.PBMC,
  Muraro_pancreas_clean.RDS = list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
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
experiments.parameters.celltype_structure.base_config <- list(train_sample_pctg=0.8,test_sample_pctg=1,
                                                            cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                                            structure_file="PBMC_celltypes.csv",marker_gene_file=NULL)

experiments.parameters.celltype_structure.base_config.PBMC <- list(train_sample_pctg=0.2,test_sample_pctg=0.2,
                                                                   cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                                                   structure_file="PBMC_celltypes.csv",marker_gene_file=NULL)

experiments.parameters.celltype_structure.base_config.pancreas <- experiments.parameters.celltype_structure.base_config

experiments.parameters.celltype_structure.base_config.ADASD <- list(train_sample_pctg=0.05,test_sample_pctg=0.05,
                                                                    cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),
                                                                    structure_file="PBMC_celltypes.csv",marker_gene_file=NULL)

experiments.parameters.celltype_structure.base_config.midbrain <- experiments.parameters.celltype_structure.base_config

experiments.parameters.celltype_structure.base_config.cellbench <- experiments.parameters.celltype_structure.base_config

experiments.parameters.celltype_structure <- list(
  PBMC_AllCells_withLabels.RDS = list(train_sample_pctg=0.05,test_sample_pctg=0.05,
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

experiments.parameters.batch_effects.base_config  <- list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                          marker_gene_file=NULL,fixed_train=T,fixed_test=T)


experiments.parameters.batch_effects.base_config.PBMC <- list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                              marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.pancreas <- list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                                  marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.ADASD <- list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
                                                               marker_gene_file=NULL,fixed_train=T,fixed_test=T)

experiments.parameters.batch_effects.base_config.midbrain <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects.base_config.cellbench <- experiments.parameters.batch_effects.base_config

experiments.parameters.batch_effects <- list(
  PBMC_AllCells_withLabels.RDS = list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),
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
experiments.parameters.inter_diseases.base_config = list(cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)


experiments.parameters.inter_diseases.base_config.PBMC <- list(cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.pancreas <- list(cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.ADASD <- list(cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL)

experiments.parameters.inter_diseases.base_config.midbrain <- experiments.parameters.inter_diseases.base_config

experiments.parameters.inter_diseases.base_config.cellbench <- experiments.parameters.inter_diseases.base_config

experiments.parameters.inter_diseases <- list(
  PBMC_AllCells_withLabels.RDS = list(cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL),
  GSE96583_batch1_3_samples.RDS = experiments.parameters.inter_diseases.base_config.PBMC,
  GSE96583_8_Stim_Pats.RDS = experiments.parameters.inter_diseases.base_config.PBMC,
  GSE96583_8_Ctrl_Pats.RDS = experiments.parameters.inter_diseases.base_config.PBMC,
  Muraro_pancreas_clean.RDS = experiments.parameters.inter_diseases.base_config.pancreas,
  Segerstolpe_pancreas_clean.RDS = experiments.parameters.inter_diseases.base_config.pancreas,
  Xin_pancreas_clean.RDS = experiments.parameters.inter_diseases.base_config.pancreas,
  ADASD_AD.RDS =  experiments.parameters.inter_diseases.base_config.ADASD,
  ADASD_autism.RDS = list(cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL),
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
                       use_inter_dataset=F,inter_dataset=NULL,target_train_num=1200, target_test_num=NULL),
  
  cell_number=list( cv=F,cv_fold=NULL, metrics=c('ARI','AMI','FMI'), batch_free=F,fixed_train=T,fixed_test=F,
                   marker_gene_file=NULL,trained=F,target_train_num=1200, target_test_num=800,
                   train_sample_start=200, test_sample_start=100,train_sample_increment=400,test_sample_increment=250,
                   test_num=4, use_intra_dataset=F,intra_dataset=list(),
                   use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                          dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC5,
                                                          dataset.interdatasets$PBMC6,dataset.interdatasets$PBMC8,
                                                          dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                          dataset.interdatasets$pancreas5,dataset.interdatasets$ADASD2
                                                          )),

  
  sequencing_depth=list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI'),fixed_train=T,fixed_test=F,
                        marker_gene_file=NULL,batch_free=F,target_train_num=1200, target_test_num=800,test_num=5,
                        use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                               dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC5,
                                                               dataset.interdatasets$PBMC6,dataset.interdatasets$PBMC8,
                                                               dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                               dataset.interdatasets$pancreas5,dataset.interdatasets$ADASD2)),
  
  # celltype_structure=list(train_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_pctg,
  #                         train_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_num,
  #                         test_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_pctg,
  #                         test_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_num,
  #                         cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),batch_free=F,
  #                         structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  
  batch_effects=list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI'),target_train_num=1200, target_test_num=800,
                     marker_gene_file=NULL,fixed_train=T,fixed_test=T,batch_correct_method="combat_seq",use_intra_dataset=F,intra_dataset=list(),
                     use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC2,dataset.interdatasets$PBMC4,
                                                            dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                            dataset.interdatasets$pancreas3,dataset.interdatasets$pancreas4,
                                                            dataset.interdatasets$pancreas5,dataset.interdatasets$pancreas6)),
  
  imbalance_impacts = list(batch_free=F,target_train_num=1200, target_test_num=1000,fixed_train=T,fixed_test=F,
                           type_pctgs = list(list(0.2,0.2,0.2,0.2,0.2),
                                             list(0.4,0.3,0.15,0.1,0.05),
                                             list(0.6,0.2,0.1,0.05,0.05),
                                             list(0.8,0.1,0.05,0.03,0.02),
                                             list(0.9,0.04,0.03,0.02,0.01)),
                            cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL,use_intra_dataset=F,intra_dataset=list(),
                            use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                                   dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC5,
                                                                   dataset.interdatasets$PBMC6,dataset.interdatasets$PBMC8,
                                                                   dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                                   dataset.interdatasets$pancreas5,dataset.interdatasets$ADASD2)),
  
  inter_diseases = list(batch_free=F,target_train_num=1200, target_test_num=800,
                        cv=FALSE,metrics=c('ARI','AMI','FMI'),marker_gene_file=NULL,use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC5,dataset.interdatasets$PBMC6,
                                                               dataset.interdatasets$PBMC8,dataset.interdatasets$PBMC9,
                                                               dataset.interdatasets$PBMC11,dataset.interdatasets$PBMC12,
                                                               dataset.interdatasets$ADASD1,dataset.interdatasets$ADASD2)),
  
  celltype_number=list( cv=F,cv_fold=NULL, metrics=c('ARI','AMI','FMI'), batch_free=F,fixed_train=T,fixed_test=F,
                        marker_gene_file=NULL,trained=F,target_train_num=1200, target_test_num=800,
                        test_num=4, use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                               dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC4,
                                                               dataset.interdatasets$PBMC7,dataset.interdatasets$PBMC10,
                                                               dataset.interdatasets$pancreas3,dataset.interdatasets$pancreas4)
  
                        ),

  celltype_detection=list( cv=F,cv_fold=NULL, metrics=c('ARI','AMI','FMI'), batch_free=F,fixed_train=T,fixed_test=F,
                           marker_gene_file=NULL,trained=F,target_train_num=1200, target_test_num=800,
                           test_num=4, use_intra_dataset=F,intra_dataset=list(),
                           use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC1,dataset.interdatasets$PBMC2,
                                                                  dataset.interdatasets$PBMC3,dataset.interdatasets$PBMC4,
                                                                  dataset.interdatasets$PBMC7,dataset.interdatasets$PBMC10,
                                                                  dataset.interdatasets$pancreas1,dataset.interdatasets$pancreas2,
                                                                  dataset.interdatasets$pancreas5)
                          ),
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
         imbalance_impacts = experiments.config.update.imbalance_impacts(train_dataset, test_dataset, exp_config),
         stop("Unkown experiments")
         )
}


experiments.config.update.simple_accuracy <-function(train_dataset, test_dataset=NULL, exp_config){
  # train_test_sample_pctgs <- utils.calculate_sampling_pctg(train_dataset, test_dataset,exp_config)
  # train_sample_pctg <- train_test_pctgs[[1]]
  # test_sample_pctg <- train_test_pctgs[[2]]
  # exp_config$train_sample_pctg <- experiments.parameters.simple_accuracy[[train_dataset]]$train_sample_pctg
  # exp_config$test_sample_pctg <- experiments.parameters.simple_accuracy[[test_dataset]]$test_sample_pctg
  # exp_config
}

experiments.config.update.cell_number <-function(train_dataset, test_dataset=NULL, exp_config){
  calculate_increment <- function(dataset, increment, exp_config, if_train){
    total_dataset_num <- dataset.properties[[dataset]]$total_num
    if(exp_config$use_intra_dataset){
      print("calculate sampling increment for intra-dataset")
      if(if_train){
        increment <- min(increment, 0.3*0.7*total_dataset_num)
        print(str_glue("{dataset} train sample increment={increment}"))
      }else{
        increment <- min(increment, 0.3*0.3*total_dataset_num)
        print(str_glue("{dataset} test sample increment={increment}"))
      }
    
    }else{
      print("calculate sampling increment for inter-dataset")
      increment <- min(increment, 0.3*total_dataset_num)
      print(str_glue("{dataset} sample increment={increment}"))
    }
    increment
  }
  exp_config$trained <- F
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
    test_dataset_name <- str_split(test_dataset,"\\.")[[1]][[1]]
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- calculate_increment(test_dataset_name, exp_config$test_sample_increment, exp_config, if_train=F)
    }
    exp_config$target_test_num <- exp_config$test_sample_start + exp_config$current_increment_index*exp_config$increment
  }
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==0) F else T
    train_dataset_name <- str_split(train_dataset,"\\.")[[1]][[1]]
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- calculate_increment(train_dataset_name, exp_config$train_sample_increment, exp_config, if_train=T)
    }
    exp_config$target_train_num <- exp_config$train_sample_start + exp_config$current_increment_index*exp_config$increment
  }
  exp_config
}


experiments.config.update.sequencing_depth <-function(train_dataset, test_dataset=NULL,exp_config){
  
  calculate_quantile_increment <- function(dataset, exp_config, if_train){
    total_dataset_num <- dataset.properties[[dataset]]$total_num
    print("calculate sampling increment for intra-dataset")
    if(if_train){
      target_train_num <- exp_config$target_train_num
      increment <- min(1/(exp_config$test_num+1), target_train_num/total_dataset_num)
      print(str_glue("{dataset} train sample quantile increment={increment}"))
    }else{
      target_test_num <- exp_config$target_test_num
      increment <- min(1/(exp_config$test_num+1), target_test_num/total_dataset_num)
      print(str_glue("{dataset} test sample quantile increment={increment}"))
    }
    increment
  }
  exp_config$trained <- F
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
    test_dataset_name <- str_split(test_dataset,"\\.")[[1]][[1]]
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- calculate_quantile_increment(test_dataset_name, exp_config, if_train=F)
    }
    print("test dataset={test_dataset_name}")
  }
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==0) F else T
    train_dataset_name <- str_split(train_dataset,"\\.")[[1]][[1]]
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- calculate_increment(train_dataset_name, exp_config$train_sample_increment, exp_config, if_train=T)
    }
    print("train dataset={train_dataset_name}")
  }
  exp_config$high_quantile <- exp_config$current_increment_index*(1/(exp_config$test_num))+exp_config$increment
  exp_config$low_quantile <- exp_config$current_increment_index*(1/(exp_config$test_num))
  print(str_glue("high quantile is {exp_config$high_quantile}"))
  print(str_glue("low quantile is {exp_config$low_quantile}"))
  exp_config
}


experiments.config.update.batch_effects <- function(train_dataset, test_dataset=NULL,exp_config){
  # exp_config$train_sample_pctg <- experiments.parameters.batch_effects[[train_dataset]]$train_sample_pctg
  # exp_config$test_sample_pctg <- experiments.parameters.batch_effects[[test_dataset]]$test_sample_pctg
  # exp_config
}


experiments.config.update.inter_diseases <- function(train_dataset, test_dataset=NULL, exp_config){
  # exp_config$train_sample_pctg <- experiments.parameters.inter_diseases[[train_dataset]]$train_sample_pctg
  # exp_config$test_sample_pctg <- experiments.parameters.inter_diseases[[test_dataset]]$test_sample_pctg
  # exp_config
}

experiments.config.update.imbalance_impacts <- function(train_dataset, test_dataset=NULL, exp_config){
  exp_config$type_pctg <- exp_config$type_pctgs[[exp_config$current_increment_index+1]]
  exp_config$trained <- F
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
  }
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==0) F else T
  }

  exp_config
}


experiments.config.update.celltype_number <-function(train_dataset, test_dataset=NULL, exp_config){
  
  get_type_order <- function(data){
    type_pctg <- as.data.frame(table(colData(data)$label)/dim(data)[2])
    celltype_order <- as.vector(type_pctg[order(type_pctg[,2],decreasing = T),][,1])
    celltype_order
  }
  
  train_data <- utils.load_datasets(train_dataset) 
  test_data <- utils.load_datasets(test_dataset) 
  train_type <- unique(colData(train_data)$label)
  common_type <- unique(intersect(colData(train_data)$label,colData(test_data)$label))
  
  calculate_type_increment <- function(train_dataset,test_dataset, exp_config, if_train){
    if(exp_config$fixed_train) increment <- length(common_type)/5
    if(exp_config$fixed_test){
      increment <- ifelse(train_dataset %in% "Segerstolpe_pancreas_clean.RDS",round((length(train_type)-length(common_type))/4),length(common_type)/5)
    } 
    return(increment)
  }
  
  exp_config$trained <- F
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
    test_dataset_name <- str_split(test_dataset,"\\.")[[1]][[1]]
    if(purrr::is_null(exp_config$increment)){
      increment <- calculate_type_increment(train_dataset,test_dataset, exp_config, if_train=F)
    }
    exp_config$type_num <- floor(exp_config$current_increment_index*increment)+round(length(common_type)*0.2)
    celltype_order <- intersect(get_type_order(test_data),common_type)
    exp_config$sample_type <- celltype_order[1:(exp_config$type_num)]  
    print("test dataset={test_dataset_name}")
  }
  if(exp_config$fixed_test){
    train_dataset_name <- str_split(train_dataset,"\\.")[[1]][[1]]
    increment <- calculate_type_increment(train_dataset,test_dataset, exp_config, if_train=T)
    type_num <- (exp_config$current_increment_index-1)*increment
    exp_config$type_num <- type_num+length(common_type)
    exp_config$sample_type <- c(setdiff(train_type,common_type)[0:type_num],common_type) 
    print("train dataset={train_dataset_name}")
  }
  exp_config
}



experiments.config.update.celltype_detection <-function(train_dataset, test_dataset=NULL, exp_config){
  exp_config$trained <- F
  train_data <- utils.load_datasets(train_dataset) 
  train_type <- unique(colData(train_data)$label)
  test_data <- utils.load_datasets(test_dataset)
  test_type <- unique(colData(test_data)$label)
  
  if(exp_config$fixed_train){
    exp_config$trained <- T
    test_dataset_name <- str_split(test_dataset,"\\.")[[1]][[1]]
    print("test dataset={test_dataset_name}")
  }
  if(exp_config$fixed_test){
    train_dataset_name <- str_split(train_dataset,"\\.")[[1]][[1]]
    increment <- calculate_type_increment(train_dataset, exp_config, if_train=T)
    print("train dataset={train_dataset_name}")
  }
  unknown_type <- setdiff(test_type,train_type)
  exp_config$unknown_type <- unknown_type[!is.na(unknown_type)]
  exp_config
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