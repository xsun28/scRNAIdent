experiment <- "imbalance_impacts"


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
experiments.methods.base_config <- list(cluster=c('raceID3','seurat_clustering','tscan','sc3','liger','cidr','monocle3','pcaReduce'),assign=c("seurat_mapping","cellassign",'scmap_cluster','scmap_cell','chetah','singlecellnet','garnett','singleR'))
experiments.methods <- list(
  simple_accuracy=experiments.methods.base_config, 
  cell_number=experiments.methods.base_config,
  celltype_number = experiments.methods.base_config,
  sequencing_depth=experiments.methods.base_config,
  celltype_structure=experiments.methods.base_config,
  batch_effects=experiments.methods.base_config, 
                       # list(cluster_batch_free=c('seurat',"sc3",'tscan','liger'), assign_batch_free=c('chetah','garnett'), marker_gene_assign_batch_free=c("cellassign"))),
  sample_bias = experiments.methods.base_config,
  celltype_complexity = experiments.methods.base_config,
  inter_species = experiments.methods.base_config,
  random_noise = experiments.methods.base_config, 
  inter_protocol = experiments.methods.base_config,
  imbalance_impacts = experiments.methods.base_config,
  unknown_types = experiments.methods.base_config
)

##########
experiments.parameters <- list(
  simple_accuracy=list(cv=TRUE, cv_fold=5,metrics=c('ARI','AMI','FMI','v_measure'),batch_free=F,
                       marker_gene_file=NULL,use_intra_dataset=T,intra_dataset=dataset.datasets,
                       use_inter_dataset=F,inter_dataset=NULL,target_train_num=1200, target_test_num=NULL),
  
  cell_number=list( cv=F,cv_fold=NULL, metrics=c('ARI','AMI','FMI','v_measure'), batch_free=F,fixed_train=T,fixed_test=F,
                   marker_gene_file=NULL,trained=F,target_train_num=1200, target_test_num=800,
                   train_sample_start=200, test_sample_start=100,train_sample_increment=400,test_sample_increment=250,
                   test_num=4, use_intra_dataset=F,intra_dataset=list(),
                   use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC9,dataset.interdatasets$PBMC24,
                                                          dataset.interdatasets$PBMC17,dataset.interdatasets$PBMC21,
                                                          dataset.interdatasets$PBMC13,dataset.interdatasets$PBMC29,
                                                          dataset.interdatasets$PBMC5,dataset.interdatasets$PBMC1,
                                                          dataset.interdatasets$pancreas2,dataset.interdatasets$pancreas3,
                                                          dataset.interdatasets$pancreas6,dataset.interdatasets$ADASD2
                                                          )),

  
  sequencing_depth=list(cv=F,cv_fold=5,metrics=c('ARI','AMI','FMI','v_measure'),fixed_train=T,fixed_test=F,
                        marker_gene_file=NULL,batch_free=F,target_train_num=1200, target_test_num=800,test_num=5,
                        use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC9,dataset.interdatasets$PBMC24,
                                                               dataset.interdatasets$PBMC17,dataset.interdatasets$PBMC21,
                                                               dataset.interdatasets$PBMC13,dataset.interdatasets$PBMC29,
                                                               dataset.interdatasets$PBMC5,dataset.interdatasets$PBMC1,
                                                               dataset.interdatasets$pancreas2,dataset.interdatasets$pancreas3,
                                                               dataset.interdatasets$pancreas6,dataset.interdatasets$ADASD2
                                                               )),
  
  # celltype_structure=list(train_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_pctg,
  #                         train_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$train_dataset$celltype_structure]]$train_sample_num,
  #                         test_sample_pctg=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_pctg,
  #                         test_sample_num=experiments.parameters.celltype_structure[[experiments.assign.data$test_dataset$celltype_structure]]$test_sample_num,
  #                         cv=TRUE,cv_fold=5,metrics=c('wNMI','wRI'),batch_free=F,
  #                         structure_file="PBMC_celltypes.csv",marker_gene_file=NULL),
  
  batch_effects=list(cv=FALSE,remove_batch=TRUE,metrics=c('ARI','AMI','FMI','v_measure'),target_train_num=1200, target_test_num=800,
                     marker_gene_file=NULL,fixed_train=T,fixed_test=T,batch_correct_method="MNN",use_intra_dataset=F,intra_dataset=list(),
                     use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC9,dataset.interdatasets$PBMC24,
                                                            dataset.interdatasets$PBMC17, dataset.interdatasets$PBMC1,
                                                            dataset.interdatasets$pancreas2,dataset.interdatasets$pancreas3,
                                                            dataset.interdatasets$pancreas6,dataset.interdatasets$ADASD2)),
  
  imbalance_impacts = list(batch_free=F,target_train_num=1200, target_test_num=1000,fixed_train=T,fixed_test=F,
                           all_type_pctgs = list(   "3"=list( list(0.33,0.33,0.34),
                                                            list(0.6,0.35,0.05),
                                                            list(0.9,0.095,0.005)
                                                            ),
                                                    "4"=list( list(0.25,0.25,0.25,0.25),
                                                            list(0.4,0.3,0.2,0.1),
                                                            list(0.6,0.2,0.15,0.05),
                                                            list(0.9,0.08,0.015,0.005)
                                                            ),
                                                    "5"=list( list(0.2,0.2,0.2,0.2,0.2),
                                                            list(0.4,0.3,0.15,0.1,0.05),
                                                            list(0.65,0.15,0.1,0.08,0.02),
                                                            list(0.8,0.1,0.05,0.04,0.01),
                                                            list(0.9,0.04,0.03,0.025,0.005)
                                                            )
                                             ),
                            cv=FALSE,metrics=c('ARI','AMI','FMI','v_measure',"BCubed"),marker_gene_file=NULL,use_intra_dataset=F,intra_dataset=list(),
                            use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC9,dataset.interdatasets$PBMC24,
                                                                   dataset.interdatasets$PBMC17,dataset.interdatasets$PBMC21,
                                                                   dataset.interdatasets$PBMC13,dataset.interdatasets$PBMC29,
                                                                   dataset.interdatasets$PBMC5,dataset.interdatasets$PBMC1,
                                                                   dataset.interdatasets$pancreas2,dataset.interdatasets$pancreas3,
                                                                   dataset.interdatasets$pancreas6,dataset.interdatasets$ADASD2
                                                                   )
                           ),
  
  sample_bias = list(batch_free=F,target_train_num=1200, target_test_num=1000, fixed_train=F,fixed_test=T,sample_seed = 10,
                        cv=FALSE,metrics=c('ARI','AMI','FMI','v_measure'),marker_gene_file=NULL,use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(list(datasets=dataset.interdatasets$PBMC34,inds=list(train_ind=1015,test_ind=1488)),
                                                               list(datasets=dataset.interdatasets$PBMC35,inds=list(train_ind=1015,test_ind=1488)),
                                                               list(datasets=dataset.interdatasets$PBMC5,inds=list(train_ind=1488,test_ind=1488)),
                                                               list(datasets=dataset.interdatasets$PBMC11,inds=list(train_ind=1488,test_ind=1488)),
                                                               list(datasets=dataset.interdatasets$PBMC23,inds=list(train_ind=1511,test_ind=1488)),
                                                               list(datasets=dataset.interdatasets$PBMC24,inds=list(train_ind=1511,test_ind=1488))
                                                               )
                     ),
  
  celltype_number=list( cv=F,cv_fold=NULL, metrics=c('ARI','AMI','FMI',"BCubed"), batch_free=F,fixed_train=T,fixed_test=F,fixed_total_cell_num=F,max_total_cell_num=1500,
                        marker_gene_file=NULL,trained=F,target_train_num=1200,target_test_num=1000,test_num=3, use_intra_dataset=F,intra_dataset=list(),
                        use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC9,dataset.interdatasets$PBMC24,
                                                               dataset.interdatasets$PBMC17,dataset.interdatasets$PBMC21,
                                                               dataset.interdatasets$PBMC13,dataset.interdatasets$PBMC29,
                                                               dataset.interdatasets$pancreas2,dataset.interdatasets$pancreas3,
                                                               dataset.interdatasets$pancreas6)
                        
                        ),
  

  unknown_types=list(cv=F,cv_fold=NULL, metrics=c('ARI','AMI','FMI',"BCubed"), batch_free=F,fixed_train=T,fixed_test=T,
                           marker_gene_file=NULL,trained=F,target_train_num=1200, target_test_num=800,
                           test_num=4, use_intra_dataset=F,intra_dataset=list(),
                           use_inter_dataset=T,inter_dataset=list(dataset.interdatasets$PBMC9,dataset.interdatasets$PBMC24,
                                                                  dataset.interdatasets$PBMC17,dataset.interdatasets$PBMC21,
                                                                  dataset.interdatasets$PBMC13,dataset.interdatasets$PBMC29,
                                                                  dataset.interdatasets$PBMC5,dataset.interdatasets$PBMC1,
                                                                  dataset.interdatasets$pancreas2,dataset.interdatasets$pancreas3,
                                                                  dataset.interdatasets$pancreas6,dataset.interdatasets$ADASD2
                                                                  )
                           ),
  celltype_complexity = list(),
  inter_species = list(),
  random_noise = list(),
  inter_protocol = list()
)
##########################
experiments.config.init <- function(experiment, train_dataset, test_dataset=NULL, exp_config){
  switch(experiment,
         simple_accuracy=experiments.config.init.simple_accuracy(train_dataset, test_dataset, exp_config),
         cell_number = experiments.config.init.cell_number(train_dataset, test_dataset, exp_config),
         celltype_number = experiments.config.init.celltype_number(train_dataset, test_dataset, exp_config),
         sequencing_depth = experiments.config.init.sequencing_depth(train_dataset, test_dataset, exp_config),
         celltype_structure = experiments.config.init.celltype_structure(train_dataset, test_dataset, exp_config),
         batch_effects = experiments.config.init.batch_effects(train_dataset, test_dataset, exp_config),
         sample_bias = experiments.config.init.sample_bias(train_dataset, test_dataset, exp_config),
         celltype_complexity = experiments.config.init.celltype_complexity(train_dataset, test_dataset, exp_config),
         inter_species = experiments.config.init.inter_species(train_dataset, test_dataset, exp_config),
         random_noise = experiments.config.init.random_noise(train_dataset, test_dataset, exp_config),
         inter_protocol = experiments.config.init.inter_protocol(train_dataset, test_dataset, exp_config),
         imbalance_impacts = experiments.config.init.imbalance_impacts(train_dataset, test_dataset, exp_config),
         unknown_types = experiments.config.init.unknown_types(train_dataset, test_dataset, exp_config),
         stop("Unkown experiments")
  )
}

experiments.config.init.base <- function(exp_config){
  exp_config$have_test <- if(exp_config$cv) F else T
  exp_config$trained <- F
  exp_config$clustered <- F
  if(purrr::is_null(exp_config$fixed_train)) exp_config$fixed_train <- F
  if(purrr::is_null(exp_config$fixed_test)) exp_config$fixed_test <- F
  exp_config
}


experiments.config.init.simple_accuracy <-function(train_dataset, test_dataset=NULL, exp_config){
  exp_config <- experiments.config.init.base(exp_config)
  exp_config
}

experiments.config.init.cell_number <-function(train_dataset, test_dataset=NULL, exp_config){
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
  exp_config <- experiments.config.init.base(exp_config)
  if(exp_config$fixed_train){
    if(purrr::is_null(exp_config$increment)){
      test_dataset_name <- dataset.name.map1[[test_dataset]]
      exp_config$increment <- calculate_increment(test_dataset_name, exp_config$test_sample_increment, exp_config, if_train=F)
    }
  }
  
  if(exp_config$fixed_test){
    if(purrr::is_null(exp_config$increment)){
      train_dataset_name <- dataset.name.map1[[train_dataset]]
      exp_config$increment <- calculate_increment(train_dataset_name, exp_config$train_sample_increment, exp_config, if_train=T)
    }
  }
  exp_config
}

experiments.config.init.sequencing_depth <-function(train_dataset, test_dataset=NULL,exp_config){
  
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
  exp_config <- experiments.config.init.base(exp_config)
  if(exp_config$fixed_train){
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- calculate_quantile_increment(test_dataset, exp_config, if_train=F)
    }
  }
  if(exp_config$fixed_test){
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- calculate_quantile_increment(train_dataset,exp_config, if_train=T)
    }
  }
  exp_config
}

experiments.config.init.imbalance_impacts <- function(train_dataset, test_dataset=NULL, exp_config){
  exp_config <- experiments.config.init.base(exp_config)
  train_test_types <- utils.get_train_test_types(train_dataset,test_dataset)
  train_type <- train_test_types$train_type
  test_type <- train_test_types$test_type
  train_test_common_type <- intersect(train_type,test_type)
  common_type_num <- length(train_test_common_type)
  if(common_type_num > 5){
    common_type_num <- 5
    train_data <- utils.load_datasets(train_dataset)
    train_test_common_type <- (dplyr::group_by(as.data.frame(colData(train_data)),label) %>% dplyr::summarize(type_number=n()) %>% dplyr::arrange(desc(type_number)))[['label']][1:common_type_num]
  }
  exp_config$train_kept_types <- train_test_common_type
  exp_config$test_kept_types <- train_test_common_type
  exp_config$type_pctgs <- exp_config$all_type_pctgs[[as.character(length(train_test_common_type))]]
  exp_config$test_num <- length(exp_config$type_pctgs)
  exp_config
}

experiments.config.init.batch_effects <- function(train_dataset, test_dataset=NULL,exp_config){
  exp_config <- experiments.config.init.base(exp_config)
  exp_config
}

experiments.config.init.sample_bias <- function(train_dataset, test_dataset=NULL, exp_config){
  exp_config <- experiments.config.init.base(exp_config)
  exp_config$cluster_results <- list()
  exp_config
}

experiments.config.init.unknown_types <-function(train_dataset, test_dataset=NULL, exp_config){
  
  dataset.properties$Muraro_pancreas$cell_types <<- c('beta','alpha','delta','epsilon',
                                                      'acinar','ductal','gamma','mesenchymal',
                                                      'endothelial')
  dataset.properties$Segerstolpe_pancreas$cell_types <<- c('beta','alpha','delta','epsilon','mast','MHC class II',
                                                           'acinar','ductal','gamma','mesenchymal','PSC',
                                                           'endothelial')
  exp_config <- experiments.config.init.base(exp_config)
  train_test_types <- utils.get_train_test_types(train_dataset,test_dataset)
  train_type <- train_test_types$train_type
  test_type <- train_test_types$test_type
  exp_config$common_type <- intersect(train_type,test_type)
  exp_config$train_type <- exp_config$common_type
  if(exp_config$fixed_test){
    exp_config$test_num <- min(length(exp_config$common_type),exp_config$test_num)
    exp_config$test_type <- exp_config$common_type
  }
  exp_config
}


experiments.config.init.celltype_number <- function(train_dataset, test_dataset=NULL, exp_config){
  get_type_order <- function(dataset){
    data <- utils.load_datasets(dataset) %>% utils.filter(filter_gene=F, filter_cells=TRUE, filter_cell_type=TRUE)
    type_pctg <- as.data.frame(table(colData(data)$label)/dim(data)[2])
    celltype_order <- as.vector(type_pctg[order(type_pctg[,2],decreasing = T),][,1])
    celltype_order
  }
  exp_config <- experiments.config.init.base(exp_config)
  exp_config$train_type <- get_type_order(train_dataset)
  exp_config$test_type <- get_type_order(test_dataset)
  exp_config$common_type <- intersect(exp_config$train_type,exp_config$test_type)
  exp_config$train_uniq_type <- setdiff(exp_config$train_type,exp_config$common_type)
  
  if(exp_config$fixed_train){
    test_dataset_name <- dataset.name.map1[[test_dataset]]
    if(!exp_config$fixed_total_cell_num){
      exp_config$test_sample_pctg <- exp_config$max_total_cell_num/dataset.properties[[test_dataset_name]]$total_num
    }
    
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- floor(length(exp_config$common_type)/exp_config$test_num)
      while(exp_config$increment==0){  
        exp_config$test_num <- exp_config$test_num-1
        stopifnot(exp_config$test_num>=2)
        exp_config$increment <- floor(length(exp_config$common_type)/exp_config$test_num)
      }
    }
    exp_config$train_number <- length(exp_config$train_type)
    exp_config$test_number <- 0
    print(str_glue("test dataset={test_dataset_name}"))
  }
  if(exp_config$fixed_test){
    exp_config$test_number <- ceiling(min(max(length(exp_config$common_type)*0.4,3),length(exp_config$common_type)))
    exp_config$test_type <- exp_config$common_type[1:exp_config$test_number]
    train_dataset_name <- dataset.name.map1[[train_dataset]]
    if(!exp_config$fixed_total_cell_num){
      exp_config$train_sample_pctg <- exp_config$max_total_cell_num/dataset.properties[[train_dataset_name]]$total_num
    }
    if(purrr::is_null(exp_config$increment)){
      exp_config$increment <- floor((length(exp_config$train_type)-exp_config$test_number)/(exp_config$test_num-1))
      while(exp_config$increment==0){  
        exp_config$test_num <- exp_config$test_num-1
        stopifnot(exp_config$test_num>=2)
        exp_config$increment <- floor((length(exp_config$train_type)-length(exp_config$test_number))/exp_config$test_num)
      }
    }

    print(str_glue("test dataset={train_dataset_name}"))
  }
  exp_config
}


########################
experiments.config.update <- function(experiment, train_dataset, test_dataset=NULL, exp_config,...){
  switch(experiment,
         simple_accuracy=experiments.config.update.simple_accuracy(train_dataset, test_dataset, exp_config),
         cell_number = experiments.config.update.cell_number(train_dataset, test_dataset, exp_config),
         celltype_number = experiments.config.update.celltype_number(train_dataset, test_dataset, exp_config),
         sequencing_depth = experiments.config.update.sequencing_depth(train_dataset, test_dataset, exp_config),
         celltype_structure = experiments.config.update.celltype_structure(train_dataset, test_dataset, exp_config),
         batch_effects = experiments.config.update.batch_effects(train_dataset, test_dataset, exp_config),
         sample_bias = experiments.config.update.sample_bias(train_dataset, test_dataset, exp_config,...),
         celltype_complexity = experiments.config.update.celltype_complexity(train_dataset, test_dataset, exp_config),
         inter_species = experiments.config.update.inter_species(train_dataset, test_dataset, exp_config),
         random_noise = experiments.config.update.random_noise(train_dataset, test_dataset, exp_config),
         inter_protocol = experiments.config.update.inter_protocol(train_dataset, test_dataset, exp_config),
         imbalance_impacts = experiments.config.update.imbalance_impacts(train_dataset, test_dataset, exp_config),
         unknown_types = experiments.config.update.unknown_types(train_dataset, test_dataset, exp_config),
         stop("Unkown experiments")
         )
}


experiments.config.update.simple_accuracy <-function(train_dataset, test_dataset=NULL, exp_config){
  exp_config
}

experiments.config.update.cell_number <-function(train_dataset, test_dataset=NULL, exp_config){
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
    test_dataset_name <- dataset.name.map1[[test_dataset]]
    exp_config$target_test_num <- exp_config$test_sample_start + exp_config$current_increment_index*exp_config$increment
  }
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==0) F else T
    train_dataset_name <- dataset.name.map1[[train_dataset]]
    exp_config$target_train_num <- exp_config$train_sample_start + exp_config$current_increment_index*exp_config$increment
  }
  exp_config
}


experiments.config.update.sequencing_depth <-function(train_dataset, test_dataset=NULL,exp_config){
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
    test_dataset_name <- dataset.name.map1[[test_dataset]]
    print("test dataset={test_dataset_name}")
  }
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==0) F else T
    train_dataset_name <- dataset.name.map1[[train_dataset]]
    print("train dataset={train_dataset_name}")
  }
  exp_config$high_quantile <- exp_config$current_increment_index*(1/(exp_config$test_num))+exp_config$increment
  exp_config$low_quantile <- exp_config$current_increment_index*(1/(exp_config$test_num))
  print(str_glue("high quantile is {exp_config$high_quantile}"))
  print(str_glue("low quantile is {exp_config$low_quantile}"))
  exp_config
}


experiments.config.update.batch_effects <- function(train_dataset, test_dataset=NULL,exp_config){
  exp_config
}


experiments.config.update.sample_bias <- function(train_dataset, test_dataset=NULL, exp_config,...){
  extra_args <- list(...)
  exp_config$train_ind <- extra_args$train_ind
  exp_config$test_ind <- extra_args$test_ind
  if(!purrr::is_null(exp_config$cluster_results[[test_dataset]])){
    exp_config$clustered <- T
    print(str_glue("clustered results for {test_dataset} already generated"))
  }else{
    print(str_glue("clustered results for {test_dataset} not generated"))
    exp_config$clustered <- F
  }
  exp_config
}

experiments.config.update.imbalance_impacts <- function(train_dataset, test_dataset=NULL, exp_config){
  exp_config$type_pctg <- exp_config$type_pctgs[[exp_config$current_increment_index+1]]
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==0) F else T
  }
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==0) F else T
  }

  exp_config
}


experiments.config.update.celltype_number <- function(train_dataset, test_dataset=NULL, exp_config){
  if(exp_config$fixed_train){
    exp_config$trained <- if(exp_config$current_increment_index==1) F else T
    test_dataset_name <- dataset.name.map1[[test_dataset]]
    exp_config$test_number <- exp_config$increment + exp_config$test_number
    exp_config$test_type <- (exp_config$common_type)[1:(exp_config$test_number)]
    print(str_glue("fixed train with test types={exp_config$test_type}"))
  }
  if(exp_config$fixed_test){
    train_dataset_name <- dataset.name.map1[[train_dataset]]
    exp_config$clustered <- if(exp_config$current_increment_index==1) F else T
    exp_config$train_number <- ifelse(exp_config$current_increment_index==1,exp_config$test_number,exp_config$train_number+exp_config$increment)
    exp_config$train_type <- c(exp_config$common_type,exp_config$train_uniq_type)[1:exp_config$train_number]

    print(str_glue("fixed test with train types={exp_config$train_type}"))
  }
  exp_config
}




experiments.config.update.unknown_types <-function(train_dataset, test_dataset=NULL, exp_config){
  
  if(exp_config$fixed_test){
    exp_config$clustered <- if(exp_config$current_increment_index==1) F else T
    test_dataset_name <- dataset.name.map1[[test_dataset]]
    exp_config$unknown_type <- exp_config$common_type[[exp_config$current_increment_index]]
    # exp_config$train_type <- exp_config$common_type[-exp_config$current_increment_index]
    print("test dataset={test_dataset_name}")
  }
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