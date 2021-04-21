
####Wrapper
constructor.data_constructor <- function(data,config,experiment,if_train=TRUE,sample_seed=NULL){
  require(SingleCellExperiment)
  if(purrr::is_list(data)&length(data)>1)
    stopifnot(all(map_lgl(data,~is(.,class2="SingleCellExperiment"))))
  else
    stopifnot(is(data,"SingleCellExperiment"))
  switch(experiment,
         simple_accuracy = constructor.simple_accuracy(data,config,if_train,sample_seed),
         cell_number = constructor.cell_number(data,config,if_train,sample_seed),
         sequencing_depth = constructor.sequencing_depth(data,config,if_train,sample_seed),
         celltype_structure = constructor.celltype_structure(data,config,if_train,sample_seed),
         batch_effects = constructor.batch_effects(data,config,if_train,sample_seed),
         inter_diseases = constructor.inter_diseases(data,config,if_train,sample_seed),
         celltype_complexity = constructor.celltype_complexity(data,config,if_train,sample_seed),
         inter_species = constructor.inter_species(data,config,if_train,sample_seed),
         random_noise = constructor.random_noise(data,config,if_train,sample_seed),
         inter_protocol = constructor.inter_protocol(data,config,if_train,sample_seed),
         celltype_number = constructor.celltype_number(data,config,if_train,sample_seed),
         celltype_detection = constructor.celltype_detection(data,config,if_train,sample_seed),
         rare_celltype = constructor.rare_celltype(data,config,if_train,sample_seed),
         stop("Can't constructor dataset for unkown experiments")
  )
}


constructor.base <- function(data,config,if_train,sample_seed=NULL){
  if(if_train){
    data <- utils.filter(data)
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg = NULL, types=unique(colData(data)$label), sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$train_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg = NULL, types=unique(colData(data)$label), sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$test_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
    }
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}


constructor.simple_accuracy <- function(data,config, if_train,sample_seed=NULL){
  constructor.base(data,config,if_train,sample_seed)
}

constructor.cell_number <- function(data,config, if_train,sample_seed=NULL){
  constructor.base(data,config,if_train,sample_seed)
}

constructor.sequencing_depth <- function(data,config,if_train,sample_seed=NULL){
  quantile <- config$quantile
  if_right <- config$right ###if take the deeper sequencing part
  if(if_train){
    data <- utils.filter(data) %>%
      utils.seqDepthSelector(quantile=quantile,right=if_right)
  }else{
    data <- utils.filter(data,filter_gene=FALSE) %>%
      utils.seqDepthSelector(quantile=quantile,right=if_right)
  }
  data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}

constructor.batch_effects <- function(data,config,if_train,sample_seed=NULL){
  if(purrr::is_list(data)&length(data)>1){
    data <- utils.combine_SCEdatasets(data,if_combined = TRUE)
  }else if(purrr::is_list(data)&length(data)==1){
    data <- data[[1]]
  }
  data <- constructor.base(data,config,if_train,sample_seed)
  data
}

constructor.rare_celltype <- function(data,config,if_train,sample_seed=NULL){
  if(if_train){
    data <- utils.filter(data)
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg = NULL, types=unique(colData(data)$label), sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$train_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg = NULL, types=unique(colData(data)$label), sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$test_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
    }
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}



constructor.celltype_detection <- function(data,config,if_train,sample_seed=NULL){
  if(if_train){
    data <- utils.filter(data)
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg = NULL, types=config$train_type, sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$train_sample_pctg, types=config$train_type, sample_seed=sample_seed)
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg = NULL, types=unique(colData(data)$label), sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$test_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
    }
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}


constructor.celltype_number <- function(data,config,if_train,sample_seed=NULL){
  if(if_train){
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg = NULL, types=config$sample_celltype, sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$train_sample_pctg, types=config$sample_celltype, sample_seed=sample_seed)
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg = NULL, types=config$sample_celltype, sample_seed=sample_seed)
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type, sample_seed is {sample_seed}"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$test_sample_pctg, types=config$sample_celltype, sample_seed=sample_seed)
    }
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}



constructor.celltype_structure <- function(data,config,if_train,sample_seed=NULL){
  data <- constructor.base(data,config,if_train,sample_seed)
  colData(data)$label <- colData(data)[[config$level]]
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}


constructor.type_architecturer <- function(config,dataset){
  data <- utils.load_datasets(dataset) 
  structure_file <- str_glue("{type_home}/{config$structure_file}")
  rule <- read.csv(structure_file,stringsAsFactors = FALSE)
  for(i in 2:ncol(rule)){
    missing_idx <- which(unlist(purrr::map(rule[i],str_length))==0)
    if(length(missing_idx)>0){
      rule[missing_idx,i] <- rule[missing_idx,i-1] ####using previous level label to fill current level
    }
    if(length(unique(rule[,i]))<=2){
      print(str_glue("{colnames(rule)[[i]]} has cell types fewer than 2, ignoring this level"))
    }
  }
  
  rule <- dplyr::select(rule,-colnames(rule)[which(purrr::map_lgl(colnames(rule),~{return(length(unique(rule[,.]))<=2)}))])

  newColData <- merge(colData(data),rule,by.x='label',by.y='Cell') 
  rownames(newColData) <- rownames(colData(data))
  colData(data) <- newColData
  data
}

constructor.inter_diseases <- function(data,config,if_train,sample_seed=NULL){
  constructor.base(data,config,if_train,sample_seed)
}

constructor.celltype_complexity <- function(data,config,if_train,sample_seed=NULL){
  
}

constructor.inter_species <- function(data,config,if_train,sample_seed=NULL){
  
}

constructor.random_noise <- function(data,config,if_train,sample_seed=NULL){
  
}

constructor.inter_protocol <- function(data,config,if_train,sample_seed=NULL){
  
}
