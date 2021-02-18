
####Wrapper
constructor.data_constructor <- function(data,config,experiment,if_train=TRUE){
  require(SingleCellExperiment)
  #stopifnot(is(data,"SingleCellExperiment"))
  switch(experiment,
    simple_accuracy = constructor.simple_accuracy(data,config,if_train),
    cell_number = constructor.cell_number(data,config,if_train),
    sequencing_depth = constructor.sequencing_depth(data,config,if_train),
    celltype_structure = constructor.celltype_structure(data,config,if_train),
    batch_effects = constructor.batch_effects(data,config,if_train),
    inter_diseases = constructor.inter_diseases(data,config,if_train),
    celltype_complexity = constructor.celltype_complexity(data,config,if_train),
    inter_species = constructor.inter_species(data,config,if_train),
    random_noise = constructor.random_noise(data,config,if_train),
    inter_protocol = constructor.inter_protocol(data,config,if_train),
    stop("Can't constructor dataset for unkown experiments")
  )
}

constructor.base <- function(data,config,if_train){
  if(if_train){
    data <- utils.filter(data)
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg = NULL, types=unique(colData(data)$label))
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$train_sample_pctg, types=unique(colData(data)$label))
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE) 
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg = NULL, types=unique(colData(data)$label))
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$test_sample_pctg, types=unique(colData(data)$label))
    }
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}

constructor.simple_accuracy <- function(data,config, if_train){
  constructor.base(data,config,if_train)
}

constructor.cell_number <- function(data,config, if_train){
  constructor.base(data,config,if_train)
}

constructor.sequencing_depth <- function(data,config,if_train){
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

constructor.batch_effects <- function(data,config,if_train){
  if(purrr::is_list(data)&length(data)>1){
    data <- utils.combine_SCEdatasets(data,if_combined = TRUE)
  }else if(purrr::is_list(data)&length(data)==1){
    data <- data[[1]]
  }
  if(if_train){
    data <- utils.filter(data)
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg=NULL, types=unique(colData(data)$label))
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg=config$train_sample_pctg, types=unique(colData(data)$label))
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg=NULL, types=unique(colData(data)$label))
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg=config$test_sample_pctg, types=unique(colData(data)$label))
    }
  }
  data[!duplicated(rownames(data)),!duplicated(colnames(data))]
  
}

constructor.celltype_structure <- function(data,config,if_train){
  if(if_train){
    data <- utils.filter(data)
    if(!purrr::is_null(config$train_sample_num) && !is.na(config$train_sample_num)){
      print(str_glue("start sampling train data: {config$train_sample_num} per cell type"))
      data <- utils.sampler(data, sample_num=config$train_sample_num, sample_pctg = NULL, types=unique(colData(data)[[config$level]]),column = config$level)
    }else{
      print(str_glue("start sampling train data: {config$train_sample_pctg} percentage per cell type"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$train_sample_pctg, types=unique(colData(data)[[config$level]]),column = config$level)
    }
  }else{
    data <- utils.filter(data,filter_gene=FALSE) 
    if(!purrr::is_null(config$test_sample_num) && !is.na(config$test_sample_num)){
      print(str_glue("start sampling test data: {config$test_sample_num} per cell type"))
      data <- utils.sampler(data, sample_num=config$test_sample_num, sample_pctg = NULL, types=unique(colData(data)[[config$level]]),column = config$level)
    }else{
      print(str_glue("start sampling test data: {config$test_sample_pctg} percentage per cell type"))
      data <- utils.sampler(data, sample_num=NULL, sample_pctg = config$test_sample_pctg, types=unique(colData(data)[[config$level]]),column = config$level)
    }
  }
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

constructor.inter_diseases <- function(data,config,if_train){
  constructor.base(data,config,if_train)
}

constructor.celltype_complexity <- function(data,config,if_train){
  
}

constructor.inter_species <- function(data,config,if_train){
  
}

constructor.random_noise <- function(data,config,if_train){
  
}

constructor.inter_protocol <- function(data,config,if_train){
  
}
