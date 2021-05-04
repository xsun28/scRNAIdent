
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
    celltype_number = constructor.celltype_number(data,config,if_train,sample_seed),
    sequencing_depth = constructor.sequencing_depth(data,config,if_train,sample_seed),
    celltype_structure = constructor.celltype_structure(data,config,if_train,sample_seed),
    batch_effects = constructor.batch_effects(data,config,if_train,sample_seed),
    inter_diseases = constructor.inter_diseases(data,config,if_train,sample_seed),
    celltype_complexity = constructor.celltype_complexity(data,config,if_train,sample_seed),
    inter_species = constructor.inter_species(data,config,if_train,sample_seed),
    random_noise = constructor.random_noise(data,config,if_train,sample_seed),
    inter_protocol = constructor.inter_protocol(data,config,if_train,sample_seed),
    imbalance_impacts = constructor.imbalance_impacts(data,config,if_train,sample_seed),
    stop("Can't constructor dataset for unkown experiments")
  )
}

constructor.base <- function(data,config,if_train,sample_seed=NULL){
  if(if_train){
    data <- utils.filter(data)
    train_sample_pctg <- if(purrr::is_null(config$train_sample_pctg)) utils.calculate_sampling_pctg(data, config$target_train_num,config, if_train) else config$train_sample_pctg
    print(str_glue("start sampling train data: {train_sample_pctg} percentage per cell type"))
    if(!purrr::is_null(sample_seed)) print(str_glue("sample_seed is {sample_seed}"))
    data <- utils.sampler(data, sample_pctg = train_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
  
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    test_sample_pctg <- if(purrr::is_null(config$test_sample_pctg)) utils.calculate_sampling_pctg(data, config$target_test_num,config, if_train) else config$test_sample_pctg
    print(str_glue("start sampling test data: {test_sample_pctg} percentage per cell type"))
    if(!purrr::is_null(sample_seed)) print(str_glue("sample_seed is {sample_seed}"))
    data <- utils.sampler(data,  sample_pctg = test_sample_pctg, types=unique(colData(data)$label), sample_seed=sample_seed)
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}

constructor.simple_accuracy <- function(data,config, if_train,sample_seed=NULL){
  constructor.base(data,config,if_train,sample_seed)
}

constructor.cell_number <- function(data,config, if_train,sample_seed=NULL){
  constructor.base(data,config,if_train,sample_seed)
}


constructor.imbalance_impacts <- function(data,config,if_train,sample_seed=NULL){
  require(tidyverse)
  type_pctg_selector <- function(data, types, type_pctg, target_sample_num){
    data <- data[,colData(data)$label %in%max_types]
    total_sample <- length(colData(data)$label)
    if(target_sample_num > total_sample) target_sample_num <- total_sample
    params <- tibble(type=types,pctg=type_pctg)
    col_data <- colData(data)
    sampled <- purrr::pmap(params,function(type,pctg){
                                                      type_data <- rownames(col_data[(col_data$label==type),])
                                                      sample_num <- min(as.integer(pctg*target_sample_num),length(type_data))
                                                      sample(type_data,sample_num)
                                                     })
    data[,unlist(sampled)]
  }
  type_pctg <- config$type_pctg
  type_num <- length(type_pctg)
  if(if_train){
    if(config$fixed_train){
      print("in imbalance impacts, sample fixed train dataset")
      data <- constructor.base(data,config,if_train,sample_seed)
    }else{
      print("in imbalance impacts, sample unfixed train dataset")
      data <- utils.filter(data)
      total_types <- length(unique(colData(data)$label))
      if(type_num > total_types) {
        type_num <- total_types
        type_pctg <- type_pctg[1:type_num]
      }
      max_types <- (dplyr::group_by(as.data.frame(colData(data)),label) %>% dplyr::summarize(type_number=n()) %>% dplyr::arrange(desc(type_number)))[['label']][1:type_num]
      
      data <- type_pctg_selector(data, max_types, type_pctg, config$target_train_num)
    }
  }else{
    if(config$fixed_test){
      print("in imbalance impacts, sample fixed test dataset")
      data <- constructor.base(data,config,if_train,sample_seed)
    }else{
      print("in imbalance impacts, sample unfixed test dataset")
      data <- utils.filter(data,filter_gene=FALSE)
      total_types <- length(unique(colData(data)$label))
      if(type_num > total_types) {
        type_num <- total_types
        type_pctg <- type_pctg[1:type_num]
      }
      max_types <- (dplyr::group_by(as.data.frame(colData(data)),label) %>% dplyr::summarize(type_number=n()) %>% dplyr::arrange(desc(type_number)))[['label']][1:type_num]
      data <- type_pctg_selector(data, max_types, type_pctg, config$target_test_num)
    }
  }
  data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}

constructor.celltype_number <- function(data,config,if_train,sample_seed=NULL){
  if(if_train){
    data <- utils.filter(data)
    train_sample_pctg <- if(purrr::is_null(config$train_sample_pctg)) utils.calculate_sampling_pctg(data, config$target_train_num,config, if_train) else config$train_sample_pctg
    print(str_glue("start sampling train data: {train_sample_pctg} percentage per cell type"))
    if(!purrr::is_null(sample_seed)) print(str_glue("sample_seed is {sample_seed}"))
    data <- utils.sampler(data, sample_pctg = train_sample_pctg, types=config$sample_type, sample_seed=sample_seed)
    
  }else{
    data <- utils.filter(data,filter_gene=FALSE)
    test_sample_pctg <- if(purrr::is_null(config$test_sample_pctg)) utils.calculate_sampling_pctg(data, config$target_test_num,config, if_train) else config$test_sample_pctg
    print(str_glue("start sampling test data: {test_sample_pctg} percentage per cell type"))
    if(!purrr::is_null(sample_seed)) print(str_glue("sample_seed is {sample_seed}"))
    data <- utils.sampler(data,  sample_pctg = test_sample_pctg, types=config$sample_type, sample_seed=sample_seed)
  }
  data <- data[!duplicated(rownames(data)),!duplicated(colnames(data))]
}




constructor.sequencing_depth <- function(data,config,if_train,sample_seed=NULL){
  low_quantile <- config$low_quantile
  high_quantile <- config$high_quantile 
  if(if_train){
    if(config$fixed_train){
      print("in sequencing depth, sample fixed train dataset")
      data <- constructor.base(data,config,if_train,sample_seed)
    }else{
      print("in sequencing depth, sample unfixed train dataset")
      data <- utils.filter(data) %>%
        utils.seqDepthSelector(low_quantile=low_quantile,high_quantile=high_quantile)
    }
  }else{
    if(config$fixed_test){
      print("in sequencing depth, sample fixed test dataset")
      data <- constructor.base(data,config,if_train,sample_seed)
    }else{
      print("in sequencing depth, sample unfixed test dataset")
      data <- utils.filter(data,filter_gene=FALSE) %>%
        utils.seqDepthSelector(low_quantile=low_quantile,high_quantile=high_quantile)
    }
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
