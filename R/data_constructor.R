source("R/utils.R")
source("R/config.R")
####Wrapper
constructor.data_constructor <- function(data,config,experiment,if_train=TRUE){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  switch(experiment,
    simple_accuracy = constructor.simple_accuracy(data,config,if_train),
    cell_number = constructor.cell_number(data,config,if_train),
    sequencing_depth = constructor.sequencing_depth(data,config,if_train),
    cell_types = constructor.cell_types(data,config,if_train),
    batch_effects = constructor.batch_effects(data,config,if_train),
    stop("Can't constructor dataset for unkown experiments")
  )
}

constructor.simple_accuracy <- function(data,config, if_train){
  if(if_train){
    data <- utils.filter(data) %>%
              utils.sampler(sample_num=config$sample_num,types=unique(colData(data)$label))
    return(data)
  }else{
    return(utils.filter(data,filter_gene=FALSE))
  }
}

constructor.cell_number <- function(data,config, if_train){
  constructor.simple_accuracy(data,config,if_train)
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
  data
}

constructor.batch_effects <- function(data,config,if_train){
  sces <- utils.combine_SCEdatasets(data,if_combined=FALSE)
  if(config$remove_batch){
    sces_no_be <- utils.remove_batch_effects(sces)
  }
}