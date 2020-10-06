source("R/utils.R")
source("R/config.R")
####Wrapper
constructor.data_constructor <- function(data,experiment,if_train=TRUE){
  require(SingleCellExperiment)
  stopifnot(is(data,"SingleCellExperiment"))
  switch(experiment,
    simple_accuracy = constructor.simple_accuracy(experiment,data,if_train),
    cell_number = constructor.cell_number(experiment,data,if_train),
    sequencing_depth = constructor.equencing_depth(experiment,data,if_train),
    cell_types = constructor.cell_types(experiment,data,if_train),
    batch_effects = constructor.batch_effects(experiment,data,if_train),
    stop("Can't constructor dataset for unkown experiments")
  )
}

constructor.simple_accuracy <- function(experiment,data, if_train){
  if(if_train){
    data <- utils.filter(data) %>%
              utils.sampler(sample_num=experiments.parameters[[experiment]]$sample_num,types=unique(colData(data)$label))
    return(data)
  }else{
    return(utils.filter(data,filter_gene=FALSE))
  }
}