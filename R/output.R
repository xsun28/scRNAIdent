source('R/config.R')
output.sink <- function(experiment,raw_results,results){
  switch(experiment,
         simple_accuracy = output.simple_accuracy(raw_results,results),
         cell_number = output.cell_number(raw_results,results),
         sequencing_depth = output.sequencing_depth(raw_results,results),
         cell_types = output.cell_types(raw_results,results),
         batch_effects = output.batch_effects(raw_results,results),
         stop("Unkown experiments")
  )
}

output.simple_accuracy <- function(raw_results,results){
  dataset <- experiments.data$simple_accuracy
  output.write_results("simple_accuracy",dataset,raw_results,results)
}

output.cell_number <- function(raw_results,results){
  dataset <- experiments.data$cell_number
  output.write_results("cell_number",dataset,raw_results,results)
}

output.sequencing_depth <- function(raw_results,results){
  dataset <- experiments.data$sequencing_depth
  output.write_results("sequencing_depth",dataset,raw_results,results)
}

output.batch_effects <- function(raw_results,results){
  dataset <- experiments.data$batch_effects
  output.write_results("batch_effects",dataset,raw_results,results)
}

output.write_results <- function(experiment,dataset,raw_results,results){
  print(str_glue('start writing {dataset} {experiment} results'))
  write_rds(raw_results,str_glue('{result_home}{experiment}_{dataset}_raw_results.rds'))
  write_rds(results,str_glue('{result_home}{experiment}_{dataset}_results.rds'))
  write_csv(rownames_to_column(bind_rows(results),'method'),str_glue('{result_home}{experiment}_{dataset}_results.csv'))
}