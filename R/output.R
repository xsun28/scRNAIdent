source('R/config.R')
output.sink <- function(experiment,results){
  switch(experiment,
         simple_accuracy = output.simple_accuracy(results),
         cell_number = output.cell_number(results),
         sequencing_depth = output.sequencing_depth(results),
         cell_types = output.cell_types(results),
         batch_effects = output.batch_effects(results),
         stop("Unkown experiments")
  )
}

output.simple_accuracy <- function(results){
  print('start writing simple accuracy results')
  write_rds(results,paste(result_home,'simple_accuracy_results.rds',sep = ''))
  write_csv(rownames_to_column(bind_rows(results),'method'),paste(result_home,'simple_accuracy_results.csv',sep = ''))
}

output.cell_number <- function(results){
  print('start cell number results')
  write_rds(results,paste(result_home,'cell_number_results.rds',sep = ''))
  write_csv(rownames_to_column(bind_rows(results),'method'),paste(result_home,'cell_number_results.csv',sep = ''))
}

output.sequencing_depth <- function(results){
  print('start sequencing depth results')
  write_rds(results,paste(result_home,'sequencing_depth_results.rds',sep = ''))
  write_csv(rownames_to_column(bind_rows(results),'method'),paste(result_home,'sequencing_depth_results.csv',sep = ''))
}