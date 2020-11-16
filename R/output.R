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
  print('start writing simple accuracy results')
  write_rds(raw_results,paste(result_home,'simple_accuracy_raw_results.rds',sep = ''))
  write_rds(results,paste(result_home,'simple_accuracy_results.rds',sep = ''))
  write_csv(rownames_to_column(bind_rows(results),'method'),paste(result_home,'simple_accuracy_results.csv',sep = ''))
}

output.cell_number <- function(raw_results,results){
  print('start writing cell number results')
  write_rds(raw_results,paste(result_home,'cell_number_raw_results.rds',sep = ''))
  write_rds(results,paste(result_home,'cell_number_results.rds',sep = ''))
  write_csv(rownames_to_column(bind_rows(results),'method'),paste(result_home,'cell_number_results.csv',sep = ''))
}

output.sequencing_depth <- function(raw_results,results){
  print('start writing sequencing depth results')
  write_rds(raw_results,paste(result_home,'sequencing_depth_raw_results.rds',sep = ''))
  write_rds(results,paste(result_home,'sequencing_depth_results.rds',sep = ''))
  write_csv(rownames_to_column(bind_rows(results),'method'),paste(result_home,'sequencing_depth_results.csv',sep = ''))
}

output.batch_effects <- function(raw_results,results){
  print('start writing batch effects results')
  write_rds(raw_results,paste(result_home,'batch_effects_raw_results.rds',sep = ''))
  write_rds(results,paste(result_home,'batch_effects_results.rds',sep = ''))
  write_csv(rownames_to_column(results,'method'),paste(result_home,'batch_effects_results.csv',sep = ''))
}