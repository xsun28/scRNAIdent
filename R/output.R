
output.sink <- function(experiment,raw_results,results,exp_config=NULL){
  dataset <- output.dataset_name[[experiment]]
  switch(experiment,
         simple_accuracy = output.simple_accuracy(raw_results,results,dataset,exp_config),
         cell_number = output.cell_number(raw_results,results,dataset,exp_config),
         sequencing_depth = output.sequencing_depth(raw_results,results,dataset,exp_config),
         celltype_structure = output.celltype_structure(raw_results,results,dataset,exp_config),
         batch_effects = output.batch_effects(raw_results,results,dataset,exp_config),
         inter_diseases = output.inter_diseases(raw_results,results,dataset,exp_config),
         celltype_complexity = output.celltype_complexity(raw_results,results,dataset,exp_config),
         inter_species = output.inter_species(raw_results,results,dataset,exp_config),
         random_noise = output.random_noise(raw_results,results,dataset,exp_config),
         inter_protocol = output.inter_protocol(raw_results,results,dataset,exp_config),
         stop("Unkown experiments")
  )
}

output.simple_accuracy <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("simple_accuracy",dataset,raw_results,results)
}

output.cell_number <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("cell_number",dataset,raw_results,results)
}

output.sequencing_depth <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("sequencing_depth",dataset,raw_results,results)
}

output.celltype_structure <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("celltype_structure",dataset,raw_results,results)
}

output.batch_effects <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("batch_effects",dataset,raw_results,results)
}

output.inter_diseases <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("inter_diseases",dataset,raw_results,results)
}

output.celltype_complexity <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("celltype_complexity",dataset,raw_results,results)
}

output.inter_species <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("inter_species",dataset,raw_results,results)
}

output.random_noise <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("random_noise",dataset,raw_results,results)
}

output.inter_protocol <- function(raw_results,results,dataset,exp_config=NULL){
  output.write_results("inter_protocol",dataset,raw_results,results)
}

output.write_results <- function(experiment,dataset,raw_results,results){
  print(str_glue('start writing {dataset} {experiment} results'))
  output_dir <- str_glue("{result_home}{experiment}/{dataset}")
  if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive=T)
  }
  write_rds(raw_results,str_glue('{result_home}{experiment}/{dataset}/raw_results.rds'))
  write_rds(results,str_glue('{result_home}{experiment}/{dataset}/results.rds'))
  write_csv(rownames_to_column(bind_rows(results),'method'),str_glue('{result_home}{experiment}/{dataset}/results.csv'))
}

output.dataset.properties_table <- function(t, file_path_name=NULL,prefix=NULL){
  require(xlsx)
  # dataset <- output.dataset_name[[experiment]]
  if(!is_null(prefix)){
    write_rds(t,str_glue('{file_path_name}/{prefix}dataset_properties.rds'))
  }
  else{
    if(!dir.exists(data_home)) dir.create(data_home)
    write_rds(t,str_glue('{data_home}/{file_name}_dataset_properties.rds'))
    write.xlsx(t, file = str_glue('{data_home}/{file_name}_dataset_properties.xlsx'))
  }
}

####
output.cell_number.summary <- function(results,dataset,exp_config){
  
}
