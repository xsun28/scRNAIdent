
output.sink <- function(experiment,raw_results,results,exp_config=NULL,...){
  # dataset <- output.dataset_name[[experiment]]
  dataset <- unique(results$dataset)
  train_dataset <- unique(results$train_dataset)
  test_dataset <- unique(results$test_dataset)
  output_dir <- output.generate_output_path(experiment, dataset, train_dataset, test_dataset, exp_config)
  
  switch(experiment,
         simple_accuracy = output.simple_accuracy(raw_results,results,output_dir),
         cell_number = output.cell_number(raw_results,results,output_dir),
         celltype_number = output.celltype_number(raw_results,results,output_dir),
         sequencing_depth = output.sequencing_depth(raw_results,results,output_dir),
         celltype_structure = output.celltype_structure(raw_results,results,output_dir),
         batch_effects = output.batch_effects(raw_results,results,output_dir),
         inter_diseases = output.inter_diseases(raw_results,results,output_dir),
         celltype_complexity = output.celltype_complexity(raw_results,results,output_dir),
         inter_species = output.inter_species(raw_results,results,output_dir),
         random_noise = output.random_noise(raw_results,results,output_dir),
         inter_protocol = output.inter_protocol(raw_results,results,output_dir),
         imbalance_impacts = output.imbalance_impacts(raw_results,results,output_dir,...),
         unknown_types = output.unknown_types(raw_results,results,output_dir,...),
         stop("Unkown experiments")
  )
}

output.simple_accuracy <- function(raw_results,results,output_dir){
  print(str_glue('start writing {experiment} results to {output_dir}'))
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.cell_number <- function(raw_results,results,output_dir){
  print(str_glue('start writing {experiment} results to {output_dir}'))
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.celltype_number <- function(raw_results,results,output_dir){
  print(str_glue('start writing {experiment} results to {output_dir}'))
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.sequencing_depth <- function(raw_results,results, output_dir){
  print(str_glue('start writing {experiment} results to {output_dir}'))
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.celltype_structure <- function(raw_results,results, output_dir){
  print(str_glue('start writing {experiment} results to {output_dir}'))
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.batch_effects <- function(raw_results,results, output_dir){
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.inter_diseases <- function(raw_results,results, output_dir){
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.imbalance_impacts <- function(raw_results,results,output_dir,...){
  single_method_results <- list(...)$single_method_results
  output.write_results(raw_results=raw_results,results=results, single_method_results=single_method_results,output_dir=output_dir)
}

output.unknown_types <- function(raw_results,results,output_dir,...){
  single_method_results <- list(...)$single_method_results
  output.write_results(raw_results=raw_results,results=results, single_method_results=single_method_results,output_dir=output_dir)
}


output.celltype_complexity <- function(raw_results,results, output_dir){
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.inter_species <- function(raw_results,results, output_dir){
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.random_noise <- function(raw_results,results, output_dir){
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.inter_protocol <- function(raw_results,results, output_dir){
  output.write_results(raw_results=raw_results,results=results,output_dir=output_dir)
}

output.write_results <- function(...,output_dir){
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive=T)
  }
  all_results <- list(...)
  result_names <- names(all_results)
  purrr::walk(result_names,~{write_rds(all_results[[.]],str_glue('{output_dir}/{.}.rds'))})
  # write_rds(raw_results,str_glue('{output_dir}/raw_results.rds'))
  # write_rds(results,str_glue('{output_dir}/results.rds'))
  results <- all_results$results
  write_csv(rownames_to_column(bind_rows(results),'method'),str_glue('{output_dir}/results.csv'))
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
output.generate_output_path <- function(experiment, dataset, train_dataset, test_dataset, exp_config){
  if(!purrr::is_null(exp_config$fixed_train)){
    if(exp_config$fixed_train&&!exp_config$fixed_test){
      output_dir <- str_glue("{result_home}{experiment}_train_fixed")
    }
  }
  if(!purrr::is_null(exp_config$fixed_test)){
    if(exp_config$fixed_test&&!exp_config$fixed_train){
      output_dir <- str_glue("{result_home}{experiment}_test_fixed")
    }
  }
  
  if(!purrr::is_null(dataset)){
    output_dir <- str_glue("{output_dir}/{dataset}")
  }else{
    if(purrr::is_null(test_dataset)){
      output_dir <- str_glue("{output_dir}/{train_dataset}")
    }else{
      output_dir <- if(train_dataset==test_dataset) str_glue("{output_dir}/{train_dataset}") else str_glue("{output_dir}/{train_dataset}_{test_dataset}")
    }
  }

  print(str_glue("output dir is {output_dir}"))
  return(str_glue("{output_dir}/"))
}


