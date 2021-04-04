
driver.drive <- function(experiment=c("simple_accuracy","cell_number", "sequencing_depth","celltype_structure",
                                        "batch_effects","inter_diseases","celltype_complexity","inter_species",
                                        "random_noise","inter_protocol")){
  switch(experiment,
         simple_accuracy = driver.simple_accuracy(),
         cell_number = driver.cell_number(),
         sequencing_depth = driver.sequencing_depth(),
         celltype_structure = driver.celltype_structure(),
         batch_effects = driver.batch_effects(),
         inter_diseases = driver.inter_diseases(),
         celltype_complexity = driver.celltype_complexity(),
         inter_species = driver.inter_species(),
         random_noise = driver.random_noise(),
         inter_protocol = driver.inter_protocol(),
         stop("Unkown experiments")
  )
}

driver.base <- function(experiment,datasets,ouput_paths){ ###for intra-dataset cv base


  for(i in seq_along(datasets)){
    source("R/experiments.R")
    experiment  <<- experiment
    experiments.data[[experiment]] <<- datasets[[i]]
    experiments.assign.data$train_dataset[[experiment]]  <<- datasets[[i]]
    experiments.assign.data$test_dataset[[experiment]]  <<- datasets[[i]]
    output.dataset_name[[experiment]] <<- ouput_paths[[i]]
    runExperiments(experiment)
  }
  summarize_experiments(experiment)
}


driver.simple_accuracy <- function(){
  output_paths <- list("PBMC_AllCells_withLabels","GSE96583_batch1_3_samples","GSE96583_8_Stim_Pats","GSE96583_8_Ctrl_Pats",
                       "Muraro_pancreas_clean", "Segerstolpe_pancreas_clean",
                       "ADASD_AD","ADASD_autism",
                       "midbrain_human","midbrain_mouse"
  )
  datasets <- list("PBMC_AllCells_withLabels.RDS","GSE96583_batch1_3_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS",
                   "Muraro_pancreas_clean.RDS", "Segerstolpe_pancreas_clean.RDS",
                   "ADASD_AD.RDS","ADASD_autism.RDS",
                   "midbrain_human.RDS","midbrain_mouse.RDS"
  )
  driver.base("simple_accuracy", datasets, output_paths)
}


driver.cell_number <- function(){
  output_paths <- list("PBMC_AllCells_withLabels","GSE96583_batch1_3_samples","GSE96583_8_Stim_Pats","GSE96583_8_Ctrl_Pats",
                       "Muraro_pancreas_clean", "Segerstolpe_pancreas_clean",
                       "ADASD_AD","ADASD_autism",
                       "midbrain_human","midbrain_mouse"
                       )
  datasets <- list("PBMC_AllCells_withLabels.RDS","GSE96583_batch1_3_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS",
                   "Muraro_pancreas_clean.RDS", "Segerstolpe_pancreas_clean.RDS",
                   "ADASD_AD.RDS","ADASD_autism.RDS",
                   "midbrain_human.RDS","midbrain_mouse.RDS"
                   )
  driver.base("cell_number", datasets, output_paths)
}


driver.sequencing_depth <- function(){
  output_paths <- list("PBMC_AllCells_withLabels","GSE96583_batch1_3_samples","GSE96583_8_Stim_Pats","GSE96583_8_Ctrl_Pats",
                       "Muraro_pancreas_clean", "Segerstolpe_pancreas_clean",
                       "ADASD_AD","ADASD_autism",
                       "midbrain_human","midbrain_mouse"
  )
  datasets <- list("PBMC_AllCells_withLabels.RDS","GSE96583_batch1_3_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS",
                   "Muraro_pancreas_clean.RDS", "Segerstolpe_pancreas_clean.RDS",
                   "ADASD_AD.RDS","ADASD_autism.RDS",
                   "midbrain_human.RDS","midbrain_mouse.RDS"
  )
  driver.base("sequencing_depth", datasets, output_paths)
}

driver.batch_effects <- function(){
  output_paths <- list("pancreas","PBMC")
  datasets <- list(list(muraro="Muraro_pancreas_clean.RDS",seger="Segerstolpe_pancreas_clean.RDS"),list("PBMC_AllCells_withLabels.RDS","GSE96583_8_Ctrl_Pats.RDS"))
  for(i in seq_along(datasets)){
    source("R/experiments.R")
    experiment  <<- "batch_effects"
    experiments.data[[experiment]] <<- datasets[[i]]
    output.dataset_name[[experiment]] <<- ouput_paths[[i]]
    runExperiments(experiment)
  }
  summarize_experiments(experiment)
}

driver.inter_diseases <- function(){
  output_paths <- list("ADASD","PBMC")
  datasets <- list(list("ADASD_AD.RDS","ADASD_autism.RDS"),
                   list("GSE96583_8_Ctrl_Pats.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_batch1_3_samples.RDS"))
  for(i in seq_along(datasets)){
    source("R/experiments.R")
    experiment  <<- "inter_diseases"
    experiments.data[[experiment]] <<- datasets[[i]]
    output.dataset_name[[experiment]] <<- ouput_paths[[i]]
    runExperiments(experiment)
  }
  summarize_experiments(experiment)
}