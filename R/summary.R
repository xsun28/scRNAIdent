inter_dataset_properties_list <- c(
                                   "train_dataset","train_complexity","train_entropy",
                                   "train_seq_depth_med","train_seq_depth_IQR","train_cell_types",
                                   "train_cell_type_max_pctg","train_sparsity","train_sample_num","train_gene_num",
                                   
                                   "test_dataset","test_complexity","test_entropy","test_seq_depth_med",
                                   "test_seq_depth_IQR","test_cell_types","test_cell_type_max_pctg",
                                   "test_sparsity","test_sample_num","test_gene_num"
                                   )         
                                             
intra_dataset_properties_list <- c("dataset","complexity","entropy","seq_depth_med","seq_depth_IQR",
                                   "cell_types","cell_type_max_pctg","sparsity","sample_num","gene_num")                                             
                                   
summarize_experiments <- function(experiment){
  require(reshape2)
  switch(experiment,
         simple_accuracy = summary.simple_accuracy(),
         cell_number = summary.cell_number(),
         sequencing_depth = summary.sequencing_depth(),
         celltype_structure = summary.celltype_structure(),
         batch_effects = summary.batch_effects(),
         inter_diseases = summary.inter_diseases(),
         celltype_complexity = summary.celltype_complexity(),
         inter_species = summary.inter_species(),
         random_noise = summary.random_noise(),
         inter_protocol = summary.inter_protocol(),
         stop("Unkown experiments")
  )
  
}


summary.cell_number <- function(){
  experiment <- "cell_number"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods",'cell_number'),measure.vars=c("ARI")) %>%
    dcast(train_dataset+test_dataset+assigned+batch_effects_removed ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,batch_effects_removed,.keep=T) %>% group_map(~{.[1,append(inter_dataset_properties_list,c("batch_effects_removed","batch_effects_amount"))]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","batch_effects_removed"),measure.vars=c("unlabeled_pctg","cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+batch_effects_removed ~ methods)
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  summary.output_excel(experiment,combined_prop_metric_results,combined_prop_other_results)
  
}

summary.sequencing_depth <- function(){
  experiment <- "sequencing_depth"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("dataset","assigned","methods",'quantile'),measure.vars=c("ARI")) %>%
    dcast(dataset+assigned+quantile ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,dataset,quantile,.keep=T) %>% group_map(~{.[1,append(intra_dataset_properties_list,"quantile")]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("dataset","methods","quantile"),measure.vars=c("unlabeled_pctg","cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics+quantile ~ methods)
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("dataset","quantile"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("dataset","quantile"))
  summary.output_excel(experiment,combined_prop_metric_results,combined_prop_other_results)
  
}


summary.simple_accuracy <- function(){
  
  experiment <- "simple_accuracy"
  exp_results <- summary.collect_experiment_results(experiment)
  metric_results <- melt(exp_results,id.vars = c("dataset","assigned","methods"),measure.vars=c("ARI")) %>%
    dcast(dataset+assigned ~ methods)
  dataset_prop <- bind_rows(group_by(exp_results,dataset,.keep=T) %>% group_map(~{.[1,intra_dataset_properties_list]},.keep=T))
  
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("dataset","methods"),measure.vars=c("unlabeled_pctg","cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("dataset"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("dataset"))
  summary.output_excel(experiment,combined_prop_metric_results,combined_prop_other_results)
}



summary.inter_diseases <- function(){
  
  experiment <- "inter_diseases"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods"),measure.vars=c("ARI")) %>%
    dcast(train_dataset+test_dataset+assigned ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,.keep=T) %>% group_map(~{.[1,inter_dataset_properties_list]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods"),measure.vars=c("unlabeled_pctg","cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics ~ methods)
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset"))
  summary.output_excel(experiment,combined_prop_metric_results,combined_prop_other_results)
}

summary.batch_effects <- function(){
  
  experiment <- "batch_effects"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods",'batch_effects_removed'),measure.vars=c("ARI")) %>%
    dcast(train_dataset+test_dataset+assigned+batch_effects_removed ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,batch_effects_removed,.keep=T) %>% group_map(~{.[1,append(inter_dataset_properties_list,c("batch_effects_removed","batch_effects_amount"))]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","batch_effects_removed"),measure.vars=c("unlabeled_pctg","cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+batch_effects_removed ~ methods)
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  summary.output_excel(experiment,combined_prop_metric_results,combined_prop_other_results)
}





summary.collect_experiment_results <- function(experiment){
  exp_results_dir <- str_glue({"{result_home}{experiment}"})
  exp_results <- summary.read_results_from_dir(exp_results_dir) %>% utils.get_methods()
  exp_results
}


# summary.read_exp_dataset_properties <- function(experiment,exp_config){
#   exp_results_dir <- str_glue({"{result_home}{experiment}"})
#   exp_dataset_props <- summary.read_dataset_prop_from_dir(exp_results_dir,exp_config)
#   exp_dataset_props
# }
# 
# summary.read_dataset_prop_from_dir <- function(current_dir, exp_config){
#   sub_dirs <- list.dirs(current_dir,recursive = F)
#   if(length(sub_dirs)==0){
#     
#     if(exp_config$cv){
#       dataset_prop <- read_rds(str_glue('{current_dir}/alldataset_properties.rds'))
#       names(dataset_prop)[[1]] <- "dataset"
#       return(dataset_prop)
#     }else{
#       train_dataset_prop <- read_rds(str_glue('{current_dir}/traindataset_properties.rds'))
#       test_dataset_prop <- read_rds(str_glue('{current_dir}/testdataset_properties.rds'))
#       combined_prop <- list(train_dataset=train_dataset_prop$name,train_complexity=train_dataset_prop$complexity,
#                             train_entropy=train_dataset_prop$entropy,train_seq_depth_med=train_dataset_prop$seq_depth_med,
#                             train_seq_depth_IQR=train_dataset_prop$seq_depth_IQR,train_cell_types=train_dataset_prop$cell_types,
#                             train_sparsity=train_dataset_prop$sparsity,train_sample_num=train_dataset_prop$sample_num,
#                             train_gene_num=train_dataset_prop$gene_num,
#                             test_dataset=test_dataset_prop$name,test_complexity=test_dataset_prop$complexity,
#                             test_entropy=test_dataset_prop$entropy,test_seq_depth_med=test_dataset_prop$seq_depth_med,
#                             test_seq_depth_IQR=test_dataset_prop$seq_depth_IQR,test_cell_types=test_dataset_prop$cell_types,
#                             test_sparsity=test_dataset_prop$sparsity,test_sample_num=test_dataset_prop$sample_num,
#                             test_gene_num=test_dataset_prop$gene_num)
#       return(combined_prop)
#     }
#     
#   }else{
#     sub_props <- vector("list",length=length(sub_dirs))
#     for(i in seq_along(sub_dirs)){
#       sub_props[[i]] <- summary.read_dataset_prop_from_dir(sub_dirs[[i]], exp_config)
#     }
#     return(bind_rows(sub_props))
#   }
# }

summary.read_results_from_dir <- function(current_dir){
  files <- list.files(current_dir,recursive = F)
  if("results.rds" %in% files){
    return(read_rds(str_glue("{current_dir}/results.rds")))
  }else{
    sub_dirs <- list.dirs(current_dir,recursive = F)
    sub_results <- vector("list",length=length(sub_dirs))
    for(i in seq_along(sub_dirs)){
      sub_results[[i]] <- summary.read_results_from_dir(sub_dirs[[i]])
    }
    return(bind_rows(sub_results))
  }
}

summary.output_excel <- function(experiment,combined_prop_metric_results,combined_prop_unlabeled_pctg_results){
  require(xlsx)
  excel_file_name <- str_glue("{result_home}{experiment}/{experiment}_summarized_results.xlsx")
  wb <- createWorkbook()
  title_style <- CellStyle(wb) +
    Font(wb, heightInPoints = 16,
         isBold = TRUE)
  rowname_style <- CellStyle(wb) +
    Font(wb, isBold = TRUE,heightInPoints=12)
  colname_style <- CellStyle(wb) +
    Font(wb, isBold = TRUE,heightInPoints=12) +
    Alignment(wrapText = TRUE, horizontal = "ALIGN_CENTER") +
    Border(color = "black",
           position =c("TOP", "BOTTOM"),
           pen =c("BORDER_THIN", "BORDER_THIN"))
  colWidth = 20
  
  #####all_results sheet
  all_results_sheet <- dplyr::filter(combined_prop_metric_results,is.na(assigned))
  all_results_ws <- createSheet(wb, sheetName = "All")
  rows <- createRow(all_results_ws, rowIndex = 1)
  sheetTitle1 <- createCell(rows, colIndex = 1)
  setCellValue(sheetTitle1[[1,1]], "All results(assigned/unassigned combined)")
  setCellStyle(sheetTitle1[[1,1]], title_style)
  addDataFrame(all_results_sheet, sheet = all_results_ws, startRow = 3, startColumn = 1,
               colnamesStyle = colname_style,
               rownamesStyle = rowname_style,
               row.names = FALSE)
  setColumnWidth(sheet = all_results_ws, colIndex = 1:dim(all_results_sheet)[[2]],colWidth = colWidth)
  
  #####assigned results sheet
  assigned_results_sheet <- dplyr::filter(combined_prop_metric_results,assigned==T)
  assigned_results_ws <- createSheet(wb, sheetName = "Assigned")
  rows <- createRow(assigned_results_ws, rowIndex = 1)
  sheetTitle2 <- createCell(rows, colIndex = 1)
  setCellValue(sheetTitle2[[1,1]], "Assigned only")
  setCellStyle(sheetTitle2[[1,1]], title_style)
  addDataFrame(assigned_results_sheet, sheet = assigned_results_ws, startRow = 3, startColumn = 1,
               colnamesStyle = colname_style,
               rownamesStyle = rowname_style,
               row.names = FALSE)
  setColumnWidth(sheet = assigned_results_ws, colIndex = 1:dim(assigned_results_sheet)[[2]], colWidth = colWidth)
  
  #####unassigned results sheet
  unassigned_results_sheet <- dplyr::filter(combined_prop_metric_results,assigned==F)
  unassigned_results_ws <- createSheet(wb, sheetName = "Unassigned")
  rows <- createRow(unassigned_results_ws, rowIndex = 1)
  sheetTitle3 <- createCell(rows, colIndex = 1)
  setCellValue(sheetTitle3[[1,1]], "Unassigned only")
  setCellStyle(sheetTitle3[[1,1]], title_style)
  addDataFrame(unassigned_results_sheet, sheet = unassigned_results_ws, startRow = 3, startColumn = 1,
               colnamesStyle = colname_style,
               rownamesStyle = rowname_style,
               row.names = FALSE)
  setColumnWidth(sheet = unassigned_results_ws, colIndex = 1:dim(unassigned_results_sheet)[[2]],colWidth = colWidth)
  
  ######unlabeled pctg sheet
  unlabeled_pctg_ws <- createSheet(wb, sheetName = "Other metrics(unlabeled pctg, cluster num...)")
  rows <- createRow(unlabeled_pctg_ws, rowIndex = 1)
  sheetTitle4 <- createCell(rows, colIndex = 1)
  setCellValue(sheetTitle4[[1,1]], "Other metrics(unlabeled pctg, cluster num...)")
  setCellStyle(sheetTitle4[[1,1]], title_style)
  addDataFrame(combined_prop_unlabeled_pctg_results, sheet = unlabeled_pctg_ws, startRow = 3, startColumn = 1,
               colnamesStyle = colname_style,
               rownamesStyle = rowname_style,
               row.names = FALSE)
  setColumnWidth(sheet = unlabeled_pctg_ws, colIndex = 1:dim(combined_prop_unlabeled_pctg_results)[[2]], colWidth = colWidth)
  saveWorkbook(wb, file = excel_file_name)
  
}

