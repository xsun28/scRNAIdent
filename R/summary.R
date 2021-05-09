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
                                   
summarize_experiments <- function(experiment,exp_config){
  require(reshape2)
  switch(experiment,
         simple_accuracy = summary.simple_accuracy(exp_config),
         cell_number = summary.cell_number(exp_config),
         celltype_number = summary.celltype_number(exp_config),
         imbalance_impacts = summary.imbalance_impacts(exp_config),
         sequencing_depth = summary.sequencing_depth(exp_config),
         celltype_structure = summary.celltype_structure(exp_config),
         batch_effects = summary.batch_effects(exp_config),
         inter_diseases = summary.inter_diseases(exp_config),
         celltype_complexity = summary.celltype_complexity(exp_config),
         inter_species = summary.inter_species(exp_config),
         random_noise = summary.random_noise(exp_config),
         inter_protocol = summary.inter_protocol(exp_config),
         stop("Unkown experiments")
  )
  
}


summary.cell_number <- function(exp_config){
  experiment <- "cell_number"
  sample_id_var <- "sample_num"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods",sample_id_var),measure.vars=c("ARI")) 

  metric_results <- dcast(metric_results, train_dataset+test_dataset+assigned+sample_num ~ methods)
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,sample_num,.keep=T) %>% group_map(~{.[1,append(inter_dataset_properties_list,"sample_num")]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","sample_num"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+sample_num ~ methods)
  
  unsup_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","sample_num"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+sample_num ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset","sample_num"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset","sample_num"))
  combined_prop_unsup_other_results <- left_join(unsup_other_results,dataset_prop,by=c("train_dataset","test_dataset","sample_num"))
  summary.output_excel(experiment,exp_config,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results)
}

summary.sequencing_depth <- function(exp_config){
  experiment <- "sequencing_depth"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("dataset","assigned","methods",'quantile'),measure.vars=c("ARI")) %>%
    dcast(dataset+assigned+quantile ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,dataset,quantile,.keep=T) %>% group_map(~{.[1,append(intra_dataset_properties_list,"quantile")]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("dataset","methods","quantile"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics+quantile ~ methods)
  
  unsupervised_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>%
    melt(id.vars = c("dataset","methods","quantile"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics+quantile ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("dataset","quantile"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("dataset","quantile"))
  combined_prop_unsup_other_results <- left_join(unsupervised_other_results,dataset_prop,by=c("dataset","quantile"))
  summary.output_excel(experiment,exp_config,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results)
}


summary.simple_accuracy <- function(exp_config){
  
  experiment <- "simple_accuracy"
  exp_results <- summary.collect_experiment_results(experiment)
  metric_results <- melt(exp_results,id.vars = c("dataset","assigned","methods"),measure.vars=c("ARI")) %>%
    dcast(dataset+assigned ~ methods)
  dataset_prop <- bind_rows(group_by(exp_results,dataset,.keep=T) %>% group_map(~{.[1,intra_dataset_properties_list]},.keep=T))
  
  supervised_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("dataset","methods"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics ~ methods)
  
  unsupervised_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>%
    melt(id.vars = c("dataset","methods"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("dataset"))
  combined_prop_sup_other_results <- left_join(supervised_other_results,dataset_prop,by=c("dataset"))
  combined_prop_unsup_other_results <- left_join(unsupervised_other_results,dataset_prop,by=c("dataset"))
  summary.output_excel(experiment,exp_config,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_sup_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results)
}



summary.inter_diseases <- function(exp_config){
  
  experiment <- "inter_diseases"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods"),measure.vars=c("ARI")) %>%
    dcast(train_dataset+test_dataset+assigned ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,.keep=T) %>% group_map(~{.[1,inter_dataset_properties_list]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics ~ methods)
  
  unsup_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset"))
  combined_prop_unsup_other_results <- left_join(unsup_other_results,dataset_prop,by=c("train_dataset","test_dataset"))
  summary.output_excel(experiment,exp_config,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results)
}

summary.batch_effects <- function(exp_config){
  
  experiment <- "batch_effects"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods",'batch_effects_removed'),measure.vars=c("ARI")) %>%
    dcast(train_dataset+test_dataset+assigned+batch_effects_removed ~ methods)
  
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,batch_effects_removed,.keep=T) %>% group_map(~{.[1,append(inter_dataset_properties_list,c("batch_effects_removed","batch_effects_amount"))]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","batch_effects_removed"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+batch_effects_removed ~ methods)
  
  unsup_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","batch_effects_removed"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+batch_effects_removed ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  combined_prop_unsup_other_results <- left_join(unsup_other_results,dataset_prop,by=c("train_dataset","test_dataset","batch_effects_removed"))
  summary.output_excel(experiment,exp_config,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results)
}

summary.imbalance_impacts <- function(experiment){
  experiment <- "imbalance_impacts"
  sample_id_var <- "imbl_entropy"
  # exp_config <- experiments.parameters[[experiment]]
  exp_results <- summary.collect_experiment_results(experiment)
  # exp_dataset_props <-summary.read_exp_dataset_properties(experiment, exp_config)
  metric_results <- melt(exp_results,id.vars = c("train_dataset","test_dataset","assigned","methods",sample_id_var),measure.vars=c("ARI")) 
  
  metric_results <- dcast(metric_results, train_dataset+test_dataset+assigned+imbl_entropy ~ methods)
  dataset_prop <- bind_rows(group_by(exp_results,train_dataset,test_dataset,imbl_entropy,.keep=T) %>% group_map(~{.[1,append(inter_dataset_properties_list,"imbl_entropy")]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","imbl_entropy"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+imbl_entropy ~ methods)
  
  unsup_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>% 
    melt(id.vars = c("train_dataset","test_dataset","methods","imbl_entropy"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(train_dataset+test_dataset+metrics+imbl_entropy ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("train_dataset","test_dataset","imbl_entropy"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("train_dataset","test_dataset","imbl_entropy"))
  combined_prop_unsup_other_results <- left_join(unsup_other_results,dataset_prop,by=c("train_dataset","test_dataset","imbl_entropy"))
  single_type_results <- summary.collect_experiment_results(experiment,"single_method_results.rds") %>% dplyr::select(-methods)
  if("train_dataset" %in% colnames(single_type_results)){
    single_type_results <- melt(single_type_results, id.vars=c("train_dataset","test_dataset","method","type","type_pctg","supervised","imbl_entropy"))
  }else{
    single_type_results <- melt(single_type_results, id.vars=c("dataset","method","type","type_pctg","supervised","imbl_entropy"))
  }
  
  summary.output_excel(experiment,exp_config,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results,
                       single_type_results=single_type_results)
  
}

summary.celltype_number <- function(){
  experiment <- "celltype_number"
  exp_config <- experiments.parameters[[experiment]]
  type_pctgs <- exp_config$type_pctg
  exp_results <- summary.collect_experiment_results(experiment)
  metric_results <- melt(exp_results,id.vars = c("dataset","assigned","methods"),measure.vars=c("ARI"),type_pctgs) 
  
  metric_results <- dcast(metric_results, dataset+assigned+type_pctgs ~ methods)
  dataset_prop <- bind_rows(group_by(exp_results,dataset,dataset,type_pctgs,.keep=T) %>% group_map(~{.[1,append(inter_dataset_properties_list,"type_pctgs")]},.keep=T))
  other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==TRUE) %>% 
    melt(id.vars = c("dataset","methods","type_pctgs"),measure.vars=c("unlabeled_pctg","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics+type_pctgs ~ methods)
  
  unsup_other_results <- dplyr::filter(exp_results,is.na(assigned),supervised==F) %>% 
    melt(id.vars = c("dataset","methods","type_pctgs"),measure.vars=c("cluster_num","pred_type_max_pctg"),variable.name="metrics") %>%
    dcast(dataset+metrics+sample_pctg ~ methods)
  
  combined_prop_metric_results <- left_join(metric_results,dataset_prop,by=c("dataset","methods","type_pctgs"))
  combined_prop_other_results <- left_join(other_results,dataset_prop,by=c("dataset","methods","type_pctgs"))
  combined_prop_unsup_other_results <- left_join(unsup_other_results,dataset_prop,by=c("dataset","methods","type_pctgs"))
  
  summary.output_excel(experiment,combined_prop_metric_results=combined_prop_metric_results,
                       combined_prop_unlabeled_pctg_results=combined_prop_other_results,
                       combined_prop_unsup_other_results=combined_prop_unsup_other_results)
}




summary.collect_experiment_results <- function(experiment,file_name="results.rds"){
  exp_results_dir <- str_glue({"{result_home}{experiment}"})
  exp_results <- summary.read_results_from_dir(exp_results_dir,file_name) %>% utils.get_methods()
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

summary.read_results_from_dir <- function(current_dir,file_name="results.rds"){
  files <- list.files(current_dir,recursive = F)
  if(file_name %in% files){
    return(read_rds(str_glue("{current_dir}/{file_name}")))
  }else{
    sub_dirs <- list.dirs(current_dir,recursive = F)
    sub_results <- vector("list",length=length(sub_dirs))
    for(i in seq_along(sub_dirs)){
      sub_results[[i]] <- summary.read_results_from_dir(sub_dirs[[i]],file_name)
    }
    return(bind_rows(purrr::map(sub_results,bind_rows)))
  }
}


summary.output_excel <- function(experiment,exp_config,...){
  require(xlsx)
  require(pivottabler)
  results <- list(...)
  combined_prop_metric_results <- results$combined_prop_metric_results
  combined_prop_unlabeled_pctg_results <- results$combined_prop_unlabeled_pctg_results
  combined_prop_unsup_other_results <- results$combined_prop_unsup_other_results
  if(exp_config$fixed_train&&(!exp_config$fixed_test)){
    suffix <- "_train_fixed"
  }else if(exp_config$fixed_test&&(!exp_config$fixed_train)){
    suffix <- "_test_fixed"
  }else{
    suffix <- ""
  }
  
  excel_file_name <- str_glue("{result_home}{experiment}/{experiment}{suffix}_summarized_results.xlsx")
  wb <- xlsx::createWorkbook()
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
  unlabeled_pctg_ws <- createSheet(wb, sheetName = "Supervised methods other metrics(unlabeled pctg...)")
  rows <- createRow(unlabeled_pctg_ws, rowIndex = 1)
  sheetTitle4 <- createCell(rows, colIndex = 1)
  setCellValue(sheetTitle4[[1,1]], "Other metrics(unlabeled pctg...)")
  setCellStyle(sheetTitle4[[1,1]], title_style)
  addDataFrame(combined_prop_unlabeled_pctg_results, sheet = unlabeled_pctg_ws, startRow = 3, startColumn = 1,
               colnamesStyle = colname_style,
               rownamesStyle = rowname_style,
               row.names = FALSE)
  setColumnWidth(sheet = unlabeled_pctg_ws, colIndex = 1:dim(combined_prop_unlabeled_pctg_results)[[2]], colWidth = colWidth)

  ######unsupervised methods other metrics
  unsupervised_other_ws <- createSheet(wb, sheetName = "Unsupervised methods other metrics(cluster number...)")
  rows <- createRow(unsupervised_other_ws, rowIndex = 1)
  sheetTitle5 <- createCell(rows, colIndex = 1)
  setCellValue(sheetTitle5[[1,1]], "Other metrics(cluster number...)")
  setCellStyle(sheetTitle5[[1,1]], title_style)
  addDataFrame(combined_prop_unsup_other_results, sheet = unsupervised_other_ws, startRow = 3, startColumn = 1,
               colnamesStyle = colname_style,
               rownamesStyle = rowname_style,
               row.names = FALSE)
  setColumnWidth(sheet = unsupervised_other_ws, colIndex = 1:dim(combined_prop_unsup_other_results)[[2]], colWidth = colWidth)
  xlsx::saveWorkbook(wb, file = excel_file_name)
  
  if(experiment %in% c("imbalance_impacts")){
    require(openxlsx)
    wb1 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb1,str_glue("Single type performance metrics"))
    single_type_results <- results$single_type_results
    pt <- PivotTable$new()
    pt$addData(single_type_results)
    pt$addRowDataGroups("variable",addTotal=FALSE,header = "metric")
    pt$addRowDataGroups("type",addTotal=FALSE,header = "cell type")
    if("train_dataset" %in% colnames(single_type_results)){
      pt$addRowDataGroups("train_dataset",addTotal=FALSE,header="train dataset")
      pt$addRowDataGroups("test_dataset",addTotal=FALSE,header="test dataset")
    }else{
      pt$addRowDataGroups("dataset",addTotal=FALSE,header="dataset")
    }
    pt$addColumnDataGroups("imbl_entropy", addTotal=FALSE,caption = "Entropy={value}",dataSortOrder="desc")
    pt$defineCalculation(calculationName = "type_pctg", caption = "type_pctg",type="value",valueName="type_pctg")
    pt$addColumnDataGroups("supervised", addTotal=FALSE,caption = "Supervised={value}")
    pt$addColumnDataGroups("method", addTotal=FALSE)
    pt$defineCalculation(calculationName = "value", caption = "metric",type="value",valueName="value",format="%.2f",cellStyleDeclarations=list("xl-value-format"="##0.0"))
    pt$evaluatePivot()
    pt$writeToExcelWorksheet(wb=wb1, topRowNumber = 2, leftMostColumnNumber = 1, wsName =  str_glue("Single type performance metrics"),outputValuesAs="rawValue",
                             applyStyles=TRUE, mapStylesFromCSS=TRUE)
  }
  excel_file_name <- str_glue("{result_home}{experiment}/{experiment}_summarized_single_type_results.xlsx")
  openxlsx::saveWorkbook(wb1, file = excel_file_name,overwrite = T)
}

