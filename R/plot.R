
plot.plot <- function(experiment,results,raw_results){
  dataset <- output.dataset_name[[experiment]]
  
  if(purrr::is_null(results)){
    results <- read_rds(str_glue("{result_home}/{experiment}/{dataset}/results.rds"))
  }
  if(purrr::is_null(raw_results)){
    raw_results <- read_rds(str_glue("{result_home}/{experiment}/{dataset}/raw_results.rds"))
  }
  switch(experiment,
    simple_accuracy = plot.simple_accuracy(results,raw_results,dataset),
    cell_number = plot.cell_number(results,raw_results,dataset),
    sequencing_depth = plot.sequencing_depth(results,raw_results,dataset),
    celltype_structure = plot.celltype_structure(results,raw_results,dataset),
    batch_effects = plot.batch_effects(results,raw_results,dataset),
    inter_diseases = plot.inter_diseases(results,raw_results,dataset),
    celltype_complexity = plot.celltype_complexity(results,raw_results,dataset),
    inter_species = plot.inter_species(results,raw_results,dataset),
    random_noise = plot.random_noise(results,raw_results,dataset),
    inter_protocol = plot.inter_protocol(results,raw_results,dataset),
    stop("Unkown experiments")
  )
}

plot.simple_accuracy <- function(results,raw_results,dataset){
  require(data.table)
  fig_path <- str_glue("{result_home}{experiment}/{dataset}/")
  if(!dir.exists(fig_path)){
    dir.create(fig_path,recursive=T)
  }
  # figure_all_name <- str_glue("simple_accuracy_{dataset}_all")
  # figure_assigned_name <- str_glue("simple_accuracy_{dataset}_assigned")
  # figure_unassigned_name <- str_glue("simple_accuracy_{dataset}_unassigned")
  
  all_results <- bind_rows(results[grepl("^all_.*",names(results))])
  all_results[,'methods'] <- rownames(all_results)
  results_assigned <- bind_rows(results[grepl("^assigned_.*",names(results))])
  results_assigned[,'methods'] <- rownames(results_assigned)
  results_unassigned <- bind_rows(results[grepl("^unassigned_.*",names(results))])
  results_unassigned[,'methods'] <- rownames(results_unassigned)

  all_results_table <- data.table(melt(all_results,
                                          id.vars = c('methods','supervised'),    
                                          measure.vars=c('ARI','AMI','FMI')))
  all_results_table$methods <- factor(all_results_table$methods)
  all_results_table$label <- round(all_results_table$value,2)
  plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                      facet_var="supervised",dodged=T,facet_wrap=F,width=10,height=7,nrow=1)
  plot_params$ylabel <- "ARI"
  plot.bar_plot(all_results_table[,all_results_table[variable=="ARI"]],plot_params,fig_path,"all_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(all_results_table[,all_results_table[variable=="AMI"]],plot_params,fig_path,str_glue("all_AMI.pdf"))
  plot_params$ylabel <- "FMI"
  plot.bar_plot(all_results_table[,all_results_table[variable=="FMI"]],plot_params,fig_path,str_glue("all_FMI.pdf"))
  
  results_assign_table <- data.table(melt(results_assigned,
                                   id.vars = c('methods','supervised'),    
                                   measure.vars=c('ARI','AMI','FMI')))
  results_assign_table$methods <- factor(results_assign_table$methods)
  results_assign_table$label <- round(results_assign_table$value,2)
  plot.bar_plot(results_assign_table[,results_assign_table[variable=="ARI"]],plot_params,fig_path,"assigned_ARI.pdf")
  plot.bar_plot(results_assign_table[,results_assign_table[variable=="AMI"]],plot_params,fig_path,"assigned_AMI.pdf")
  plot.bar_plot(results_assign_table[,results_assign_table[variable=="FMI"]],plot_params,fig_path,"assigned_FMI.pdf")
  
  
  results_unassign_table <- data.table(melt(results_unassigned,
                                          id.vars = c('methods','supervised'),    
                                          measure.vars=c('ARI','AMI','FMI')))
  results_unassign_table$methods <- factor(results_unassign_table$methods)
  results_unassign_table$label <- round(results_unassign_table$value,2)
  plot.bar_plot(results_unassign_table[,results_unassign_table[variable=="ARI"]],plot_params,fig_path,"unassigned_ARI.pdf")
  plot.bar_plot(results_unassign_table[,results_unassign_table[variable=="AMI"]],plot_params,fig_path,"unassigned_AMI.pdf")
  plot.bar_plot(results_unassign_table[,results_unassign_table[variable=="FMI"]],plot_params,fig_path,"unassigned_FMI.pdf")
  
  #####unlabeled pctg 
  # figure_name <- str_glue("unlabeled_pctg.pdf")
  results_unlabeled_pctg <- results[["all_assign_results"]]['unlabeled_pctg']
  results_unlabeled_pctg[,'methods'] <- factor(rownames(results_unlabeled_pctg))
  results_unlabeled_pctg[,'label'] <- round(results_unlabeled_pctg$unlabeled_pctg,2)
  plot_params <- list(x="methods",y="unlabeled_pctg",label="label",xlabel="methods",ylabel="unlabeled_pctg",
                      fill="methods",dodged=F,facet_wrap=F,width=6,height=4)
  plot.bar_plot(results_unlabeled_pctg,plot_params,fig_path,"unlabeled_pctg.pdf")
  
 ####sankey plot
 # figure_name <- str_glue("simple_accuracy_{dataset}_sankey.pdf")
 methods <- colnames(dplyr::select(raw_results,-label))

 raw_results1 <- gather(raw_results,"methods","pred",-label) %>%
                   group_by(methods,label,pred) %>%
                     summarize(freq=n()) %>%
                       filter(freq>0)
 
 sankey_plot_params <- list(y="freq", axis1="label",axis2="pred",fill="label",label="label",xlabel="methods",
                            scale_x_discrete_limits="c(\"label\",\"value\")",facet_var="methods",facet_wrap=T,width=10,height=7,nrow=1)
 
 plot.sankey_plot(raw_results1,sankey_plot_params,fig_path,"sankey.pdf")
}

plot.cell_number <- function(results,raw_results,dataset){
  # dataset <- experiments.data$cell_number
  fig_path <- str_glue("{result_home}{experiment}/{dataset}/")
  if(!dir.exists(fig_path)){
    dir.create(fig_path,recursive=T)
  }
  
  cell_numbers <- experiments.parameters[[experiment]]$sample_num
  sampling_by_pctg <- purrr::is_null(cell_numbers)||length(cell_numbers)==0
  
  results <- utils.get_methods(results)
  if(sampling_by_pctg){
    all_results <- dplyr::select(dplyr::filter(results,is.na(assigned)),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,sample_pctg,supervised))
    
    plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                        facet_wrap=T,dodged=T,facet_var='sample_pctg',width=10,height=7,nrow=3)
  }else{
    all_results <- dplyr::select(dplyr::filter(results,is.na(assigned)),-c(unlabeled_pctg,assigned)) %>%
                    gather("metric","value",-c(methods,sample_num,supervised))
    
    plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                        facet_wrap=T,dodged=T,facet_var='sample_num',width=10,height=7,nrow=3)
  }
  all_results$label <- round(all_results$value,2)
  # plot_params <-
  # list(x="sample_num",y="value",group="methods",line_color="methods",point_color="methods",
  # facet_wrap=T,facet_var="metric",width=10,height=10)
  # plot.line_plot(all_results[all_results$metric=="ARI",],plot_params,result_home,figure_all_name)
  

  plot_params$ylabel <- "ARI"
  plot.bar_plot(all_results[all_results$metric=="ARI",],plot_params,fig_path,"all_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(all_results[all_results$metric=="AMI",],plot_params,fig_path,"all_AMI.pdf")
  plot_params$ylabel <- "FMI"
  plot.bar_plot(all_results[all_results$metric=="FMI",],plot_params,fig_path,"all_FMI.pdf")
  
  if(sampling_by_pctg){
    results_assigned <- dplyr::select(dplyr::filter(results,assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,sample_pctg,supervised))
  }
  else{
    results_assigned <- dplyr::select(dplyr::filter(results,assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,sample_num,supervised))
  }
  results_assigned$label <- round(results_assigned$value,2)
  plot_params$ylabel <- "ARI"
  plot.bar_plot(results_assigned[results_assigned$metric=="ARI",],plot_params,fig_path,"assigned_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(results_assigned[results_assigned$metric=="AMI",],plot_params,fig_path,"assigned_AMI.pdf")
  plot_params$ylabel <- "FMI"
  plot.bar_plot(results_assigned[results_assigned$metric=="FMI",],plot_params,fig_path,"assigned_FMI.pdf")
  # plot.line_plot(results_assigned,plot_params,result_home,figure_assigned_name)
  
  if(sampling_by_pctg){
    results_unassigned <- dplyr::select(dplyr::filter(results,!assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,sample_pctg,supervised))
  }
  else{
    results_unassigned <- dplyr::select(dplyr::filter(results,!assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,sample_num,supervised))
  }
  results_unassigned$label <- round(results_unassigned$value,2)
  plot_params$ylabel <- "ARI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="ARI",],plot_params,fig_path,"unassigned_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="AMI",],plot_params,fig_path,"unassigned_AMI.pdf")
  plot_params$ylabel <- "FMI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="FMI",],plot_params,fig_path,"unassigned_FMI.pdf")
  
  # plot.line_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
  
  ######
  # figure_name <- str_glue("cell_number_{dataset}_unlabeled_pctg.pdf")
  if(sampling_by_pctg){
    results_unlabeled_pctg <- dplyr::filter(dplyr::select(results[results$supervised,],methods,sample_pctg,unlabeled_pctg),!is.na(unlabeled_pctg))
    
    plot_params <- list(x="sample_pctg",y="unlabeled_pctg",group="methods",line_color="methods",point_color="methods",
                        facet_wrap=F,width=2*length(unique(results_unlabeled_pctg$methods)),
                        height=2*length(unique(results_unlabeled_pctg$methods)))
  }
  else{
    results_unlabeled_pctg <- dplyr::filter(dplyr::select(results[results$supervised,],methods,sample_num,unlabeled_pctg),!is.na(unlabeled_pctg))
    
    plot_params <- list(x="sample_num",y="unlabeled_pctg",group="methods",line_color="methods",point_color="methods",
                        facet_wrap=F,width=2*length(unique(results_unlabeled_pctg$methods)),
                        height=2*length(unique(results_unlabeled_pctg$methods)))
  }
  plot.line_plot(results_unlabeled_pctg,plot_params,fig_path,"unlabeled_pctg.pdf")
  
  ###sankey_plot
  sample_plan <- if(sampling_by_pctg) "sample_pctg" else "sample_num"
  uniq_sample_plan <- unique(raw_results[[sample_plan]])
  for(plan in uniq_sample_plan){
    raw_results1 <- raw_results[raw_results[[sample_plan]]==plan,]
    methods <- colnames(dplyr::select(raw_results1,-c(label,sample_plan)))
    
    raw_results1 <- gather(raw_results1,"methods","pred",-c(label,sample_plan)) %>%
      group_by(methods,label,pred) %>%
      summarize(freq=n()) %>%
      filter(freq>0)
    
    sankey_plot_params <- list(y="freq", axis1="label",axis2="pred",fill="label",label="label",xlabel="methods",
                               scale_x_discrete_limits="c(\"label\",\"value\")",facet_var="methods",facet_wrap=T,width=10,height=7,nrow=1)
    
    plot.sankey_plot(raw_results1,sankey_plot_params,fig_path,str_glue("sankey_{sample_plan}_{plan}.pdf"))
  }
  
}


plot.celltype_structure <- function(results,raw_results,dataset){
  # dataset <- experiments.data$celltype_structure
  fig_path <- str_glue("{result_home}{experiment}/{dataset}/")
  if(!dir.exists(fig_path)){
    dir.create(fig_path,recursive=T)
  }
  # figure_all_name <- str_glue("celtype_structure_{dataset}_all")
  # figure_assigned_name <- str_glue("celltype_structure_{dataset}_assigned")
  # figure_unassigned_name <- str_glue("celltype_structure_{dataset}_unassigned")
  
  results <- utils.get_methods(results)
  all_results <- dplyr::select(dplyr::filter(results,is.na(assigned)),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,level,supervised))
  all_results$level <- factor(all_results$level)
  all_results$label <- round(all_results$value,2)
  
  plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                      facet_wrap=T,dodged=T,facet_var='level',width=10,height=7,nrow=length(unique(all_results$level)))
  plot_params$ylabel <- "wRI"
  plot.bar_plot(all_results[all_results$metric=="wRI",],plot_params,fig_path,"all_wRI.pdf")
  plot_params$ylabel <- "wNMI"
  plot.bar_plot(all_results[all_results$metric=="wNMI",],plot_params,fig_path,"all_wNMI.pdf")
  
  
  # plot_params <- list(x="level",y="value",group="methods",line_color="methods",point_color="methods",
  #                     facet_wrap=T,facet_var="metric",width=10,height=10)
  # plot.line_plot(all_results,plot_params,result_home,figure_all_name)
  
  results_assigned <- dplyr::select(dplyr::filter(results,assigned),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,level,supervised))
  results_assigned$level <- factor(results_assigned$level)
  results_assigned$label <- round(results_assigned$value,2)
  # plot.line_plot(results_assigned,plot_params,result_home,figure_assigned_name)
  plot_params$ylabel <- "wRI"
  plot.bar_plot(results_assigned[results_assigned$metric=="wRI",],plot_params,fig_path,"assigned_wRI.pdf")
  plot_params$ylabel <- "wNMI"
  plot.bar_plot(results_assigned[results_assigned$metric=="wNMI",],plot_params,fig_path,"assigned_wNMI.pdf")
  
  
  results_unassigned <- dplyr::select(dplyr::filter(results,!assigned),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,level,supervised))
  results_unassigned$level <- factor(results_unassigned$level)
  results_unassigned$label <- round(results_unassigned$value,2)
  plot_params$ylabel <- "wRI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="wRI",],plot_params,fig_path,"unassigned_wRI.pdf")
  plot_params$ylabel <- "wNMI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="wNMI",],plot_params,fig_path,"unassigned_wNMI.pdf")
  
  # plot.line_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
  
  ######
  # figure_name <- str_glue("celltype_structure_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- dplyr::filter(dplyr::select(results[results$supervised,],methods,level,unlabeled_pctg),!is.na(unlabeled_pctg))
  results_unlabeled_pctg$level <- factor(results_unlabeled_pctg$level)
  plot_params <- list(x="level",y="unlabeled_pctg",group="methods",line_color="methods",point_color="methods",
                      facet_wrap=F,width=2*length(unique(results_unlabeled_pctg$methods)),
                      height=2*length(unique(results_unlabeled_pctg$methods)))
  plot.line_plot(results_unlabeled_pctg,plot_params,fig_path,"unlabeled_pctg.pdf")
}


plot.sequencing_depth <- function(results,raw_results,dataset){
  # dataset <- experiments.data$sequencing_depth
  
  results <- utils.get_methods(results)
  results$quantile <- factor(results$quantile)
  
  #######
  fig_path <- str_glue("{result_home}{experiment}/{dataset}/")
  if(!dir.exists(fig_path)){
    dir.create(fig_path,recursive=T)
  }
  
  # figure_all_name <- str_glue("sequencing_depth_{dataset}_all")
  # figure_assigned_name <- str_glue("sequencing_depth_{dataset}_assigned")
  # figure_unassigned_name <- str_glue("sequencing_depth_{dataset}_unassigned")
  
  
  all_results <- dplyr::select(dplyr::filter(results,is.na(assigned)),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,quantile,supervised))
  all_results$label <- round(all_results$value,2)
  plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                      facet_wrap=T,dodged=T,facet_var="quantile",width=10,height=7,nrow=length(unique(all_results$quantile)))
  plot_params$ylabel <- "ARI"
  plot.bar_plot(all_results[all_results$metric=="ARI",],plot_params,fig_path,"all_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(all_results[all_results$metric=="AMI",],plot_params,fig_path,"all_AMI.pdf")
  plot_params$ylabel <- "FMI"
  plot.bar_plot(all_results[all_results$metric=="FMI",],plot_params,fig_path,"all_FMI.pdf")
  
  # plot.heatmap_plot(all_results,"quantile","methods",result_home,figure_all_name)
  
  results_assigned <- dplyr::select(dplyr::filter(results,assigned),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,quantile,supervised))
  results_assigned$label <- round(results_assigned$value,2)
  plot_params$ylabel <- "ARI"
  plot.bar_plot(results_assigned[results_assigned$metric=="ARI",],plot_params,fig_path,"assigned_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(results_assigned[results_assigned$metric=="AMI",],plot_params,fig_path,"assigned_AMI.pdf")
  plot_params$ylabel <- "FMI"
  plot.bar_plot(results_assigned[results_assigned$metric=="FMI",],plot_params,fig_path,"assigned_FMI.pdf")
  
  
  results_unassigned <- dplyr::select(dplyr::filter(results,!assigned),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,quantile,supervised))
  results_unassigned$label <- round(results_unassigned$value,2)
  plot_params$ylabel <- "ARI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="ARI",],plot_params,fig_path,"unassigned_ARI.pdf")
  plot_params$ylabel <- "AMI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="AMI",],plot_params,fig_path,"unassigned_AMI.pdf")
  plot_params$ylabel <- "FMI"
  plot.bar_plot(results_unassigned[results_unassigned$metric=="FMI",],plot_params,fig_path,"unassigned_FMI.pdf")
  
  #######
  # figure_name <- str_glue("sequencing_depth_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- dplyr::filter(dplyr::select(results[results$supervised,],methods,quantile,unlabeled_pctg),!is.na(unlabeled_pctg))
  results_unlabeled_pctg[,'label'] <- round(results_unlabeled_pctg$unlabeled_pctg,2)
  plot_params <- list(x="methods",y="unlabeled_pctg",label="label",xlabel="methods",ylabel="unlabeled_pctg",
                      fill="quantile",dodged=T,facet_wrap=F,width=6,height=4)
  plot.bar_plot(results_unlabeled_pctg,plot_params,fig_path,"unlabeled_pctg.pdf")
  # ggplot(results_unlabeled_pctg, aes(quantile,methods, fill=unlabeled_pctg)) + 
  #   geom_tile() +
  #   scale_fill_distiller(palette = "RdPu") +
  #   theme_ipsum()
  # ggsave(
  #   figure_name,
  #   plot = last_plot(),
  #   device = 'pdf',
  #   path = result_home,
  #   width = 2*length(unique(results_unlabeled_pctg$methods)),
  #   height = 2*length(unique(results_unlabeled_pctg$methods)),
  #   units = "in"
  # )
####sankey plot
  
  uniq_quantile <- unique(raw_results[["quantile"]])
  for(quantile in uniq_quantile){
    raw_results1 <- raw_results[raw_results$quantile==quantile,]
    methods <- colnames(dplyr::select(raw_results1,-c("label","quantile")))
    
    raw_results1 <- gather(raw_results1,"methods","pred",-c(label,"quantile")) %>%
      group_by(methods,label,pred) %>%
      summarize(freq=n()) %>%
      filter(freq>0)
    
    sankey_plot_params <- list(y="freq", axis1="label",axis2="pred",fill="label",label="label",xlabel="methods",
                               scale_x_discrete_limits="c(\"label\",\"value\")",facet_var="methods",facet_wrap=T,width=10,height=7,nrow=1)
    
    plot.sankey_plot(raw_results1,sankey_plot_params,fig_path,str_glue("sankey_quantile_{quantile}.pdf"))
  }
}

plot.batch_effects <- function(results,raw_results,dataset){
  # dataset <- experiments.data$batch_effects
  require(data.table)
  results <- utils.get_methods(results)
  # figure_all_name <- str_glue("batch_effects_{dataset}_all")
  # figure_assigned_name <- str_glue("batch_effects_{dataset}_assigned")
  # figure_unassigned_name <- str_glue("batch_effects_{dataset}_unassigned")

  for(grouped_results in dplyr::group_split(results,train_dataset,test_dataset)){
    train_dataset <- unique(grouped_results$train_dataset)
    test_dataset <- unique(grouped_results$test_dataset)
    fig_path <- str_glue("{result_home}{experiment}/{dataset}/{train_dataset}_{test_dataset}/")
    if(!dir.exists(fig_path)){
      dir.create(fig_path,recursive=T)
    }
    print(str_glue("Generating figures for batch effect pairs {fig_path}"))
    
    all_results <- dplyr::select(dplyr::filter(grouped_results,is.na(assigned)),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,batch_effects_removed,supervised))
    all_results$value <- as.double(all_results$value)
    all_results$label <- round(all_results$value,2)
    plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                        facet_wrap=T,dodged=T,facet_var="batch_effects_removed",width=10,height=7,nrow=2)
    plot_params$ylabel <- "ARI"
    plot.bar_plot(all_results[all_results$metric=="ARI",],plot_params,fig_path,"all_ARI.pdf")
    plot_params$ylabel <- "AMI"
    plot.bar_plot(all_results[all_results$metric=="AMI",],plot_params,fig_path,"all_AMI.pdf")
    plot_params$ylabel <- "FMI"
    plot.bar_plot(all_results[all_results$metric=="FMI",],plot_params,fig_path,"all_FMI.pdf")
    
    # plot.heatmap_plot(all_results,"batch_effects_removed","methods",result_home,figure_all_name)
    
    results_assigned <- dplyr::select(dplyr::filter(grouped_results,assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,batch_effects_removed,supervised))
    results_assigned$value <- as.double(results_assigned$value)
    results_assigned$label <- round(results_assigned$value,2)
    plot_params$ylabel <- "ARI"
    plot.bar_plot(results_assigned[results_assigned$metric=="ARI",],plot_params,fig_path,"assigned_ARI.pdf")
    plot_params$ylabel <- "AMI"
    plot.bar_plot(results_assigned[results_assigned$metric=="AMI",],plot_params,fig_path,"assigned_AMI.pdf")
    plot_params$ylabel <- "FMI"
    plot.bar_plot(results_assigned[results_assigned$metric=="FMI",],plot_params,fig_path,"assigned_FMI.pdf")
    
    
    results_unassigned <- dplyr::select(dplyr::filter(grouped_results,!assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,batch_effects_removed,supervised))
    results_unassigned$value <- as.double(results_unassigned$value)
    results_unassigned$label <- round(results_unassigned$value,2)
    plot_params$ylabel <- "ARI"
    plot.bar_plot(results_unassigned[results_unassigned$metric=="ARI",],plot_params,fig_path,"unassigned_ARI.pdf")
    plot_params$ylabel <- "AMI"
    plot.bar_plot(results_unassigned[results_unassigned$metric=="AMI",],plot_params,fig_path,"unassigned_AMI.pdf")
    plot_params$ylabel <- "FMI"
    plot.bar_plot(results_unassigned[results_unassigned$metric=="FMI",],plot_params,fig_path,"unassigned_FMI.pdf")
    
    
    # figure_name <- str_glue("batch_effects_{dataset}_unlabeled_pctg.pdf")
    results_unlabeled_pctg <- dplyr::filter(dplyr::select(grouped_results[grouped_results$supervised,],methods,batch_effects_removed,unlabeled_pctg),!is.na(unlabeled_pctg))
    results_unlabeled_pctg$label <- round(results_unlabeled_pctg$unlabeled_pctg,2)
    plot_params <- list(x="methods",y="unlabeled_pctg",fill="batch_effects_removed",label="label",xlabel="methods",ylabel="unlabeled_pctg",
                        facet_wrap=F,dodged=T,width=10,height=7)
    plot.bar_plot(results_unlabeled_pctg,plot_params,fig_path,"unlabeled_pctg.pdf")
    # ggplot(results_unlabeled_pctg, aes(batch_effects_removed,methods, fill=unlabeled_pctg)) + 
    #   geom_tile() +
    #   scale_fill_distiller(palette = "RdPu") +
    #   theme_ipsum()
    # ggsave(
    #   figure_name,
    #   plot = last_plot(),
    #   device = 'pdf',
    #   path = result_home,
    #   width = 2*length(unique(results_unlabeled_pctg$methods)),
    #   height = 2*length(unique(results_unlabeled_pctg$methods)),
    #   units = "in"
    # )
  }
  for(grouped_raw_results in dplyr::group_split(raw_results,train_dataset,test_dataset)){
    ###sankey plot
    train_dataset <- unique(grouped_raw_results$train_dataset)
    test_dataset <- unique(grouped_raw_results$test_dataset)
    batch_effects_vals <-  unique(grouped_raw_results$batch_effects)
    for(unique_val in batch_effects_vals){
      fig_path <- str_glue("{result_home}{experiment}/{dataset}/{train_dataset}_{test_dataset}/")
      
      grouped_raw_results1 <- grouped_raw_results[grouped_raw_results$batch_effects==unique_val,]
      methods <- colnames(dplyr::select(grouped_raw_results1,-c("label","batch_effects","train_dataset","test_dataset")))
      
      grouped_raw_results1 <- gather(grouped_raw_results1,"methods","pred",-c(label,"batch_effects","train_dataset","test_dataset")) %>%
        group_by(methods,label,pred) %>%
        summarize(freq=n()) %>%
        filter(freq>0)
      
      sankey_plot_params <- list(y="freq", axis1="label",axis2="pred",fill="label",label="label",xlabel="methods",
                                 scale_x_discrete_limits="c(\"label\",\"value\")",facet_var="methods",facet_wrap=T,width=10,height=7,nrow=1)
      
      plot.sankey_plot(grouped_raw_results1,sankey_plot_params,fig_path,str_glue("sankey_batch_effects_{unique_val}.pdf"))
    }
  }
}
  
plot.inter_diseases <- function(results,raw_results,dataset){
    
  require(data.table)
  results <- utils.get_methods(results)
  for(grouped_results in dplyr::group_split(results,train_dataset,test_dataset)){
    train_dataset <- unique(grouped_results$train_dataset)
    test_dataset <- unique(grouped_results$test_dataset)
    fig_path <- str_glue("{result_home}{experiment}/{dataset}/{train_dataset}_{test_dataset}/")
    if(!dir.exists(fig_path)){
      dir.create(fig_path,recursive=T)
    }
    print(str_glue("Generating figures for inter diseases pairs {fig_path}"))
    
    # figure_all_name <- str_glue("inter_diseases_{dataset}_all")
    # figure_assigned_name <- str_glue("cell_number_{dataset}_assigned")
    # figure_unassigned_name <- str_glue("cell_number_{dataset}_unassigned")
    

    all_results <- dplyr::select(dplyr::filter(grouped_results,is.na(assigned)),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,train_dataset,test_dataset,supervised))
    all_results$value <- as.double(all_results$value)
    all_results$label <- round(all_results$value,2)
    # plot_params <-
    # list(x="sample_num",y="value",group="methods",line_color="methods",point_color="methods",
    # facet_wrap=T,facet_var="metric",width=10,height=10)
    # plot.line_plot(all_results[all_results$metric=="ARI",],plot_params,result_home,figure_all_name)
    
    plot_params <- list(x="methods",y="value",fill="supervised",label="label",xlabel="methods",
                        facet_grid=T,dodged=T,facet_grid_x='train_dataset',facet_grid_y='test_dataset'
                        ,width=10,height=7)
    plot_params$ylabel <- "ARI"
    plot.bar_plot(all_results[all_results$metric=="ARI",],plot_params,fig_path,"all_ARI.pdf")
    plot_params$ylabel <- "AMI"
    plot.bar_plot(all_results[all_results$metric=="AMI",],plot_params,fig_path,"all_AMI.pdf")
    plot_params$ylabel <- "FMI"
    plot.bar_plot(all_results[all_results$metric=="FMI",],plot_params,fig_path,"all_FMI.pdf")
    
    
    results_assigned <- dplyr::select(dplyr::filter(grouped_results,assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,train_dataset,test_dataset,supervised))
    results_assigned$value <- as.double(results_assigned$value)
    results_assigned$label <- round(results_assigned$value,2)
    plot_params$ylabel <- "ARI"
    plot.bar_plot(results_assigned[results_assigned$metric=="ARI",],plot_params,fig_path,"assigned_ARI.pdf")
    plot_params$ylabel <- "AMI"
    plot.bar_plot(results_assigned[results_assigned$metric=="AMI",],plot_params,fig_path,"assigned_AMI.pdf")
    plot_params$ylabel <- "FMI"
    plot.bar_plot(results_assigned[results_assigned$metric=="FMI",],plot_params,fig_path,"assigned_FMI.pdf")
    # plot.line_plot(results_assigned,plot_params,result_home,figure_assigned_name)
    
    results_unassigned <- dplyr::select(dplyr::filter(grouped_results,!assigned),-c(unlabeled_pctg,assigned)) %>%
      gather("metric","value",-c(methods,train_dataset,test_dataset,supervised))
    results_unassigned$value <- as.double(results_unassigned$value)
    results_unassigned$label <- round(results_unassigned$value,2)
    plot_params$ylabel <- "ARI"
    plot.bar_plot(results_unassigned[results_unassigned$metric=="ARI",],plot_params,fig_path,"unassigned_ARI.pdf")
    plot_params$ylabel <- "AMI"
    plot.bar_plot(results_unassigned[results_unassigned$metric=="AMI",],plot_params,fig_path,"unassigned_AMI.pdf")
    plot_params$ylabel <- "FMI"
    plot.bar_plot(results_unassigned[results_unassigned$metric=="FMI",],plot_params,fig_path,"unassigned_FMI.pdf")
    
    # plot.line_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
    
    ######
    # figure_name <- str_glue("inter_diseases_{dataset}_unlabeled_pctg.pdf")
    results_unlabeled_pctg <- dplyr::filter(dplyr::select(grouped_results[grouped_results$supervised,],methods,train_dataset,test_dataset,unlabeled_pctg),!is.na(unlabeled_pctg))%>%
      gather("metric","value",-c(methods,train_dataset,test_dataset))
    results_unlabeled_pctg$label <- round(results_unlabeled_pctg$value,2)
    plot_params <- list(x="methods",y="value",label="label",xlabel="methods",ylabel="unlabeled_pctg",
                        facet_grid=T,dodged=T,facet_grid_x='train_dataset',facet_grid_y='test_dataset'
                        ,width=10,height=7)
    
    plot.bar_plot(results_unlabeled_pctg,plot_params,fig_path,"unlabeled_pctg.pdf")
  }
  
  for(grouped_raw_results in dplyr::group_split(raw_results,train_dataset,test_dataset)){
    ###sankey plot
    train_dataset <- unique(grouped_raw_results$train_dataset)
    test_dataset <- unique(grouped_raw_results$test_dataset)
    fig_path <- str_glue("{result_home}{experiment}/{dataset}/{train_dataset}_{test_dataset}/")
    
    methods <- colnames(dplyr::select(grouped_raw_results,-c("label","train_dataset","test_dataset")))
    
    grouped_raw_results1 <- gather(grouped_raw_results,"methods","pred",-c(label,"train_dataset","test_dataset")) %>%
      group_by(methods,label,pred) %>%
      summarize(freq=n()) %>%
      filter(freq>0)
    
    sankey_plot_params <- list(y="freq", axis1="label",axis2="pred",fill="label",label="label",xlabel="methods",
                               scale_x_discrete_limits="c(\"label\",\"value\")",facet_var="methods",facet_wrap=T,width=10,height=7,nrow=1)
    
    plot.sankey_plot(grouped_raw_results1,sankey_plot_params,fig_path,str_glue("sankey.pdf"))
  }
  
}
  
  


plot.celltype_complexity <- function(results,raw_results,dataset){
  
}

plot.inter_species <- function(results,raw_results,dataset){
  
}

plot.random_noise <- function(results,raw_results,dataset){
  
}

plot.inter_protocol <- function(results,raw_results,dataset){
  
}


plot.bar_plot <- function(results,params,fig_path,fig_name){
  x <- params$x
  y <- params$y
  label <- params$label
  plot_str <- "ggplot(results,aes_string(x= x,y = y,fill = {params$fill})) +
                geom_bar("
  if(params$dodged){
    plot_str <- str_glue(
    "{plot_str} position=\"dodge\",")
  } 
  plot_str <- str_glue("
    {plot_str} width = 0.8, stat = \"identity\") +
    scale_fill_brewer(palette = \"Pastel1\") +
    geom_text(aes_string(x = x,label = label),
              vjust = 0.1,size=2,hjust=0.7,    
              colour = \"brown\"")
  if(params$dodged){
    plot_str <- str_glue("{plot_str},
              position = position_dodge(width = 0.5),
              show.legend = F")
    }
  plot_str <- str_glue("{plot_str})+
                        labs(x=\"{params$xlabel}\",y = \"{params$ylabel}\") +
                        scale_x_discrete(guide = guide_axis(n.dodge=3))")
  
  if(!is_null(params$facet_wrap)&&params$facet_wrap){
    plot_str <- str_glue("{plot_str} + 
                          facet_wrap(~ {params$facet_var},nrow={params$nrow},labeller = label_both)")
  }
  
  if(!is_null(params$facet_grid)&&params$facet_grid){
    plot_str <- str_glue("{plot_str} + 
                          facet_grid({params$facet_grid_x} ~ {params$facet_grid_y},labeller = label_both)")
  }
  eval(parse(text = plot_str))
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'pdf',
    path = fig_path,
    width = params$width,
    height = params$height,
    units = "in"
  )
  
}


plot.sankey_plot <- function(raw_results,plot_params,fig_path,fig_name){
  require(ggalluvial)
  library(RColorBrewer)
  colors <- colorRampPalette(brewer.pal(8, "Set2"))(40)
  n_methods <- length(unique(raw_results$methods))
  plot_str <- str_glue("ggplot(as.data.frame(raw_results),aes(y={plot_params$y},axis1={plot_params$axis1},axis2 = {plot_params$axis2})) +
                       geom_alluvium(aes(fill = {plot_params$fill}), width = 1/12) +
                       geom_stratum(width = 1/12, fill = \"black\", color = \"grey\") +
                       geom_label(stat = \"stratum\", aes(label = after_stat(stratum))) +
                       scale_x_discrete(limits = {plot_params$scale_x_discrete_limits}, expand = c(.05, .05)) +
                       scale_fill_manual(values = colors)")
  if(!is_null(plot_params$facet_wrap)&&plot_params$facet_wrap){
    plot_str <- str_glue("{plot_str} + 
                          facet_wrap(~ {plot_params$facet_var},nrow={n_methods%/%2},labeller = label_both)")
  }
  
  if(!is_null(plot_params$facet_grid)&&plot_params$facet_grid){
    plot_str <- str_glue("{plot_str} + 
                          facet_grid({plot_params$facet_grid_x} ~ {plot_params$facet_grid_y},labeller = label_both)")
  }
  
  plot_str <- str_glue("{plot_str} + 
                        theme(strip.text = element_text(size={3*n_methods}),legend.title = element_text(color = \"blue\", size = 3*{n_methods}),
                        legend.text = element_text(size={3*n_methods}))")
  eval(parse(text = plot_str))
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'pdf',
    path = fig_path,
    width = 3.5*n_methods,
    height = 3.5*n_methods,
    units = "in"
  )
}

plot.line_plot <- function(results,params,fig_path,fig_name){
  x <- params$x
  y <- params$y
  group <- params$group
  line_color <- params$line_color
  point_color <- params$point_color
  width <- params$width
  height <- params$height
  plot_str <- 
  "ggplot(results, aes_string(x=x, y=y, group=group)) +
    geom_line(aes_string(color=line_color))+
    geom_point(aes_string(color=point_color))"
  if(params$facet_wrap){   
    plot_str <- str_glue("{plot_str}+
      facet_wrap(~{params$facet_var})"
    )
  }
  print(plot_str)
  eval(parse(text = plot_str))
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'pdf',
    path = fig_path,
    width = width,
    height = height,
    units = "in"
  )
}

plot.heatmap_plot <- function(results,params,fig_path,fig_name){
  require(hrbrthemes)
  x <- params$x
  y <- params$y
  fill <- params$fill
  
  plot_str <- "
  ggplot(results, aes_string(x,y, fill=fill)) + 
    geom_tile() +
    scale_fill_distiller(palette = \"RdPu\") +
    theme_ipsum()"
  if(params$facet_wrap){
    plot_str <- str_glue("{plot_str}+
      facet_wrap(~{params$facet_var},nrow={length(unique(results$metric))})")
  }
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'pdf',
    path = fig_path,
    width = params$width,
    height = params$height,
    units = "in"
  )
}

