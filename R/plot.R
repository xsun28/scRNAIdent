source("R/utils.R")

plot.plot <- function(experiment,results,raw_results){
  dataset <- experiments.data[[experiment]]
  if(purrr::is_null(results)){
    results <- read_rds(str_glue("{result_home}/{experiment}_{dataset}_results.rds"))
  }
  if(purrr::is_null(raw_results)){
    raw_results <- read_rds(str_glue("{result_home}/{experiment}_{dataset}_raw_results.rds"))
  }
  switch(experiment,
    simple_accuracy = plot.simple_accuracy(results,raw_results),
    cell_number = plot.cell_number(results,raw_results),
    sequencing_depth = plot.sequencing_depth(results,raw_results),
    celltype_structure = plot.celltype_structure(results,raw_results),
    batch_effects = plot.batch_effects(results,raw_results)
  )
}

plot.simple_accuracy <- function(results,raw_results){
  require(data.table)
  dataset <- experiments.data$simple_accuracy
  figure_all_name <- str_glue("simple_accuracy_{dataset}_all_metrics.pdf")
  figure_assigned_name <- str_glue("simple_accuracy_{dataset}_assigned_metrics.pdf")
  figure_unassigned_name <- str_glue("simple_accuracy_{dataset}_unassigned_metrics.pdf")
  
  all_results <- bind_rows(results[grepl("^all_.*",names(results))])
  all_results[,'methods'] <- rownames(all_results)
  results_assigned <- bind_rows(results[grepl("^assigned_.*",names(results))])
  results_assigned[,'methods'] <- rownames(results_assigned)
  results_unassigned <- bind_rows(results[grepl("^unassigned_.*",names(results))])
  results_unassigned[,'methods'] <- rownames(results_unassigned)

  all_results_table <- data.table(melt(all_results,
                                          id.vars = c('methods'),    
                                          measure.vars=c('ARI','AMI','FMI')))
  all_results_table$methods <- factor(all_results_table$methods)
  all_results_table$label <- round(all_results_table$value,2)
  plot_params <- list(x="methods",y="value",fill="variable",label="label",
                      facet_var="variable",dodged=F,facet_wrap=T,width=10,height=7,nrow=3)
  plot.bar_plot(all_results_table,plot_params,result_home,figure_all_name)
  
  results_assign_table <- data.table(melt(results_assigned,
                                   id.vars = c('methods'),    
                                   measure.vars=c('ARI','AMI','FMI')))
  results_assign_table$methods <- factor(results_assign_table$methods)
  results_assign_table$label <- round(results_assign_table$value,2)
  plot.bar_plot(results_assign_table,plot_params,result_home,figure_assigned_name)
  
  results_unassign_table <- data.table(melt(results_unassigned,
                                          id.vars = c('methods'),    
                                          measure.vars=c('ARI','AMI','FMI')))
  results_unassign_table$methods <- factor(results_unassign_table$methods)
  results_unassign_table$label <- round(results_unassign_table$value,2)
  plot.bar_plot(results_unassign_table,plot_params,result_home,figure_unassigned_name)
  
  #####unlabeled pctg 
  figure_name <- str_glue("simple_accuracy_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- results[["all_assign_results"]]['unlabeled_pctg']
  results_unlabeled_pctg[,'methods'] <- factor(rownames(results_unlabeled_pctg))
  results_unlabeled_pctg[,'label'] <- round(results_unlabeled_pctg$unlabeled_pctg,2)
  plot_params <- list(x="methods",y="unlabeled_pctg",label="label",fill="methods",dodged=F,facet_wrap=F,width=6,height=4)
  plot.bar_plot(results_unlabeled_pctg,plot_params,result_home,figure_name)
  
 #####sankey plot
  figure_name <- str_glue("simple_accuracy_{dataset}_sankey.pdf")
  methods <- colnames(dplyr::select(raw_results,-label))
  raw_results1 <- gather(raw_results,"methods","pred",-label) %>%
                    group_by(methods,label,pred) %>%
                      summarize(freq=n()) %>%
                        filter(freq>0)
  plot.sankey_plot(raw_results1,"label","pred",result_home,figure_name)                   
}

plot.cell_number <- function(results,raw_results){
  dataset <- experiments.data$cell_number
  figure_all_name <- str_glue("cell_number_{dataset}_all_metrics.pdf")
  figure_assigned_name <- str_glue("cell_number_{dataset}_assigned_metrics.pdf")
  figure_unassigned_name <- str_glue("cell_number_{dataset}_unassigned_metrics.pdf")
  
  
  results <- purrr::map(results,utils.get_methods)
  all_results <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,is.na(assigned)))),-c(unlabeled_pctg,assigned)) %>%
                gather("metric","value",-c(methods,sample_num))
  
  plot_params <- list(x="sample_num",y="value",group="methods",line_color="methods",point_color="methods",
                      facet_wrap=T,facet_var="metric",width=10,height=10)
  plot.line_plot(all_results,plot_params,result_home,figure_all_name)
  
  results_assigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,sample_num))
  plot.line_plot(results_assigned,plot_params,result_home,figure_assigned_name)
  
  results_unassigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,!assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,sample_num))
  
  plot.line_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
  
  ######
  figure_name <- str_glue("cell_number_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- dplyr::filter(dplyr::select(results$assign_results,methods,sample_num,unlabeled_pctg),!is.na(unlabeled_pctg))
  
  plot_params <- list(x="sample_num",y="unlabeled_pctg",group="methods",line_color="methods",point_color="methods",
                      facet_wrap=F,width=2*length(unique(results_unlabeled_pctg$methods)),
                      height=2*length(unique(results_unlabeled_pctg$methods)))
  plot.line_plot(results_unlabeled_pctg,plot_params,result_home,figure_name)
}


plot.celltype_structure <- function(results,raw_results){
  dataset <- experiments.data$celltype_structure
  figure_all_name <- str_glue("celtype_structure_{dataset}_all_metrics.pdf")
  figure_assigned_name <- str_glue("celtype_structure_{dataset}_assigned_metrics.pdf")
  figure_unassigned_name <- str_glue("celtype_structure_{dataset}_unassigned_metrics.pdf")
  
  results <- purrr::map(results,utils.get_methods)
  all_results <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,is.na(assigned)))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,level))
  all_results$level <- factor(all_results$level)
  
  plot_params <- list(x="level",y="value",group="methods",line_color="methods",point_color="methods",
                      facet_wrap=T,facet_var="metric",width=10,height=10)
  plot.line_plot(all_results,plot_params,result_home,figure_all_name)
  
  results_assigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,level))
  results_assigned$level <- factor(results_assigned$level)
  plot.line_plot(results_assigned,plot_params,result_home,figure_assigned_name)
  
  results_unassigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,!assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,level))
  results_unassigned$level <- factor(results_unassigned$level)
  plot.line_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
  
  ######
  figure_name <- str_glue("celltype_structure_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- dplyr::filter(dplyr::select(results$assign_results,methods,level,unlabeled_pctg),!is.na(unlabeled_pctg))
  results_unlabeled_pctg$level <- factor(results_unlabeled_pctg$level)
  plot_params <- list(x="level",y="unlabeled_pctg",group="methods",line_color="methods",point_color="methods",
                      facet_wrap=F,width=2*length(unique(results_unlabeled_pctg$methods)),
                      height=2*length(unique(results_unlabeled_pctg$methods)))
  plot.line_plot(results_unlabeled_pctg,plot_params,result_home,figure_name)
}


plot.sequencing_depth <- function(results,raw_results){
  dataset <- experiments.data$sequencing_depth
  
  results <- purrr::map(results,utils.get_methods)%>%
    purrr::map(~{.$quantile <- factor(.$quantile)
    return(.)
    })
  #######
  figure_all_name <- str_glue("sequencing_depth_{dataset}_all_metrics.pdf")
  figure_assigned_name <- str_glue("sequencing_depth_{dataset}_assigned_metrics.pdf")
  figure_unassigned_name <- str_glue("sequencing_depth_{dataset}_unassigned_metrics.pdf")
  
  
  all_results <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,is.na(assigned)))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,quantile))
  all_results$label <- round(all_results$value,2)
  plot_params <- list(x="methods",y="value",fill="quantile",label="label",
                      facet_wrap=T,dodged=T,facet_var="metric",width=10,height=7,nrow=2)
  plot.bar_plot(all_results,plot_params,result_home,figure_all_name)
  # plot.heatmap_plot(all_results,"quantile","methods",result_home,figure_all_name)
  
  results_assigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,quantile))
  results_assigned$label <- round(results_assigned$value,2)
  plot.bar_plot(results_assigned,plot_params,result_home,figure_assigned_name)
  
  results_unassigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,!assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,quantile))
  results_unassigned$label <- round(results_unassigned$value,2)
  plot.bar_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
  
  #######
  figure_name <- str_glue("sequencing_depth_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- dplyr::filter(dplyr::select(results$assign_results,methods,quantile,unlabeled_pctg),!is.na(unlabeled_pctg))
  results_unlabeled_pctg[,'label'] <- round(results_unlabeled_pctg$unlabeled_pctg,2)
  plot_params <- list(x="methods",y="unlabeled_pctg",label="label",fill="quantile",dodged=T,facet_wrap=F,width=6,height=4)
  plot.bar_plot(results_unlabeled_pctg,plot_params,result_home,figure_name)
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
  #####
  # figure_name <- str_glue("sequencing_depth_{dataset}_sankey.pdf")
  # n_methods <- length(colnames(dplyr::select(raw_results,-c(label,quantile))))
  # raw_results1 <- gather(raw_results,"methods","pred",-c(label,quantile)) %>%
  #   group_by(methods,label,pred,quantile) %>%
  #   summarize(freq=n()) %>%
  #   filter(freq>0)
  # raw_results1$quantile <- factor(raw_results1$quantile)
  # ggplot(as.data.frame(raw_results1),
  #        aes(y="freq",axis1 = label, axis2 = pred)) +
  #   geom_alluvium(aes(fill = label), width = 1/12) +
  #   geom_stratum(width = 1/12, fill = "black", color = "grey") +
  #   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  #   scale_x_discrete(limits = c("label", "value"), expand = c(.05, .05)) +
  #   scale_fill_brewer(type = "qual", palette = "Set1") +
  #   facet_grid(quantile ~ methods,labeller = label_both)+
  #   theme(strip.text = element_text(size=n_methods*3),legend.title = element_text(color = "blue", size = n_methods*3),
  #         legend.text = element_text(size=n_methods*3))
  # ggsave(
  #   figure_name,
  #   plot = last_plot(),
  #   device = 'pdf',
  #   path = result_home,
  #   width = n_methods*8,
  #   height = n_methods*4,
  #   units = "in"
  # )
}

plot.batch_effects <- function(results,raw_results){
  dataset <- experiments.data$batch_effects
  
  results <- purrr::map(results,utils.get_methods)
  
  figure_all_name <- str_glue("batch_effects_{dataset}_all_metrics.pdf")
  figure_assigned_name <- str_glue("batch_effects_{dataset}_assigned_metrics.pdf")
  figure_unassigned_name <- str_glue("batch_effects_{dataset}_unassigned_metrics.pdf")

  all_results <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,is.na(assigned)))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,batch_effects_removed))
  
  all_results$label <- round(all_results$value,2)
  plot_params <- list(x="methods",y="value",fill="batch_effects_removed",label="label",
                      facet_wrap=T,dodged=T,facet_var="metric",width=10,height=7,nrow=3)
  plot.bar_plot(all_results,plot_params,result_home,figure_all_name)
  
  # plot.heatmap_plot(all_results,"batch_effects_removed","methods",result_home,figure_all_name)
  
  results_assigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,batch_effects_removed))
  results_assigned$label <- round(results_assigned$value,2)
  plot.bar_plot(results_assigned,plot_params,result_home,figure_assigned_name)
  

  results_unassigned <- dplyr::select(bind_rows(purrr::map(results,~dplyr::filter(.,!assigned))),-c(unlabeled_pctg,assigned)) %>%
    gather("metric","value",-c(methods,batch_effects_removed))
  results_unassigned$label <- round(results_unassigned$value,2)
  plot.bar_plot(results_unassigned,plot_params,result_home,figure_unassigned_name)
  

  figure_name <- str_glue("batch_effects_{dataset}_unlabeled_pctg.pdf")
  results_unlabeled_pctg <- dplyr::filter(dplyr::select(results$assign_results,methods,batch_effects_removed,unlabeled_pctg),!is.na(unlabeled_pctg))
  results_unlabeled_pctg$label <- round(results_unlabeled_pctg$unlabeled_pctg,2)
  plot_params <- list(x="methods",y="unlabeled_pctg",fill="batch_effects_removed",label="label",
                      facet_wrap=F,dodged=T,width=10,height=7)
  plot.bar_plot(results_unlabeled_pctg,plot_params,result_home,figure_name)
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
                        scale_x_discrete(guide = guide_axis(n.dodge=3))")
  
  if(params$facet_wrap){
    plot_str <- str_glue("{plot_str} + 
                          facet_wrap(~ {params$facet_var},nrow={params$nrow})")
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


plot.sankey_plot <- function(raw_results,label,pred,fig_path,fig_name){
  require(ggalluvial)
  n_methods <- length(unique(raw_results$methods))
  ggplot(as.data.frame(raw_results),
         aes_string(y="freq",axis1 = label, axis2 = pred)) +
    geom_alluvium(aes_string(fill = label), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("label", "value"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    facet_wrap(~ methods,nrow=n_methods%/%2) +
    theme(strip.text = element_text(size=3*n_methods),legend.title = element_text(color = "blue", size = 3*n_methods),
        legend.text = element_text(size=3*n_methods))
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