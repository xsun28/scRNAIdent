source("R/utils.R")

plot.plot <- function(experiment,results,raw_results){
  if(purrr::is_null(results)){
    results <- read_rds(str_glue("{result_home}/{experiment}_results.rds"))
  }
  if(purrr::is_null(raw_results)){
    raw_results <- read_rds(str_glue("{result_home}/{experiment}_raw_results.rds"))
  }
  switch(experiment,
    simple_accuracy = plot.simple_accuracy(results,raw_results),
    cell_number = plot.cell_number(results,raw_results),
    sequencing_depth = plot.sequencing_depth(results,raw_results),
    cell_type = plot.cell_type(results,raw_results),
    batch_effects = plot.batch_effects(results,raw_results)
  )
}

plot.simple_accuracy <- function(results,raw_results){
  require(data.table)
  figure_name <- "simple_accuracy_metrics.png"
  results_assigned <- bind_rows(results[grepl("^assigned_.*",names(results))])
  results_assigned[,'methods'] <- rownames(results_assigned)
  results_unassigned <- bind_rows(results[grepl("^unassigned_.*",names(results))])
  results_unassigned[,'methods'] <- rownames(results_unassigned)
  results1 <- bind_rows(results_assigned,results_unassigned)
  results_table <- data.table(melt(results1,
                id.vars = c('methods','assigned'),    
                measure.vars=c('ARI','AMI','FMI')))
  results_table$methods <- factor(results_table$methods)
  results_table$label <- round(results_table$value,2)
  plot.dodge_bar_plot(results_table,"methods","value",result_home,figure_name)
  
  #####unlabeled pctg 
  figure_name <- "simple_accuracy_unlabeled_pctg.png"
  results_unlabeled_pctg <- results[["all_assign_results"]]['unlabeled_pctg']
  results_unlabeled_pctg[,'methods'] <- factor(rownames(results_unlabeled_pctg))
  results_unlabeled_pctg[,'label'] <- round(results_unlabeled_pctg$unlabeled_pctg,2)
  plot.bar_plot(results_unlabeled_pctg,'methods','unlabeled_pctg',result_home,figure_name)
  
 #####sankey plot
  figure_name <- "simple_accuracy_sankey.png"
  methods <- colnames(select(raw_results,-label))
  raw_results1 <- gather(raw_results,"methods","pred",-label) %>%
                    group_by(methods,label,pred) %>%
                      summarize(freq=n()) %>%
                        filter(freq>0)
  plot.sankey_plot(raw_results1,"label","pred",result_home,figure_name)                   
}

plot.cell_number <- function(results,raw_results){
  figure_name <- "cell_number_metrics.png"

  results <- purrr::map(results,utils.get_methods)
  results1 <- select(bind_rows(purrr::map(results,~dplyr::filter(.,!is.na(assigned)))),-unlabeled_pctg) %>%
                gather("metric","value",-c(methods,sample_num,assigned))

  plot.line_plot(results1,"sample_num","value",result_home,figure_name)
  
  ######
  figure_name <- "cell_number_unlabeled_pctg.png"
  results_unlabeled_pctg <- dplyr::filter(select(results$assign_results,methods,sample_num,unlabeled_pctg),!is.na(unlabeled_pctg))
  ggplot(results_unlabeled_pctg, aes(x=sample_num, y=unlabeled_pctg, group=methods)) +
    geom_line(aes(color=methods))+
    geom_point(aes(color=methods))
  ggsave(
    figure_name,
    plot = last_plot(),
    device = 'png',
    path = result_home,
    width = 2*length(unique(results_unlabeled_pctg$methods)),
    height = 2*length(unique(results_unlabeled_pctg$methods)),
    units = "in"
  )
}

plot.sequencing_depth <- function(results,raw_results){
  results <- purrr::map(results,utils.get_methods)%>%
    purrr::map(~{.$quantile <- factor(.$quantile)
    return(.)
    })
  #######
  figure_name <- "sequencing_depth_metrics.png"
  results1 <- select(bind_rows(purrr::map(results,~dplyr::filter(.,!is.na(assigned)))),-unlabeled_pctg) %>%
    gather("metric","value",-c(methods,quantile,assigned))
  plot.heatmap_plot(results1,"quantile","methods",result_home,figure_name)
  #######
  figure_name <- "sequencing_depth_unlabeled_pctg.png"
  results_unlabeled_pctg <- dplyr::filter(select(results$assign_results,methods,quantile,unlabeled_pctg),!is.na(unlabeled_pctg))
  ggplot(results_unlabeled_pctg, aes(quantile,methods, fill=unlabeled_pctg)) + 
    geom_tile() +
    scale_fill_distiller(palette = "RdPu") +
    theme_ipsum()
  ggsave(
    figure_name,
    plot = last_plot(),
    device = 'png',
    path = result_home,
    width = 2*length(unique(results_unlabeled_pctg$methods)),
    height = 2*length(unique(results_unlabeled_pctg$methods)),
    units = "in"
  )
  #####
  figure_name <- "sequencing_depth_sankey.png"
  n_methods <- length(colnames(select(raw_results,-c(label,quantile))))
  raw_results1 <- gather(raw_results,"methods","pred",-c(label,quantile)) %>%
    group_by(methods,label,pred,quantile) %>%
    summarize(freq=n()) %>%
    filter(freq>0)
  raw_results1$quantile <- factor(raw_results1$quantile)
  ggplot(as.data.frame(raw_results1),
         aes(y="freq",axis1 = label, axis2 = pred)) +
    geom_alluvium(aes(fill = label), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("label", "value"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    facet_grid(quantile ~ methods,labeller = label_both)+
    theme(strip.text = element_text(size=n_methods*3),legend.title = element_text(color = "blue", size = n_methods*3),
          legend.text = element_text(size=n_methods*3))
  ggsave(
    figure_name,
    plot = last_plot(),
    device = 'png',
    path = result_home,
    width = n_methods*8,
    height = n_methods*4,
    units = "in"
  )
}

plot.batch_effects <- function(results,raw_results){
  results <- purrr::map(results,utils.get_methods)
  figure_name <- "batch_effects_metrics.png"
  figure_name <- "batch_effects_unlabeled_pctg.png"
  figure_name <- "batch_effects_sankey.png"
}

plot.dodge_bar_plot <- function(results,x,y,fig_path,fig_name){
  ggplot(results,aes_string(x=x,y = y,fill = "variable")) + 
    geom_bar(position = "dodge",width = 0.8, stat = "identity") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_text(aes_string(x = x,label = "label"),
              vjust = 0.1,size=2,hjust=0.7,    
              colour = "brown",
              position = position_dodge(width = 0.5),
              show.legend = F)+
    scale_x_discrete(guide = guide_axis(n.dodge=3))+
    facet_wrap(~ assigned,nrow=1,labeller = label_both)
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'png',
    path = fig_path,
    width = 10,
    height = 7,
    units = "in"
  )
}

plot.bar_plot <- function(results,x,y,fig_path,fig_name){
  ggplot(results,aes_string(x = x,y = y)) +
    geom_bar(width = 0.8, stat = "identity") +
    scale_fill_brewer(palette = "Blues") +
    geom_text(aes_string(x = x,label = "label"),
              vjust = 0.1,size=6,hjust=0.7,    
              colour = "brown")+
    scale_x_discrete(guide = guide_axis(n.dodge=3))
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'png',
    path = fig_path,
    width = 10,
    height = 7,
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
    device = 'png',
    path = fig_path,
    width = 5*n_methods,
    height = 5*n_methods,
    units = "in"
  )
}

plot.line_plot <- function(results,x,y,fig_path,fig_name){
  ggplot(results, aes_string(x=x, y=y, group="methods")) +
    geom_line(aes(color=methods))+
    geom_point(aes(color=methods))+
    facet_grid(metric~assigned,labeller = label_both)
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'png',
    path = fig_path,
    width = 10,
    height = 10,
    units = "in"
  )
}

plot.heatmap_plot <- function(results,x,y,fig_path,fig_name){
  require(hrbrthemes)
  ggplot(results, aes_string(x,y, fill= "value")) + 
    geom_tile() +
    scale_fill_distiller(palette = "RdPu") +
    theme_ipsum()+
    facet_wrap(~metric,nrow=length(unique(results$metric)))
  ggsave(
    fig_name,
    plot = last_plot(),
    device = 'png',
    path = fig_path,
    width = 8,
    height = 10,
    units = "in"
  )
}