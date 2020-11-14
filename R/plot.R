plot.plot <- function(experiment,results){
  switch(experiment,
    simple_accuracy = plot.simple_accuracy(results),
    cell_number = plot.cell_number(results),
    sequencing_depth = plot.sequencing_depth(results),
    cell_type = plot.cell_type(results),
    batch_effects = plot.batch_effects(results)
  )
}

plot.simple_accuracy <- function(results){
  
  figure_path <- str_glue("{results_home}/{simple_accuracy}.png")
  
  
}