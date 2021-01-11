raw_datasets <- list(PBMC="PBMC_AllCells_withLabels.RData",
                     pancreas=c("Muraro_pancreas.RData","Segerstolpe_pancreas.RData","Xin_pancreas.RData"))

datasets <- list(PBMC="PBMC_AllCells_withLabels.RDS",
                 pancreas=list(muraro="Muraro_pancreas_clean.RDS",
                               seger="Segerstolpe_pancreas_clean.RDS",
                               xin="Xin_pancreas_clean.RDS"))
batch_effects_free_datasets <-list(pancreas=list(muraro="Muraro_pancreas_batch_effects_free.RDS",
                                                 seger="Muraro_pancreas_batch_effects_free.RDS",
                                                 xin="Xin_pancreas_batch_effects_free.RDS"
))

dataset.properties <- list(pancreas=list(sample_threshold=100,cell_types=c('beta','alpha','delta',
                                                                           'acinar','ductal','gamma',
                                                                           'endothelial')),
                           PBMC=list(sample_threshold=100)
                           )