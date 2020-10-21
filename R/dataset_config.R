raw_datasets <- list(PBMC="PBMC_AllCells_withLabels.Rdata",
                     pancreas=c("Muraro_pancreas.Rdata","Segerstolpe_pancreas.Rdata","Xin_pancreas.Rdata"))

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