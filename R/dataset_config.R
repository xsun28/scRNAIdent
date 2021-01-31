raw_datasets <- list(PBMC=c("PBMC_AllCells_withLabels.RData","GSE96583_batch1_3_samples.RData","GSE96583_8_Stim_Pats.RData","GSE96583_8_Ctrl_Pats.RData"),
                     pancreas=c("Muraro_pancreas.RData","Segerstolpe_pancreas.RData","Xin_pancreas.RData"),
                     adasd="AD_autism_data.RData",
                     midbrain="midbrain.RData",
                     lung_cancer="lung_cancer.RData")

datasets <- list(PBMC=list("PBMC_AllCells_withLabels.RDS","GSE96583_batch1_3_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS"),
                 pancreas=list(muraro="Muraro_pancreas_clean.RDS",
                               seger="Segerstolpe_pancreas_clean.RDS",
                               xin="Xin_pancreas_clean.RDS"),
                 ADASD=list(AD="ADASD_AD.RDS",autism="ADASD_autism.RDS"),
                 midbrain=list(human="midbrain_human.RDS",mouse="midbrain_mouse.RDS"))

batch_effects_free_datasets <-list(pancreas=list(muraro="Muraro_pancreas_batch_effects_free.RDS",
                                                 seger="Muraro_pancreas_batch_effects_free.RDS",
                                                 xin="Xin_pancreas_batch_effects_free.RDS"
))

dataset.properties <- list(pancreas=list(sample_threshold=100,cell_types=c('beta','alpha','delta',
                                                                           'acinar','ductal','gamma',
                                                                           'endothelial')),
                           PBMC=list(sample_threshold=100),
                           ADASD=list(sample_threshold=100),
                           midbrain=list(sample_threshold=100),
                           lung_cancer=list(sample_threshold=100)
                           )