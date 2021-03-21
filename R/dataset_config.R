raw_datasets <- list(PBMC=c("PBMC_AllCells_withLabels.RData","GSE96583_batch1_3_samples.RData","GSE96583_8_Stim_Pats.RData","GSE96583_8_Ctrl_Pats.RData"),
                     # pancreas=c("Muraro_pancreas.RData","Segerstolpe_pancreas.RData","Xin_pancreas.RData"),
                     pancreas=c("Muraro_pancreas.RData","Segerstolpe_pancreas.RData"),
                     ADASD="AD_autism_data.RData",
                     midbrain="midbrain.RData",
                     cellbench="cellbench.RData")

datasets <- list(PBMC=list("PBMC_AllCells_withLabels.RDS","GSE96583_batch1_3_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS"),
                 pancreas=list(muraro="Muraro_pancreas_clean.RDS",
                               seger="Segerstolpe_pancreas_clean.RDS"
                               #xin="Xin_pancreas_clean.RDS"
                               ),
                 ADASD=list(AD="ADASD_AD.RDS",autism="ADASD_autism.RDS"),
                 midbrain=list(human="midbrain_human.RDS",mouse="midbrain_mouse.RDS"),
                 cellbench=list(tenx="cellbench_10x.RDS",CELseq2="cellbench_CELseq2.RDS",Dropseq="cellbench_Dropseq.RDS")
                 )

dataset.properties <- list(
                          PBMC_AllCells_withLabels=list(sample_threshold=100,gene_name_type="SYMBOL"),
                          GSE96583_batch1_3_samples=list(sample_threshold=100,gene_name_type="SYMBOL",cell_type_map="PBMC_celltype_map.csv"),
                          GSE96583_8_Stim_Pats=list(sample_threshold=100,gene_name_type="SYMBOL",cell_type_map="PBMC_celltype_map.csv"),
                          GSE96583_8_Ctrl_Pats=list(sample_threshold=100,gene_name_type="SYMBOL",cell_type_map="PBMC_celltype_map.csv"),
                          ADASD_AD=list(sample_threshold=100,gene_name_type="SYMBOL"),
                          ADASD_autism=list(sample_threshold=100,gene_name_type="SYMBOL"),
                          Muraro_pancreas=list(sample_threshold=20,gene_name_type="SYMBOL",
                                                              cell_types=c('beta','alpha','delta',
                                                                          'acinar','ductal','gamma',
                                                                          'endothelial')),
                          Segerstolpe_pancreas=list(sample_threshold=20,gene_name_type="SYMBOL",
                                               cell_types=c('beta','alpha','delta',
                                                            'acinar','ductal','gamma',
                                                            'endothelial')),
                          Xin_pancreas=list(sample_threshold=20,gene_name_type="SYMBOL",
                                               cell_types=c('beta','alpha','delta',
                                                            'acinar','ductal','gamma',
                                                            'endothelial')),
                          midbrain_human=list(sample_threshold=20,gene_name_type="SYMBOL"),
                          midbrain_mouse=list(sample_threshold=20,gene_name_type="SYMBOL"),
                          cellbench_10x=list(sample_threshold=20,gene_name_type="ENSEMBL"),
                          cellbench_CELseq2=list(sample_threshold=20,gene_name_type="ENSEMBL"),
                          cellbench_Dropseq=list(sample_threshold=20,gene_name_type="ENSEMBL")
                      )

