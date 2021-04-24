raw_datasets <- list(PBMC=c("PBMC_AllCells_withLabels.RData","GSE96583_batch1_3_samples.RData","GSE96583_8_Stim_Pats.RData","GSE96583_8_Ctrl_Pats.RData"),
                     # pancreas=c("Muraro_pancreas.RData","Segerstolpe_pancreas.RData","Xin_pancreas.RData"),
                     pancreas=c("Muraro_pancreas.RData","Segerstolpe_pancreas.RData"),
                     ADASD="AD_autism_data.RData",
                     midbrain="midbrain.RData",
                     cellbench="cellbench.RData")

dataset.datasets <- list("PBMC_AllCells_withLabels.RDS","GSE96583_batch1_3_samples.RDS","GSE96583_8_Stim_Pats.RDS","GSE96583_8_Ctrl_Pats.RDS",
                 "Muraro_pancreas_clean.RDS","Segerstolpe_pancreas_clean.RDS","Xin_pancreas_clean.RDS",
                 "ADASD_AD.RDS","ADASD_autism.RDS",
                 "midbrain_human.RDS","midbrain_mouse.RDS"
                 # cellbench=list(tenx="cellbench_10x.RDS",CELseq2="cellbench_CELseq2.RDS",Dropseq="cellbench_Dropseq.RDS")
                 )

dataset.interdatasets <- list(
                              PBMC1 = list(train_dataset="PBMC_AllCells_withLabels.RDS", test_dataset="GSE96583_batch1_3_samples.RDS"),
                              PBMC2 = list(train_dataset="PBMC_AllCells_withLabels.RDS",  test_dataset="GSE96583_8_Ctrl_Pats.RDS"),
                              PBMC3 = list(train_dataset="PBMC_AllCells_withLabels.RDS",  test_dataset="GSE96583_8_Stim_Pats.RDS"),
                              PBMC4 = list(train_dataset="GSE96583_8_Ctrl_Pats.RDS",  test_dataset="PBMC_AllCells_withLabels.RDS"),
                              PBMC5 = list(train_dataset="GSE96583_8_Ctrl_Pats.RDS",  test_dataset="GSE96583_8_Stim_Pats.RDS"),
                              PBMC6 = list(train_dataset="GSE96583_8_Ctrl_Pats.RDS",  test_dataset="GSE96583_batch1_3_samples.RDS"),
                              PBMC7 = list(train_dataset="GSE96583_batch1_3_samples.RDS",  test_dataset="PBMC_AllCells_withLabels.RDS"),
                              PBMC8 = list(train_dataset="GSE96583_batch1_3_samples.RDS",  test_dataset="GSE96583_8_Stim_Pats.RDS"),
                              PBMC9 = list(train_dataset="GSE96583_batch1_3_samples.RDS",  test_dataset="GSE96583_8_Ctrl_Pats.RDS"),
                              PBMC10 = list(train_dataset="GSE96583_8_Stim_Pats.RDS",  test_dataset="PBMC_AllCells_withLabels.RDS"),
                              PBMC11 = list(train_dataset="GSE96583_8_Stim_Pats.RDS",  test_dataset="GSE96583_8_Ctrl_Pats.RDS"),
                              PBMC12 = list(train_dataset="GSE96583_8_Stim_Pats.RDS",  test_dataset="GSE96583_batch1_3_samples.RDS"),
                              pancreas1 = list(train_dataset="Muraro_pancreas_clean.RDS",  test_dataset="Segerstolpe_pancreas_clean.RDS"),
                              pancreas2 = list(train_dataset="Muraro_pancreas_clean.RDS",  test_dataset="Xin_pancreas_clean.RDS"),
                              pancreas3 = list(train_dataset="Segerstolpe_pancreas_clean.RDS",  test_dataset="Muraro_pancreas_clean.RDS"),
                              pancreas4 = list(train_dataset="Segerstolpe_pancreas_clean.RDS",  test_dataset="Xin_pancreas_clean.RDS"),
                              pancreas5 = list(train_dataset="Xin_pancreas_clean.RDS",  test_dataset="Muraro_pancreas_clean.RDS"),
                              pancreas6 = list(train_dataset="Xin_pancreas_clean.RDS",  test_dataset="Segerstolpe_pancreas_clean.RDS"),
                              ADASD1 = list(train_dataset="ADASD_AD.RDS",  test_dataset="ADASD_autism.RDS"),
                              ADASD2 = list(train_dataset="ADASD_autism.RDS",  test_dataset="ADASD_AD.RDS"),
                              midbrain1 = list(train_dataset="midbrain_human.RDS",  test_dataset="midbrain_mouse.RDS"),
                              midbrain2 = list(train_dataset="midbrain_mouse.RDS",  test_dataset="midbrain_human.RDS")
                              # cellbench=list(tenx="cellbench_10x.RDS",CELseq2="cellbench_CELseq2.RDS",Dropseq="cellbench_Dropseq.RDS")
                            )


dataset.properties <- list(
                          PBMC_AllCells_withLabels=list(sample_threshold=5,gene_name_type="SYMBOL",total_num=61309),
                          GSE96583_batch1_3_samples=list(sample_threshold=5,gene_name_type="SYMBOL",cell_type_map="PBMC_celltype_map.csv",total_num=14017),
                          GSE96583_8_Stim_Pats=list(sample_threshold=5,gene_name_type="SYMBOL",cell_type_map="PBMC_celltype_map.csv",total_num=14126),
                          GSE96583_8_Ctrl_Pats=list(sample_threshold=5,gene_name_type="SYMBOL",cell_type_map="PBMC_celltype_map.csv",total_num=14617),
                          ADASD_AD=list(sample_threshold=5,gene_name_type="SYMBOL",total_num=70634),
                          ADASD_autism=list(sample_threshold=5,gene_name_type="SYMBOL",total_num=104559),
                          Muraro_pancreas=list(sample_threshold=3,gene_name_type="SYMBOL",total_num=2126,
                                                              cell_types=c('beta','alpha','delta',
                                                                          'acinar','ductal','gamma',
                                                                          'endothelial')),
                          Segerstolpe_pancreas=list(sample_threshold=3,gene_name_type="SYMBOL",total_num=3514,
                                               cell_types=c('beta','alpha','delta',
                                                            'acinar','ductal','gamma',
                                                            'endothelial')),
                          Xin_pancreas=list(sample_threshold=3,gene_name_type="SYMBOL",total_num=1600,
                                               cell_types=c('beta','alpha','delta',
                                                            'acinar','ductal','gamma',
                                                            'endothelial')),
                          midbrain_human=list(sample_threshold=3,total_num=1695,gene_name_type="SYMBOL"),
                          midbrain_mouse=list(sample_threshold=3,total_num=1518,gene_name_type="SYMBOL"),
                          cellbench_10x=list(sample_threshold=3,gene_name_type="ENSEMBL"),
                          cellbench_CELseq2=list(sample_threshold=3,gene_name_type="ENSEMBL"),
                          cellbench_Dropseq=list(sample_threshold=3,gene_name_type="ENSEMBL")
                      )

