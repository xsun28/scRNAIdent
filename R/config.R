home <- '/Users/Xiaobo/git/scRNAIdent/'
data_home <- paste(home, 'data/',sep='')
result_home <- paste(home,'results/',sep='')
raw_data_home <- paste('/Users/Xiaobo/Dropbox/CelltypeIdentification/', 'data/',sep='')
raw_datasets <- list(PBMC="PBMC_AllCells_withLabels.Rdata",
                 pancreas=c("Muraro_pancreas.Rdata","Segerstolpe_pancreas.Rdata","Xin_pancreas.Rdata"))

datasets <- list(PBMC="PBMC_AllCells_withLabels.RDS",
                 pancreas=c("Muraro_pancreas_clean.RDS","Segerstolpe_pancreas_clean.RDS","Xin_pancreas_clean.RDS"))


