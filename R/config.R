home <- '/Users/Xiaobo/git/scRNAIdent/'
conda_home <- "/Users/Xiaobo/opt/anaconda3/bin/conda"
data_home <- paste(home, 'data/',sep='')
dir.create(data_home)
result_home <- paste(home,'results/',sep='')
marker_home <- paste(home,'markers/',sep='')
type_home <- str_glue("{home}types")
pretrained_classifier_home <- str_glue("{home}pretrained")
raw_data_home <- paste('/Volumes/Seagate/jobs/Dropbox/CelltypeIdentification/', 'data/',sep='')


