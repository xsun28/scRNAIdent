home <- '/home/xsun/git/scRNAIdent/'
conda_home <- "/usr/local/anaconda3/condabin/conda"
data_home <- paste(home, 'data/',sep='')
dir.create(data_home)
result_home <- paste(home,'results/',sep='')
marker_home <- paste(home,'markers/',sep='')
type_home <- str_glue("{home}types")
pretrained_classifier_home <- str_glue("{home}pretrained")
raw_data_home <- paste('/volume/scRNA/', 'data/',sep='')


