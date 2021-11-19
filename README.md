# scRNAIdent
## Description

*__scRNAIdent__* provides a modularized R pipeline tool for automating the evaluation and comparison of cell typing methods in scRNA-seq analysis. The The workflow of the pipeline is illustrated in the picture below.

![workflow image](https://github.com/xsun28/scRNAIdent/blob/master/workflow.jpg?raw=true)

## Pipeline modules
In brevity, the pipeline consists of seven modules. These modules are functionally independent and can be easily modified for adding new cell typing methods or evaluation metrics. The modules are described in the following subsections.
### 1. Configuration module
The configuration module consists of multiple configuration files. It specifies instructions for driving other modules and provides a means of communicating among modules. The components of the module is listed below.

| Name | File | Description |
|----|----|-------|
|Pipeline|*config.R*|Pipeline home path and various directories|
| Experiment|*experiments_config.R*| Parameters for driving experiment module. For example, experiment name, initial values experiment parameters, tested methods, used datasets, parameters of constructing experimental datasets, evaluation metrics, etc.|
|Method|*methods_config.R*|Methods hyperparameters|
|Dataset|*dataset_config.R*|Basic properties of datasets, lists of pairs of reference/query datasets|
|Marker gene|*markergene_config.R*|Parameters for marker gene selection method|


### 2. Data preprocessor
This module preprocesses raw scRNA-seq datasets into SingleCellExperiment objects, performs basic QC, and calculates and adds corresponding metadata.
### 3. Experiment drive
This module is the central part of the pipeline. It coordinates the functioning of other modules, sending instructions and receiving outputs from each other module. Its functionalities include:

- Parsing and updating configurations for experimental instructions
- Sending instructions to dataset constructor and receiving constructed training/testing datasets.
- Sending training/testing datasets to method drive. Receiving classification/clustering results from methods.
- Sending datasets and results to analysis module. Receiving reports of dataset properties metrics and performance metrics.
- Sending performance reports to reporting module for saving and plotting. 


### 4. Experimental dataset constructor
This module samples the reference/query datasets to construct the training/testing datasets based on received instructions from experiment drive.
### 5. Method drive
This module parses the method configurations and run the methods on received training/testing datasets based on specific experimental instructions from experiment drive.
### 6. Analysis module
This module receives training/testing datasets and classification/clustering results. It then calculates the analytical properties of the datasets, and performance metrics on the results.  
### 7. Reporting module
This module dumps the performance results and plots various analytical figures.
## Examples
To run the experiment of amount of cells

1. Specify the home, data and other relevant directories of the pipeline in *config.R*

		home <- '***'
		data_home <- '***'
		result_home <- '***'
		log_home <- '***'
		log_file <- "scRNAIdent.log"
		marker_home <- '***'
		type_home <- '***'
		pretrained_home <- '***'
		raw_data_home <- '***'

2. Specify the experimental parameters in the *experiments_config.R* 
	* Set experiment name

			experiment <- "cell_number"  

	* Set experimental methods and parameters 

			experiments.methods <- list(cell_number=list(...))
			experiments.parameters <- list(cell_number=list(...))

	* (optional) If adding new experiment, implement experimental parameters initiating and updating functions

			experiments.config.init.cell_number <- function(train_dataset, test_dataset=NULL, exp_config)
			experiments.config.update.cell_number <- function(train_dataset, test_dataset=NULL, exp_config)

3. (optional) If testing new datasets, add to and update the configurations in the *dataset_config.R*
4. (optional) If testing a new method, add the new method to *methods.R* and its configurations to *methods_config.R*
5. (optional) If testing new marker gene selection methods, add the methods to *marke_gene.R* and its configurations to *markergene_config.R*
6. (optional) If using a new evaluation metric, add the metric calculating function to *analysis.R* and add the metric to the metrics parameter in the corresponding experiment configurations in *experiments_config.R*
7. load the *experiments.R* script and run

		runExperiments(experiment)

