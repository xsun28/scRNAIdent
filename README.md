# scRNAIdent

## Description

*__scRNAIdent__* provides a modularized pipeline tool for automating the evaluation and comparison of cell typing methods in scRNA-seq analysis. The The workflow of the pipeline is illustrated in the picture below.

![workflow image](https://github.com/xsun28/scRNAIdent/blob/master/workflow.png?raw=true)

## Pipeline modules
In brevity, the pipeline consists of seven modules. These modules are functionally independent and can be easily modified for adding new cell typing methods or evaluation metrics. The modules are described in the following subsections.
### 1. Configuration module
The instructions to drive modules and communications between modules are achieved via multiple configuration file and objects. 

### 2. Data preprocessor
### 3. Experiment drive
### 4. Experimental dataset constructor
### 5. Method drive
### 6. Analysis module
### 7. Reporting module

## Examples
To run the experiment of amount of cells

1. Specifying the experimental parameters in the *exp_config.R* 
	* Set experiment name
	
			experiment <- "cell_number"  
		
	* Set experimental initiate parameters
		