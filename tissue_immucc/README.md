This repository contains a brief introduction on how to run the tissue specific computational tool
======================================================================================================

1.Dependencies
------------------------------------------------------------------------------------------------------
1.1 Software
---------
 (1).	`R` -- Install R from https://cran.r-project.org/.
 (2).	Once installed, open a terminal and at the command prompt, type R.
 (3).	At the R prompt: Install the following R packages by issuing command:
      `install.packages(c(“e1071”, “preprocessCore”))`  

 After that you can run the following commands to see whether it has been successfully installed.                       
 `library(e1071)`       
 `library(preprocessCor)`
 
1.2 Code
------------------------------------------------------------------------------------------------------
 1.	Scirpts of `ImmuCC` can be accessed at: https://github.com/wuaipinglab/ImmuCC/blob/master/Microarray_Deconvolution.R
 2.	Scirpts of `CIBERSORT.R` can be accessed from https://cibersort.stanford.edu/ upon an request from `CIBERSORT` team.
 3.	Tissue specific signature matrix can be downloaded at: https://github.com/wuaipinglab/ImmuCC/tree/master/tissue_immucc/SignatureMatrix


2 How to estimate tissue immune cell proportion
-----------------------------------------------------------------------------------------------------
2.1 Preprocess the raw RNA-Seq data
------------------------------------
 Methods to preprocess the bulk RNA-Seq data were the same as what we used in `seq_ImmuCC`. The shell scripts on how to  transform the fastq format sequencing data into the expression matrix were available at https://github.com/wuaipinglab/ImmuCC/tree/master/webserver

2.2 Estimate tissue immune cell proportion with tissue specific model
------------------------------------------------------------------------------------------------------
 E.g. when estimating the relative proportion between different immune cells from the transcriptome of lung, the lung specific signature matrix was used.

`Immune.proportion <- ImmuCC (expression, training_data = ”Lung.sig.matrix.csv”)`

 #expression: transcriptome profile of biological sample;
 traing_data： training signature matrix;


3.Output result
--------------------------------------------------------------------------------------------------------
 The calculated result will be save into a csv or txt format file. Each row of result table represents each biological sample and the column name correspondents to each immune cell. 
