This repository contains a brief introduction on how to run the tissue specific deconvolution model
======================================================================================================

1.Dependencies
------------------------------------------------------------------------------------------------------
1.1 Software
---------
 (1).	`R` -- Install R from https://cran.r-project.org/.
 (2).	Once installed, open a terminal and at the command prompt, type R.
 (3).	At the R prompt: Install the following R packages by issuing command:
      `install.packages(c(“e1071”, “preprocessCore”))`  

 You can run the following commands to see whether it has been successfully installed.                       
 `library(e1071)`       
 `library(preprocessCor)`
 
1.2 Scripts and training data
------------------------------------------------------------------------------------------------------
 1.	Scirpts for `ImmuCC` can be accessed at: https://github.com/wuaipinglab/ImmuCC/blob/master/Microarray_Deconvolution.R
 2.	Scirpts for `CIBERSORT.R` can be accessed from https://cibersort.stanford.edu/ upon an request from `CIBERSORT` team.
 3.	Tissue specific signature matrix can be downloaded at: https://github.com/wuaipinglab/ImmuCC/tree/master/tissue_immucc/SignatureMatrix


2 How to estimate tissue immune cell proportion
-----------------------------------------------------------------------------------------------------
2.1 Preprocess the raw RNA-Seq data
------------------------------------
 Methods used to preprocess the bulk RNA-Seq data are the same as what we used in `seq_ImmuCC`. The shell scripts on how to transform fastq format sequencing data into the expression matrix are available at https://github.com/wuaipinglab/ImmuCC/tree/master/webserver

2.2 Estimate tissue immune cell proportion with tissue specific model
------------------------------------------------------------------------------------------------------
 E.g. when estimating the immune cell constitution from the transcriptome of lung, the lung specific signature matrix `”Lung.sig.matrix.csv”` is used in parameter `training_data`.

`Immune.proportion <- ImmuCC (expression, training_data = ”Lung.sig.matrix.csv”)`

 #expression: matrix of sample expression profile;
 traing_data：training signature matrix;


3.Output result
--------------------------------------------------------------------------------------------------------
 The calculated result will be save into a csv or txt format file. the row represents samples and the column was correspondents to the predicted immune cells. 
