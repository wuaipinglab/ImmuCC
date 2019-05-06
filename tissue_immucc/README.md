This repository contains the dependencies to the tissue specific computational tool
=============

1.Dependencies
--------------
##software
---------
1.	R -- Install R from https://cran.r-project.org/.
2.	Once installed, open a terminal and at the command prompt, type R.
3.	At the R prompt: Install the following R packages by issuing command:
install.packages(c(“e1071”, “preprocessCore”))  

After that you can run the following commands to see whether it has been successfully installed.                       
 library(`e1071`)        
 library(`preprocessCore`) 
 
##Codes
----------------
1.	Download the code repository at: https://github.com/wuaipinglab/ImmuCC/blob/master/Microarray_Deconvolution.R
2.	Download the tissue specific signature matrix at: https://github.com/wuaipinglab/ImmuCC/tree/master/tissue_immucc/SignatureMatrix
3.	Access the scirpts of `CIBERSORT.R` from https://cibersort.stanford.edu/ upon an request from `CIBERSORT` team.

2.preprocess the raw RNA-Seq data
------------------------------------
Methods to preprocess the bulk RNA-Seq data were the same as what we used in `seq_ImmuCC`. The scripts on how to transform the fastq format sequencing data into the expression matrix were listed in https://github.com/wuaipinglab/ImmuCC/tree/master/webserver

3.estimate tissue immune cell proportion with tissue specific model
------------------------------------------------------------------------------
E.g. To estimate the immune cell proportion from the transcriptome of lung, the lung specific signature matrix was used.
`Immune.proportion <- ImmuCC (expression, training_data = ”Lung.sig.matrix.csv”)`
`expression`: transcriptome profile of biological sample
`traing_data`： training signature matrix

4.Output for the calculated result
--------
csv or txt format file will be created by function “ImmuCC” . Each row of result table represents each biological sample and the column name correspondents to each immune cell. 
