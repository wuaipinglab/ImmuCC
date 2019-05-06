This repository contains the dependencies on how to run 

Dependencies
=============

1.	R -- Install R from https://cran.r-project.org/.
2.	Once installed, open a terminal and at the command prompt, type R.
3.	At the R prompt: Install the following R packages by issuing command:
install.packages(c(“e1071”, “preprocessCore”, “”))  

After that you can run the following commands to see whether it has been successfully installed.                   
 library(mouse4302mmentrezgcdf)     
 library(e1071)        
 library(preprocessCore) 
 
How to preprocess your RNA-Seq data
Methods to preprocess the rna-Seq data were the same as what we used in Seq_ImmuCC. The scripts on how to transform the fastq format sequencing data into the expression matrix were listed in https://github.com/wuaipinglab/ImmuCC/tree/master/webserver
How to run code
1.	Download the code repository at: https://github.com/wuaipinglab/ImmuCC/blob/master/Microarray_Deconvolution.R
2.	Download the tissue specific signature matrix at: https://github.com/wuaipinglab/ImmuCC/tree/master/tissue_immucc/SignatureMatrix
3.	Access the scirpts of CIBERSORT.R from https://cibersort.stanford.edu/ upon an request from CIBERSORT team.
4.  The pipeline for preprocessing bulk RNA-Seq data were the same as the methods used in Seq-ImmuCC which was available at our web site
https://github.com/wuaipinglab/ImmuCC/tree/master/webserver. The fastq format transcriptome sequencing data can be analysed with it to produce the gene expression matrix of your sample.

How to estimate tissue immune cell proportion with this tissue specific model
E.g. To estimate the immune cell proportion from the transcriptome of lung, the lung specific signature matrix was used.
Immune.proportion <- ImmuCC (expression, training_data = ”Lung.sig.matrix.csv”)

Output
Csv or txt format file will be created by function “ImmuCC” . Each row of result table represents each biological sample and the column name correspondents to each immune cell. 