---
title: "A basic workflow for tissue specific deconvolution model"
author: "Ziyi Chen"
output: html_document

abstract:
Estimate the relative proportion of tissue immune cell from tissue transcriptome with a tissue specific model.
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Introduction
The tissue specific computational model can be used to predict the relative proportion of immune cell with a series of tissue specific training signature matrix.

It is achieved by three major steps:
  1. Raw bulk RNA-Seq data preprocessing
  2. Immune cell proportion prediction

------------------------------------------------------------------------------------------------------------------------
1. Raw bulk RNA-Seq data preprocessing
Raw fastq format sequencing data should be first preprocessed into the expression matrix.
Here, six major steps are included to preprocess the raw data. 
(1)	Quality control;
(2)	Mapping. Mapping sequencing reads to the reference genome with STAR;
(3)	Sorting. Sort the mapped reads according to their name;
(4)	Strand. Get the strand information of the library;
(5)	Quantification. Quantify Gene expression with HTSeq.
(6)	Normalization. Subtypes of T cell receptor and B cell receptor genes are merged into major families and the raw expression matrix are normalized with the quantile normalization method.

The shell scripts on how to preprocess it can be obtained from https://github.com/wuaipinglab/ImmuCC/blob/master/webserver/RNASeq_pipeline.sh.

Commands for data preprocessing are listed below:
`sh RNASeq_pipeline.sh ${Directory to the base} PE ${Directory to the software} ${Directory to the reference} ${Directory to the scripts} 24`
  Description for the arguments in this command:
  `${Directory to the base}`: This directory contains 7 files, namely, `01fastq, 02trimmed, 03mapping, 04sorted, 05htseq, raw_fastqc, new_fastqc`. 
  `${Directory to the software}`: Softwares used here including: `FastQC, STAR, samtools, RSeQC, htseq-count, R.` were all put in this directory.
  `${Directory to the reference}`: This directory contains all reference data used. In my analysis, the following reference data. including: `Mus_musculus.GRCm38.83.gtf, Mus_musculus.GRCm38.dna.primary_assembly.83.fa, GRCm38_mm10_Ensembl.bed` are used. You can download the lattest version as you want.
  `${Directory to the scripts}`: All scripts used here including: `02qc.sh, 03mapping.sh, 04samtools.sh, 05-1.strand.sh, 05-2.RSEQc.stat.R, 06htseq.sh, MouseHTSeq_counts_stat.R` are put in this directory.

------------------------------------------------------------------------------------------------------------------------
2. Immune cell proportion prediction
E.g. To estimate immune cell proportion from the transcriptome data of lung, the lung specific signature matrix is used as the training data. As the sample transcriptome data have been normaliazed during data preprocessing steps, it will be not necessary to normalize it again.

Function `ImmuCC` was used to calculated immune cell proportions from the expression matrix of sample.
Two basic arguments including `expression` and `traing_data` were needed for this function 
# expression: expression matrix of the biological sample
# traing_data: tissue specific signature matrix to be used

An examble on how to run this function was listed below
   Immune.proportion <- ImmuCC (expression, training_data = ”Lung.sig.matrix.csv”)
