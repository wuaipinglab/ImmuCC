---
title: "A basic workflow for tissue specific deconvolution model"
author: "Ziyi Chen"
output: html_document

abstract:
Estimate the relative proportion of tissue immune cell from tissue transcriptome with a tissue specific model.
---


# Introduction
The tissue specific computational model can be used to predict the relative immune cell proportion with a series of tissue specific training signature matrix.

It is achieved by two major steps:
  *  1. Raw bulk RNA-Seq data preprocessing
  *  2. Immune cell proportion prediction

------------------------------------------------------------------------------------------------------------------------
## 1. Raw bulk RNA-Seq data preprocessing
Raw fastq format sequencing data should be first preprocessed into the expression matrix.
Here, six major steps are included in this process. 
  (1)	Quality control. Removing low quality sequence and adaptor trimming et.al.
  (2)	Mapping. Mapping sequencing reads to the reference genome with STAR;
  (3)	Sorting. Sort the mapped reads according to their name;
  (4)	Strand. Get the strand information of the library;
  (5)	Quantification. Quantify gene expression with HTSeq;
  (6)	Normalization. Subtypes of T cell receptor and B cell receptor genes are merged into major families and the raw expression matrix are normalized with the quantile normalization method.

Shell scripts on how to preprocess the raw RNA-Seq data can be obtained from https://github.com/wuaipinglab/ImmuCC/blob/master/webserver/RNASeq_pipeline.sh.

* Command for data preprocessing is listed below:
  >`sh RNASeq_pipeline.sh ${Directory to the base} PE ${Directory to the software} ${Directory to the reference} ${Directory to the scripts} ${thread_number}`

* Description for the arguments in this command:

   `${Directory to the base}`: This directory contains 7 files, namely, `01fastq, 02trimmed, 03mapping, 04sorted, 05htseq, raw_fastqc, new_fastqc`. 
  
   `${Directory to the software}`: This directory contains all softwares used in this analysis including: `FastQC, STAR, samtools, RSeQC, htseq-count, R.`
  
   `${Directory to the reference}`: This directory contains the reference data required for mapping and quantification. In my analysis, the following reference data including: `Mus_musculus.GRCm38.83.gtf, Mus_musculus.GRCm38.dna.primary_assembly.83.fa, GRCm38_mm10_Ensembl.bed` are used. You can also download the lattest reference data.
  
   `${Directory to the scripts}`: All basic scripts used above including: `02qc.sh, 03mapping.sh, 04samtools.sh, 05-1.strand.sh, 05-2.RSEQc.stat.R, 06htseq.sh, MouseHTSeq_counts_stat.R` are put in this directory.
   
   `${thread_number}`: Number of thread to be used.
------------------------------------------------------------------------------------------------------------------------
## 2. Immune cell proportion prediction
E.g. To estimate immune cell proportion from the transcriptome data of lung, the lung specific signature matrix is used in parameter `training_data`. As the transcriptome data has been normaliazed during data preprocessing steps, it will be not necessary to normalize it again.

Function `ImmuCC` is used to calculated immune cell proportions from the expression matrix of samples.
Two basic arguments including `expression` and `traing_data` are needed for this function.
> expression: expression matrix of the biological sample;

> traing_data: tissue specific signature matrix to be used.

An examble on how to run this function was listed below:
> Immune.proportion <- ImmuCC (expression, ”Lung.sig.matrix.csv”)
