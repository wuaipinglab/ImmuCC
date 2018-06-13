
#!/bin/bash

# The scripts in "RNASeq_pipeline.sh" are used to preprocess the RNA-Seq data(from "fastq" to "reads"). The scripts including four major steps: 
# quality control; mapping; bam file sorting; quantification.
# After that, the scripts in "HTSeq_stat.R" were used for merging the expression of immune receptor gene families 
# and data normalization.

# For more detail about RNA-Seq data analysis, please read our paper:
# Chen Z, Quan L, Huang A, Zhao Q, Yuan Y, Yuan X, Shen Q, Shang J, Ben Y, Qin FX-F and Wu A (2018) 
# seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment
# From Mouse RNA-Seq Data. Front. Immunol. 9:1286. doi: 10.3389/fimmu.2018.01286.

###############################################################################
#                          input parameters
###############################################################################

## Directory of the base work directory
base_dir=$1

# library_layout
library_layout=$2

# Directory to the software
software_path=$3

# Directory to the reference
ref_path=4

# Directory to the scripts
script_path=5

# number of threads
thread=$6

###############################################################################
gtf=${ref_path}/New_Mus_musculus.GRCm38.83.gtf
genome=${ref_path}/Mus_musculus.GRCm38.dna.primary_assembly.83.fa
star_ref=${ref_path}/star
RSeQC_ref=${ref_path}/GRCm38_mm10_Ensembl.bed

trimmomatic=${software_path}/Trimmomatic-0.35/trimmomatic-0.35.jar
trimmomatic_adapter=${software_path}/Trimmomatic-0.35/adapters/TruSeq3-PE.fa
trimmomatic_software_para="2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10"
fastqc=${software_path}/FastQC/fastqc
star=${software_path}/STAR/bin/Linux_x86_64/STAR
samtools=${software_path}/samtools-1.3.1/samtools
infer_experiment=${software_path}/RSeQC-2.6.3/scripts/infer_experiment.py
htseq=${software_path}/htseq-count


if [ ! -d "${base_dir}/01fastq" ]; then
mkdir ${base_dir}/01fastq
fi

if [ ! -d "${base_dir}/02trimmed" ]; then
mkdir ${base_dir}/02trimmed
fi

if [ ! -d "${base_dir}/03mapping" ]; then
mkdir ${base_dir}/03mapping
fi

if [ ! -d "${base_dir}/04sorted" ]; then
mkdir ${base_dir}/04sorted
fi

if [ ! -d "${base_dir}/05htseq" ]; then
mkdir ${base_dir}/05htseq
fi

if [ ! -d "${base_dir}/raw_fastqc" ]; then
mkdir ${base_dir}/raw_fastqc
fi

if [ ! -d "${base_dir}/new_fastqc" ]; then
mkdir ${base_dir}/new_fastqc
fi


#####################################################################################################
# 
if [ ! -d "${star_ref}" ]; then
mkdir ${star_ref}
$star --runThreadN ${thread} \
      --runMode genomeGenerate \
      --genomeDir ${star_ref} \
      --genomeFastaFiles ${genome} \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang 100
fi

## gunzip the *.gz files
echo "quality control"
# sh 01gunzip.sh ${input_path}/01fastq ${result_path}/01fastq

## gunzip the files
echo "SRA to fastq\n"
# sh ${script_dir}/00sra_fastq.sh ${input_path}/00sra ${result_path}/01fastq

## quality control
echo "quality control\n"
sh ${script_dir}/02qc.sh ${library_layout} ${base_dir}/01fastq ${base_dir}/02trimmed ${fastqc} ${trimmomatic} ${trimmomatic_adapter} ${base_dir}/raw_fastqc ${base_dir}/new_fastqc ${thread}

## mapping
echo "mapping"
sh ${script_dir}/03mapping.sh ${library_layout} ${base_dir}/02trimmed ${base_dir}/03mapping ${star} ${star_ref} ${thread}

## sorting
echo "Bam sorting"
sh ${script_dir}/04samtools.sh ${samtools} ${base_dir}/03mapping ${base_dir}/04sorted ${thread}

## Get the strandness information
echo "strandness"
#sh ${script_dir}/05-1.strand.sh ${base_dir}/04sorted  ${base_dir}/StrandInfo ${infer_experiment} ${RSeQC_ref}
Rscript ${script_dir}/05-2.RSEQc.stat.R ${base_dir}/StrandInfo ${base_dir}/strand.txt

## quantification with HTSeq
echo "HTSeq" 
sh ${script_dir}/06htseq.sh ${base_dir}/04sorted ${base_dir}/05htseq ${base_dir}/strand.txt ${htseq} ${gtf} 

# quantification with RSEM
#echo "RSEM"
#sh ${script_dir}/07rsem.sh ${base_dir}/03mapping ${base_dir}/06rsem ${base_dir}/strand.txt ${rsem} ${rsem_ref} ${thread}

## immune receptor gene merging and normalization
echo "normalization"
Rscript ${script_dir}/MouseHTSeq_counts_stat.R ${script_path} ${result_path}/06htseq ${filename}



