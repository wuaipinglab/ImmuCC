
#!/bin/bash

# The scripts in "RNASeq_pipeline.sh" are used to preprocess the RNA-Seq data(from "fastq" to "reads"). The scripts including four major steps: 
# quality control; mapping; bam file sorting; quantification.
# After that, the scripts in "HTSeq_stat.R" were used for merging the expression of immune receptor gene families 
# and data normalization.

# For more detail about RNA-Seq data analysis, please read our paper:
# Chen Z, Quan L, Huang A, Zhao Q, Yuan Y, Yuan X, Shen Q, Shang J, Ben Y, Qin FX-F and Wu A (2018) 
# seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment
# From Mouse RNA-Seq Data. Front. Immunol. 9:1286. doi: 10.3389/fimmu.2018.01286.

## Directory to the software
software=""

# Directory to the refernce file
ref=""

adapter=${software}/Trimmomatic-0.35/adapters/TruSeq3-PE.fa
software_para="2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10"
gtf=${ref}/mm10/New_Mus_musculus.GRCm38.83.V2.gtf
genome=${ref}/mm10/Mus_musculus.GRCm38.dna.primary_assembly.83.fa 
star_ref=${ref}/mm10/star

fastqc=${software}/FastQC/fastqc
trim=${software}/Trimmomatic-0.35/trimmomatic-0.35.jar
star=${software}/STAR/bin/Linux_x86_64/STAR
samtools=${software}/samtools-1.3.1/samtools
htseq_count=${software}/htseq-count

#
library_layout=$1


input_path=$2

# whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse'
strand=$3

# number of threads
thread=24


if [ ! -d "${input_path}/02trimmed" ]; then  
　　mkdir ${input_path}/02trimmed  
fi

if [ ! -d "${input_path}/03mapping" ]; then  
　　mkdir ${input_path}/03mapping 
fi

if [ ! -d "${input_path}/04sorted" ]; then  
　　mkdir ${input_path}/04sorted 
fi

if [ ! -d "${input_path}/05htseq" ]; then  
　　mkdir ${input_path}/05htseq 
fi

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

############################################################################################################################### 
#                                                Paired end sequencing
###############################################################################################################################
if [ "$library_layout" = "PE" ];then

fastq_path=${input_path}/01fastq
trimm_path=${input_path}/02trimmed
mapping_path=${input_path}/03mapping
sorted_path=${input_path}/04sorted
htseq_path=${input_path}/05htseq

for i in $(ls ${fastq_path}/*_R1.fastq)
do

##sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/_R1.fastq//'`
echo ${sample_name}
echo "Quality Control begin"

## Quality stat of raw samples
$fastqc ${fastq_path}/${sample_name}_R1.fastq -t ${thread} -o ${input_path}/raw_fastqc 
$fastqc ${fastq_path}/${sample_name}_R2.fastq -t ${thread} -o ${input_path}/raw_fastqc

## Remove low quality reads and trimm adapters
java -jar $trim PE -threads ${thread} ${fastq_path}/${sample_name}_R1.fastq ${fastq_path}/${sample_name}_R2.fastq -baseout ${trimm_path}/${sample_name}_trimm ILLUMINACLIP:$adapter:${software_para}

## Quality stat after QC
$fastqc ${trimm_path}/${sample_name}_trimm_1P -t ${thread} -o ${input_path}/new_fastqc 
$fastqc ${trimm_path}/${sample_name}_trimm_2P -t ${thread} -o ${input_path}/new_fastqc

## Mapping
echo "Mapping"
mkdir ${mapping_path}/${sample_name}
$star --runMode alignReads 
      --genomeDir ${star_ref} 
      --runThreadN ${thread} 
      --readFilesIn ${trimm_path}/${sample_name}_trimm_1P ${trimm_path}/${sample_name}_trimm_2P 
      --alignIntronMin 20 
      --alignIntronMax 1000000 
      --alignMatesGapMax 1000000 
      --alignSJoverhangMin 8 
      --alignSJDBoverhangMin 1 
      --limitBAMsortRAM 99000000000000 
      --outSAMunmapped Within 
      --outFilterMultimapNmax 1000 
      --outFilterMismatchNmax 999 
      --outSAMtype BAM Unsorted 
      --quantMode TranscriptomeSAM 
      --outFileNamePrefix ${mapping_path}/${sample_name}/${sample_name}

## Sort the bam file by read name
echo "Samtools sorting"
$samtools sort -@ ${thread} -n -o ${sorted_path}/${samplename}.bam ${mapping_path}/${sample_name}/${sample_name}.bam

## Quantification with HTSeq
echo "Quantification"
python $htseq_count -f bam -t gene -m union -s $strand ${sorted_path}/${sample_name}.bam ${gtf} >${htseq_path}/${sample_name}.txt
done
echo "RNASeq finished"

############################################################################################################################### 
#                                                Single end sequencing
###############################################################################################################################
elif [ "$library_layout" = "SE" ];then
fastq_path=${input_path}/01fastq
trimm_path=${input_path}/02trimmed
mapping_path=${input_path}/03mapping
sorted_path=${input_path}/04sorted
htseq_path=${input_path}/05htseq

for i in $(ls ${fastq_path}/*.fastq)
do
#sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/.fastq//'`
echo ${sample_name}
echo "Quality Control begin"

## Quality stat of raw samples
$fastqc ${fastq_path}/${sample_name}.fastq -t ${thread} -o ${input_path}/raw_fastqc 

## Remove low quality reads and trimm adapters
java -jar $trim SE -threads ${thread} ${fastq_path}/${sample_name}.fastq ${trimm_path}/${sample_name}_trimm ILLUMINACLIP:$adapter:${software_para}

## Quality stat after QC
$fastqc ${trimm_path}/${sample_name}_trimm -t ${thread} -o ${input_path}/new_fastqc

## Mapping
echo "Mapping"

mkdir ${mapping_path}/${sample_name}
$star --runMode alignReads 
      --genomeDir ${star_ref} 
      --runThreadN ${thread} 
      --readFilesIn ${trimm_path}/${sample_name}_trimm 
      --alignIntronMin 20 
      --alignIntronMax 1000000 
      --alignMatesGapMax 1000000 
      --alignSJoverhangMin 8 
      --alignSJDBoverhangMin 1 
      --limitBAMsortRAM 99000000000 
      --outSAMunmapped Within 
      --outFilterMultimapNmax 1000 
      --outFilterMismatchNmax 999 
      --outSAMtype BAM Unsorted 
      --quantMode TranscriptomeSAM 
      --outFileNamePrefix ${mapping_path}/${sample_name}/${sample_name}

## Sort the bam file by read name
echo "Samtools sorting"
$samtools sort -@ ${thread} -n -o ${sorted_path}/${samplename}.bam ${mapping_path}/${sample_name}/${sample_name}.bam

## Quantification with HTSeq
echo "Quantification"
python $htseq_count -f bam -t gene -m union -s $strand ${sorted_path}/${sample_name}.bam ${gtf} >${htseq_path}/${sample_name}.txt

done
echo "RNASeq finished"

#####################################################################################################
else
cd ~
echo "error" > error.log
fi
