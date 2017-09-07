# PE for Paired end sequencing and SE for Singlen end sequencing
library_layout=$1

# Directory of the RNASeq result
result_path=$2

## Directory to the software
software=$3

# number of threads
thread=$4

adapter=${software}/Trimmomatic-0.35/adapters/TruSeq3-PE.fa
software_para="2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10"
fastqc=${software}/FastQC/fastqc
trim=${software}/Trimmomatic-0.35/trimmomatic-0.35.jar

#####################################################################################################
if [ "$library_layout" = "PE" ];then

fastq_path=${result_path}/01fastq
trimm_path=${result_path}/02trimmed

for i in $(ls ${fastq_path}/*_R1.fastq)
do
##sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/_R1.fastq//'`
echo ${sample_name}

## Quality stat of raw samples
$fastqc ${fastq_path}/${sample_name}_R1.fastq -t ${thread} -o ${result_path}/raw_fastqc 
$fastqc ${fastq_path}/${sample_name}_R2.fastq -t ${thread} -o ${result_path}/raw_fastqc

## Remove low quality reads and trimm adapters
java -jar $trim PE -threads ${thread} ${fastq_path}/${sample_name}_R1.fastq ${fastq_path}/${sample_name}_R2.fastq -baseout ${trimm_path}/${sample_name}_trimm ILLUMINACLIP:$adapter:${software_para}

## Quality stat after QC
$fastqc ${trimm_path}/${sample_name}_trimm_1P -t ${thread} -o ${result_path}/new_fastqc 
$fastqc ${trimm_path}/${sample_name}_trimm_2P -t ${thread} -o ${result_path}/new_fastqc
done

#####################################################################################################
elif [ "$library_layout" = "SE" ];then
fastq_path=${result_path}/01fastq
trimm_path=${result_path}/02trimmed

for i in $(ls ${fastq_path}/*.fastq)
do
#sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/.fastq//'`
echo ${sample_name}
echo "Quality Control begin"

## Quality stat of raw samples
$fastqc ${fastq_path}/${sample_name}.fastq -t ${thread} -o ${result_path}/raw_fastqc 

## Remove low quality reads and trimm adapters
java -jar $trim SE -threads ${thread} ${fastq_path}/${sample_name}.fastq ${trimm_path}/${sample_name}_trimm ILLUMINACLIP:$adapter:${software_para}

## Quality stat after QC
$fastqc ${trimm_path}/${sample_name}_trimm -t ${thread} -o ${result_path}/new_fastqc

done
echo "RNASeq finished"

#####################################################################################################
else
cd ~
echo "error" > error.log
fi
