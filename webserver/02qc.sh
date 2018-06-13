

Note: For paired end sequencing data, the raw fastq files should be ended with *_1.fastq and *_2.fastq

# PE for Paired end sequencing and SE for Singlen end sequencing
library_layout=$1

# Directory of the raw fastq files
fastq_path=$2

# Directory of the trimmed fastq file
trimm_path=$3

# Directory to the software of fastqc
fastqc=$4

# Directory to the software of trimmomatic
trimmomatic=$5

# Directory to the trimmomatic adapters
adapter=$6

# Directory to the result of raw fastqc result
raw_fastqc=$7

# Directory to the result of new fastqc result
new_fastqc=$8

# number of threads
thread=$9

software_para="2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10"

#####################################################################################################
if [ "$library_layout" = "PE" ];then

for i in $(ls ${fastq_path}/*_1.fastq)
do
sample_name=`basename $i|awk -F"_" '{print $1}'`
#sample_end=`basename $i|awk -F"_" '{print $2}'`
#sample_name=`basename $i|sed 's/_1.fastq//'`
echo ${sample_name}

## Quality stat of raw samples
$fastqc ${fastq_path}/${sample_name}_1.fastq -t ${thread} -o ${raw_fastqc}
$fastqc ${fastq_path}/${sample_name}_2.fastq -t ${thread} -o ${raw_fastqc}

## Remove low quality reads and trimm adapters
java -jar $trimmomatic PE -threads ${thread} ${fastq_path}/${sample_name}_1.fastq ${fastq_path}/${sample_name}_2.fastq -baseout ${trimm_path}/${sample_name}_trimm ILLUMINACLIP:$adapter:${software_para}

## Quality stat after QC
$fastqc ${trimm_path}/${sample_name}_trimm_1P -t ${thread} -o ${new_fastqc} 
$fastqc ${trimm_path}/${sample_name}_trimm_2P -t ${thread} -o ${new_fastqc}
done

#####################################################################################################
elif [ "$library_layout" = "SE" ];then

for i in $(ls ${fastq_path}/*.fastq)
do
#sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/.fastq//'`
echo ${sample_name}
echo "Quality Control begin"

## Quality stat of raw samples
$fastqc ${fastq_path}/${sample_name}.fastq -t ${thread} -o ${raw_fastqc}

## Remove low quality reads and trimm adapters
java -jar $trimmomatic SE -threads ${thread} ${fastq_path}/${sample_name}.fastq ${trimm_path}/${sample_name}_trimm ILLUMINACLIP:$adapter:${software_para}

## Quality stat after QC
$fastqc ${trimm_path}/${sample_name}_trimm -t ${thread} -o ${new_fastqc}

done
echo "Quality control finished"

#####################################################################################################
else
cd ~
echo "error" > error.log
fi
