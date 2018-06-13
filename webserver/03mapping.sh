
# PE for Paired end sequencing and SE for Singlen end sequencing
library_layout=$1

# Directory of the trimmed fastq file
input_path=$2

# Directory of the mapping result
output_path=$3

## Directory to STAR software
star=$4

# Directory to the STAR refernce file
star_ref=$5

# number of threads
thread=$6

trimm_path=${input_path}
mapping_path=${output_path}

#####################################################################################################
## Mapping
if [ "$library_layout" = "PE" ];then
for i in $(ls ${input_path}/*trimm_1P)
do
sample_name=`basename $i|sed 's/_trimm_1P//'`
echo ${sample_name}
mkdir ${output_path}/${sample_name}
$star --runMode alignReads \
      --genomeDir ${star_ref} \
      --runThreadN ${thread} \
      --readFilesIn ${input_path}/${sample_name}_trimm_1P ${input_path}/${sample_name}_trimm_2P \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --limitBAMsortRAM 99000000000000 \
      --outSAMunmapped Within \
      --outFilterMultimapNmax 1000 \
      --outFilterMismatchNmax 999 \
      --outSAMtype BAM Unsorted \
      --quantMode TranscriptomeSAM \
      --outFileNamePrefix ${output_path}/${sample_name}/${sample_name}
done

#####################################################################################################
elif [ "$library_layout" = "SE" ];then
for i in $(ls ${input_path}/*trimm)
do
sample_name=`basename $i|sed 's/_trimm//'`
echo ${sample_name}
mkdir ${output_path}/${sample_name}
$star --runMode alignReads \
      --genomeDir ${star_ref} \
      --runThreadN ${thread} \
      --readFilesIn ${input_path}/${sample_name}_trimm \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --limitBAMsortRAM 99000000000 \
      --outSAMunmapped Within \
      --outFilterMultimapNmax 1000 \
      --outFilterMismatchNmax 999 \
      --outSAMtype BAM Unsorted \
      --quantMode TranscriptomeSAM \
      --outFileNamePrefix ${output_path}/${sample_name}/${sample_name}
done
echo "mapping finished!"
#####################################################################################################
else
cd ~
echo "error" > error.log
fi


