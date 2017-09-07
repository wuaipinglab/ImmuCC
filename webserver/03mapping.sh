
# PE for Paired end sequencing and SE for Singlen end sequencing
library_layout=$1

# Directory of the RNASeq result
result_path=$2

## Directory to the software
software=$3

# Directory to the refernce file
ref=$4

# number of threads
thread=$5

star=${software}/STAR/bin/Linux_x86_64/STAR
star_ref=${ref}/mm10/star_v2

gtf=${ref}/mm10/New_Mus_musculus.GRCm38.83.V2.gtf
genome=${ref}/mm10/Mus_musculus.GRCm38.dna.primary_assembly.83.fa

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

#####################################################################################################
## Mapping
if [ "$library_layout" = "PE" ];then
trimm_path=${result_path}/02trimmed
mapping_path=${result_path}/03mapping

for i in $(ls ${trimm_path}/*trimm_1P)
do
sample_name=`basename $i|sed 's/_trimm_1P//'`
echo ${sample_name}
mkdir ${mapping_path}/${sample_name}
$star --runMode alignReads --genomeDir ${star_ref} --runThreadN ${thread} --readFilesIn ${trimm_path}/${sample_name}_trimm_1P ${trimm_path}/${sample_name}_trimm_2P --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --limitBAMsortRAM 99000000000000 --outSAMunmapped Within --outFilterMultimapNmax 1000 --outFilterMismatchNmax 999 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix ${mapping_path}/${sample_name}/${sample_name}
done

#####################################################################################################
elif [ "$library_layout" = "SE" ];then
trimm_path=${result_path}/02trimmed
mapping_path=${result_path}/03mapping

for i in $(ls ${trimm_path}/*.fastq)
do
sample_name=`basename $i|sed 's/_trimm//'`
echo ${sample_name}
mkdir ${mapping_path}/${sample_name}
$star --runMode alignReads --genomeDir ${star_ref} --runThreadN ${thread} --readFilesIn ${trimm_path}/${sample_name}_trimm --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --limitBAMsortRAM 99000000000 --outSAMunmapped Within --outFilterMultimapNmax 1000 --outFilterMismatchNmax 999 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix ${mapping_path}/${sample_name}/${sample_name}
done

#####################################################################################################
else
cd ~
echo "error" > error.log
fi


