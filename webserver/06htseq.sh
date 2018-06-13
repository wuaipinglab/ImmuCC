# Directory of the sorted BAM file
input_path=$1

# Directory of the HTSeq RESULT
output_path=$2

## out file name
strand_file=$3

## Directory to the software
htseq=$4

# Directory to the refernce file
gtf=$5


#####################################################################################################
for i in $(ls ${input_path}/*.bam)
do
sample_name=`basename $i | sed 's/.bam//'`
strand=`cat $strand_file|grep $sample_name|awk -F" " '{print $4}'`
#echo "python $htseq -f bam -t gene -m union -s $strand ${input_path}/${sample_name}.bam ${gtf} >${output_path}/${sample_name}.txt"
python $htseq -f bam -t gene -m union -s $strand ${input_path}/${sample_name}.bam ${gtf} >${output_path}/${sample_name}.txt
done  
