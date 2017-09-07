
# Directory of the RNASeq result
result_path=$1

## out file name
strand_file=$2

## Directory to the software
software=$3

# Directory to the refernce file
ref=$4

gtf=${ref}/mm10/New_Mus_musculus.GRCm38.83.V2.gtf


#####################################################################################################
for i in $(ls ${result_path}/04sorted/*.bam)
do
sample_name=`basename $i | sed 's/.bam//'`
strand=`cat $strand_file|grep $sample_name|awk -F" " '{print $4}'`

python ${software}/htseq-count -f bam -t gene -m union -s $strand ${result_path}/04sorted/${sample_name}.bam ${gtf} >${result_path}/06htseq/${sample_name}.txt
done  

