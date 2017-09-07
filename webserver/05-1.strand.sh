
# Directory of the RNASeq result
result_path=$1

## out file name of strand information
out_name=$2

## Directory to the software
software=$3

# Directory to the refernce file
ref=$4

ref=${ref}/GRCm38_mm10_Ensembl.bed
infer_experiment=${software}/RSeQC-2.6.3/scripts/infer_experiment.py


#####################################################################################################
for i in $(ls ${result_path}/04sorted/*bam)
do
sample_name=`basename $i|sed 's/.bam//'`
echo ${sample_name}
python ${infer_experiment} -r ${ref} -i ${result_path}/04sorted/${sample_name}.bam

done >${out_name}
