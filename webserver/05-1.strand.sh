
# Directory of the sorted bam files
input_path=$1

## out file name of strand information
out_name=$2

## Directory to the infer_experiment in RSeQC
infer_experiment=$3

# Directory to the RSeQC refernce bed file
RSeQC_ref=$4

#####################################################################################################
for i in $(ls ${input_path}/*bam)
do
sample_name=`basename $i|sed 's/.bam//'`
echo ${sample_name}
python ${infer_experiment} -r ${RSeQC_ref} -i ${input_path}/${sample_name}.bam

done >${out_name}
