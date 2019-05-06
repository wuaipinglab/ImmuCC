################################################
# Sort the mapped reads according to their name
################################################

# Directory to the software
samtools=$1

# Directory of the RNASeq result
input_path=$2

# Directory of the RNASeq result
output_path=$3

# number of threads
thread=$4

for i in $(ls -d ${input_path}/*)
do
samplename=`basename $i`
$samtools sort -@ $thread -n -o ${output_path}/${samplename}.bam ${input_path}/${samplename}/${samplename}Aligned.out.bam
done
