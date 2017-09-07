

# Directory of the RNASeq result
result_path=$1

## Directory to the software
software=$2

# number of threads
thread=$3

samtools=${software}/samtools-1.3.1/samtools

for i in $(ls -d ${result_path}/03mapping/*)
do
samplename=`basename $i`
$samtools sort -@ $thread -n -o ${result_path}/04sorted/${samplename}.bam ${result_path}/03mapping/${samplename}/${samplename}Aligned.sortedByCoord.out.bam
done

