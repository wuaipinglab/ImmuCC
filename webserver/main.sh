
# Directory to the reference data
ref=$1

# Directory to the software
software=$2

## Directory of the fastq files
input_path=$3

# Directory of the RNASeq result
result_path=$4

# name of the expression file
filename=$5

# Number of threads used
thread=$6

# Directory to the scripts
script_path=$7

cd ${script_path}

if [ ! -d "${result_path}/01fastq"]; then
mkdir ${result_path}/01fastq
fi

if [ ! -d "${result_path}/02trimmed" ]; then
mkdir ${result_path}/02trimmed
fi

if [ ! -d "${result_path}/03mapping" ]; then
mkdir ${result_path}/03mapping
fi

if [ ! -d "${result_path}/04sorted" ]; then
mkdir ${result_path}/04sorted
fi

if [ ! -d "${result_path}/05htseq" ]; then
mkdir ${result_path}/05htseq
fi

if [ ! -d "${result_path}/raw_fastqc" ]; then
mkdir ${result_path}/raw_fastqc
fi

if [ ! -d "${result_path}/new_fastqc" ]; then
mkdir ${result_path}/new_fastqc
fi

## quality control
sh ${script_path}/02qc.sh PE ${result_path} ${software} $thread

## mapping
sh ${script_path}/03mapping.sh PE ${result_path} ${software} ${ref} $thread

## sorting
sh ${script_path}/04samtools.sh ${result_path} ${software} $thread

## strandness
sh ${script_path}/05-1.strand.sh ${result_path} ${result_path}/StrandInfo ${software} ${ref}
Rscript ${script_path}/05-2.RSEQc.stat.R ${result_path}/StrandInfo ${result_path}/strand.txt

## quantification
sh ${script_path}/06htseq ${result_path} ${result_path}/strand.txt ${software} ${ref}

## normalization
Rscript ${script_path}/MouseHTSeq_counts_stat.R ${result_path}/05htseq ${result_path} ${filename}


