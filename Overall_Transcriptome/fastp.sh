#!/bin/bash
#SBATCH --account=def-dogfish
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=03:00:00
#SBATCH --mem=15000M
#SBATCH --array=1-192

module load StdEnv/2020
module load fastp/0.20.1

list=/home/biomatt/scratch/lkst_transcriptomes/raw_fastq_list.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')
set -- $var

forward=$1
reverse=$2

echo ${forward}
echo ${reverse}

sample=$(basename ${forward} -1.fastq.gz | cut -d "/" -f5)

echo ${sample}

fastp -i ${forward} -I ${reverse} -o /home/biomatt/scratch/lkst_transcriptomes/trimmed_data/${sample}.R1.fq.gz -O /home/biomatt/scratch/lkst_transcriptomes/trimmed_data/${sample}.R2.fq.gz -j /home/biomatt/scratch/lkst_transcriptomes/trimmed_data/${sample}_fastp.json -h /home/biomatt/scratch/lkst_transcriptomes/trimmed_data/${sample}_fastp.html --detect_adapter_for_pe
