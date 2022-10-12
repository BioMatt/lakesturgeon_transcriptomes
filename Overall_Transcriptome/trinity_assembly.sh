#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load bowtie2/2.4.4 trinity/2.14.0 samtools/1.13 jellyfish/2.3.0 salmon/1.7.0
module load scipy-stack
pip list | grep -a numpy
module load java/14.0.2

Trinity --seqType fq --max_memory 240G --FORCE --bflyHeapSpaceInit 1M --CPU 8 --min_kmer_cov 2 --normalize_by_read_set --no_bowtie --samples_file /home/biomatt/scratch/lkst_transcriptomes/trimmed_fastq_list_reduced.txt --output /home/biomatt/scratch/lkst_transcriptomes/trinity_tissuespecific_outputs
