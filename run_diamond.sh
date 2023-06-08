#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=3:00:00
#SBATCH --mem=20000M
#SBATCH --array=1-14

module load StdEnv/2020
module load diamond/2.1.6

transcriptome_list=/home/biomatt/scratch/sterlet_blast/lkst_transcriptome_list.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${transcriptome_list}"
transcriptome=$($string)

set -- $transcriptome

echo ${transcriptome}

# Diamond database created from the sterlet protein fasta with diamond makedb --in GCF_010645085.2_ASM1064508v2_protein.faa -d sterlet_db
diamond blastx -q /home/biomatt/scratch/sterlet_blast/lkst_transcriptomes/${transcriptome} -d /home/biomatt/scratch/sterlet_blast/sterlet_db -o /home/biomatt/scratch/sterlet_blast/lkst_blast/${transcriptome}_out.tsv --very-sensitive
