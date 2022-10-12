#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

cd /home/biomatt

module load StdEnv/2020 gcc python 
module load augustus hmmer blast+ gcc/9.3.0 openmpi/4.0.3 metaeuk/4-a0f584d prodigal r

virtualenv /home/$USER/busco_env
source /home/biomatt/scratch/overall_transcriptome/busco_env/bin/activate

cd /home/biomatt/scratch/overall_transcriptome

busco --offline --in /home/biomatt/scratch/overall_transcriptome/trinity_overall.Trinity.fasta --cpu $SLURM_CPUS_PER_TASK --out BUSCO --lineage_dataset /home/biomatt/projects/def-kmj477/biomatt/redside_dace/BUSCO/actinopterygii_odb10 --mode transcriptome

busco --offline --in /home/biomatt/scratch/overall_transcriptome/trinity_overall.Trinity.fasta --cpu 2 --out BUSCO --lineage_dataset /home/biomatt/scratch/overall_transcriptome/actinopterygii_odb10 --mode transcriptome
