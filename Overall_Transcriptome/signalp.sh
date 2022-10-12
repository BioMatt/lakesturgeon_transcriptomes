#!/bin/bash
#SBATCH --account=def-dogfish
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=50000M

module load StdEnv/2020 signalp/4.1f

perl $EBROOTSIGNALP/signalp -f short -n /scratch/ibouyouc/dogfish/annotate/signalp.out /scratch/ibouyouc/dogfish/annotate/longest_orfs.pep