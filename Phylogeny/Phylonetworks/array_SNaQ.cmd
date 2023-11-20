#!/bin/bash
#SBATCH --output=/projects/VONHOLDT/jsala/C_Phylo/E_Phylonetworks/OE/%x-N%N-J%j-%t.out
#SBATCH --error=/projects/VONHOLDT/jsala/C_Phylo/E_Phylonetworks/OE/%x-N%N-J%j-%t.err
#SBATCH -J runsnaq
#SBATCH --array=0-2
#SBATCH -c 10
#SBATCH --mem=16G

## --array: to run multiple instances of this script,
##          one for each value in the array.
##          1 instance = 1 task
## -J job name
## -c number of cores (CPUs) per task

echo "slurm task ID = $SLURM_ARRAY_TASK_ID used as hmax"
echo "start of SNaQ parallel runs on $(hostname)"
# finally: launch the julia script, using Julia executable appropriate for slurm, with full paths:
~/bin/julia --history-file=no -- runSNaQ.jl $SLURM_ARRAY_TASK_ID 30 > net$SLURM_ARRAY_TASK_ID_30runs.screenlog 2>&1
echo "end of SNaQ run ..."
