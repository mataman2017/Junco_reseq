#!/bin/bash
#SBATCH --job-name=filter_fasta
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node (usually 1 for serial jobs)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=8G                   # Memory per node (adjust as needed)
#SBATCH --time=04:00:00              # Maximum runtime in HH:MM:SS
#SBATCH --chdir=/projects/VONHOLDT/jsala/A_infile/E_fasta/D_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY_2/A_batch_raxml_scripts
#SBATCH --output=OE/%x-N%N-J%j-%t.out      # Output file
#SBATCH --error=OE/%x-N%N-J%j-%t.err       # Error file

SCRIPT="/projects/VONHOLDT/jsala/A_infile/E_fasta/D_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY_2/A_batch_raxml_scripts/loci.filtering.py"

python $SCRIPT
