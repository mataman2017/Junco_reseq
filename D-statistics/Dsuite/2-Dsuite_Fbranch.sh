#!/bin/bash
#SBATCH --job-name=Dsuite_fbranch_C
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --chdir=/home/jg2334/download/Dsuite
#SBATCH --output=/projects/VONHOLDT/jsala/B_PopGen/D_suite/C_FUL_outgroup/OE/%x-N%N-J%j-%t.out
#SBATCH --error=/projects/VONHOLDT/jsala/B_PopGen/D_suite/C_FUL_outgroup/OE/%x-N%N-J%j-%t.err

GZVCF="/projects/VONHOLDT/jsala/A_infile/B_VCFs/JUN_227_VAR.u.vcf.gz"
DIRECTORY="/projects/VONHOLDT/jsala/B_PopGen/D_suite/C_FUL_outgroup"
POP=${DIRECTORY}/pop.txt
TREE=${DIRECTORY}/tree.nwk
DTRIOS="/home/jg2334/download/Dsuite/utils/DtriosParallel"
F4RATIO=${DIRECTORY}/DTparallel_pop_run_205_ful_out_combined_tree.txt

/home/jg2334/download/Dsuite/Build/Dsuite Fbranch $TREE $F4RATIO -Z
