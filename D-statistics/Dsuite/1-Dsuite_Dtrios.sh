#!/bin/bash
#SBATCH --job-name=Dsuite
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=96:00:00
#SBATCH --chdir=/home/jg2334/download/Dsuite/Build
#SBATCH --output=/projects/VONHOLDT/jsala/B_PopGen/D_suite/C_FUL_outgroup/OE/%x-N%N-J%j-%t.out
#SBATCH --error=/projects/VONHOLDT/jsala/B_PopGen/D_suite/C_FUL_outgroup/OE/%x-N%N-J%j-%t.err

GZVCF="/projects/VONHOLDT/jsala/A_infile/B_VCFs/JUN_227_VAR.u.vcf.gz"
DIRECTORY="/projects/VONHOLDT/jsala/B_PopGen/D_suite/C_FUL_outgroup"
POP=${DIRECTORY}/pop.txt
TREE=${DIRECTORY}/tree.nwk
DTRIOS="/home/jg2334/download/Dsuite/utils/DtriosParallel"

$DTRIOS -n run_205_ful_out -t $TREE --cores 32  $POP $GZVCF
