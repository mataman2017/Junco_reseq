#!/bin/bash
#SBATCH --job-name=vcf2phy_mask-miss_auto_neutral_hetHY_2kb_20kb
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node (usually 1 for serial jobs)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=16G                   # Memory per node (adjust as needed)
#SBATCH --time=4:00:00              # Maximum runtime in HH:MM:SS
#SBATCH --chdir=/_/A_infile/A_cmd/F_vcf2fasta/B_new
#SBATCH --output=./OE/%x-N%N-J%j-%t.out      # Output file
#SBATCH --error=./OE/%x-N%N-J%j-%t.err       # Error file


# Input VCF file
INDIR="/_/A_infile/B_VCFs/C_split_2kb_100kb/D_225_2kb_20kb_masked-miss_auto_neutral"
OUT="/_/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/B_fasta_files"
CMD="~/bin/vcf2phyllip.py"

for x in $(ls $INDIR/*gz)
do

# Chromosome name
LALA=$(basename $x)
FRAGM="${LALA%.vcf.gz}"

$CMD -i $x --output-folder $OUT \
        --output-prefix $FRAGM \
        -o VUL1 \
        -f -p -m 0
done
