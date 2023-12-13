#!/bin/bash
#SBATCH --job-name=sample_rm_vcf2phy_mask-miss_auto_neutral_hetHY_2kb_20kb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=20:00:00
#SBATCH --chdir=/projects/VONHOLDT/jsala/A_infile/A_cmd/F_vcf2fasta/B_new
#SBATCH --output=./OE/%x-N%N-J%j-%t.out
#SBATCH --error=./OE/%x-N%N-J%j-%t.err

# Input VCF file
INDIR="/projects/VONHOLDT/jsala/A_infile/B_VCFs/C_split_2kb_100kb/E_225_2kb_20kb_masked-miss_auto_neutral"
OUT="/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/B_fasta_files"
CMD="/home/jg2334/bin/vcf2phyllip.py"

# Use a while loop to iterate over files
find "$INDIR" -name "*_3544_*.vcf.gz" | while read -r x; do
    # Chromosome name
    LALA=$(basename "$x")
    FRAGM="${LALA%.vcf.gz}"

    $CMD -i "$x" --output-folder "$OUT" \
            --output-prefix "$FRAGM" \
            -o VUL1 \
            -f -p -m 0

    python /projects/VONHOLDT/jsala/A_infile/A_cmd/F_vcf2fasta/B_new/remove_filt.py "$OUT/$FRAGM.min0.fasta" "$OUT/$FRAGM.fasta"

    rm "$OUT/$FRAGM.min0.fasta"
done
