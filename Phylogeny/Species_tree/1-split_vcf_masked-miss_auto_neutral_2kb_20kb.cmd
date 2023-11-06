#!/bin/bash
#SBATCH --job-name=split_vcf_4kb_100kb_masked-miss
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node (usually 1 for serial jobs)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=40G                   # Memory per node (adjust as needed)
#SBATCH --time=12:00:00              # Maximum runtime in HH:MM:SS
#SBATCH --chdir=/_/A_infile/B_VCFs/C_split_2kb_100kb/D_225_2kb_20kb_masked-miss_auto_neutral_hetHY
#SBATCH --output=./OE/%x-N%N-J%j-%t.out      # Output file
#SBATCH --error=./OE/%x-N%N-J%j-%t.err       # Error file


# Input VCF file
INDIR="/_/A_infile/B_VCFs/0_Joint/B_225_nomaf_mask-miss_auto_neutral"

for x in $(ls $INDIR/*gz)
do

bcft='/_/BIN/bcftools-1.9/bcftools'

# Output directory
output_dir="/_/A_infile/B_VCFs/C_split_2kb_100kb/D_225_2kb_20kb_masked-miss_auto_neutral"

# Chromosome name
LALA=$(basename $x)
chromosome="${LALA%.vcf.gz}"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Calculate the total number of variants
num_variants=$($bcft view -h $x  | grep -m 1 "##contig=<ID=$chromosome" | awk -F'length=' '{print $2}'  | awk -F',' '{print $1}' | sed 's/>//' )

# Window size (2kb)
window_size=2000

# Spacing between windows (100kb)
spacing=100000

# Initialize a counter for naming output files
counter=1

# Initialize the start position
start_pos=1

while ((start_pos <= num_variants)); do
    end_pos=$((start_pos + window_size - 1))

    # Check if the end position exceeds the total number of variants
    if ((end_pos > num_variants)); then
        end_pos=$num_variants
    fi

    output_vcf="$output_dir/${chromosome}_${counter}.vcf.gz"

    # Use bcftools to extract the region and save it as a new VCF
    $bcft view -r "${chromosome}:${start_pos}-${end_pos}" -O z -o "$output_vcf" $x

    # Use bcftools to count the number of variants
    variants=$(bcftools view -H -v snps "$output_vcf" | wc -l)

    # Check if the number of variants is zero
    if [ "variants" -eq 0 ]; then
        # Delete the output VCF file
        rm -f "$output_vcf"
    fi

    # Increment the counter
    ((counter++))

    # Update the start position for the next window
    start_pos=$((end_pos + 1 + spacing))
done

rm -f $output_dir/*csi
done
