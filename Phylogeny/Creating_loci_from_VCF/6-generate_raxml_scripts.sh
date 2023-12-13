# Set the directory containing the split FASTA files
fasta_batch_dir="batch_raxml_output"

# Loop through the subdirectories
for batch_dir in /projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/A_batch_raxml_scripts/batch_raxml_split/fasta_batch_*; do
    batch_name=$(basename $batch_dir)
    OUTDIR="/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/D_output_trees"
    mkdir -p $OUTDIR
    # Create a separate script for each batch
    script_name="/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/A_batch_raxml_scripts/RAxML_scripts/run_raxml_$batch_name.sh"

    cat > $script_name <<EOL
#!/bin/bash
#SBATCH --job-name=raxml_$batch_name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --chdir=$OUTDIR
#SBATCH --output=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/A_batch_raxml_scripts/OE_trees/%x-N%N-J%j-%t.out
#SBATCH --error=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/A_batch_raxml_scripts/OE_trees/%x-N%N-J%j-%t.err

LIST=$batch_dir  # Use the correct list file
OUT=$OUTDIR

for i in \$(cat \$LIST); do
    LALA=\$(basename \$i)
    NAME="\${LALA%.min0.fasta}"
    SEED=\$(echo \$RANDOM)
    raxmlHPC-PTHREADS -f a -s \$i -m GTRGAMMA -n \$NAME -o VUL1 -# 100 -x \$SEED -p \$SEED -T 2
done
EOL

    chmod +x $script_name
done
