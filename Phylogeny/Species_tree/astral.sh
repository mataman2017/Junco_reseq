#!/bin/bash
#SBATCH --job-name=ALL_astral_minisample
#SBATCH --chdir=/home/jg2334/download/ASTRAL
#SBATCH --output=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/E_Astral/%x-N%N-J%j-%t.out
#SBATCH --error=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/E_Astral/%x-N%N-J%j-%t.err       # Error file
#SBATCH --ntasks=1
#SBATCH --mem=32G                   # Memory per node (adjust as needed)
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00

INDIR=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/E_Astral
TREES=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/D_output_trees

for i in $(ls $TREES/*best*)
do
cat $i >> $INDIR/Trees_concat.txt
done

java -D"java.library.path=lib/" -jar astral.5.15.5.jar -i $INDIR/Trees_concat.txt -o $INDIR/FINAL_tree -a $INDIR/pop.txt -s $RANDOM
