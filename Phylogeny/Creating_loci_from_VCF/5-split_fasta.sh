#!/bin/bash

# Define your input and output directories
input_dir=/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/C_filt_fasta_files
output_dir="/projects/VONHOLDT/jsala/A_infile/E_fasta/C_JUN_225_2kb_20kb_masked-miss_auto_neutral_hetHY/A_batch_raxml_scripts/batch_raxml_split"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Split the files into subdirectories with a maximum of 50 files per directory
find $input_dir -maxdepth 1 -type f -name "*.fasta" | split -l 50 -d - $output_dir/fasta_batch_
