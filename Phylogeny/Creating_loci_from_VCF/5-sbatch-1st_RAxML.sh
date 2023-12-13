#!/bin/bash

for i in {01..89}; 
do 
sbatch  RAxML_scripts/run_raxml_fasta_batch_$i.sh
done
