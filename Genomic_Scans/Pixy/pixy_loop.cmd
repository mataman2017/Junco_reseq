#!/bin/bash

# Define variables
pop=/mnt/DATA/B_JUN_reseq_variants/B_pixy/A_cmd/pops.txt
files="/mnt/DATA/B_JUN_reseq_variants/A_chr_split/B_VCFs/chromosomes/*"
log=log_file.txt

# create log file or overrite if already present
printf "Log File - " > $log

# append date to log file
date >> $log
# activate  env where pixy is located

for i in $files

do

echo "runnning pixy for $i"

filename=$(basename "$i" .vcf.gz) 
tabix $i
pixy --stats  fst \
--vcf $i \
--populations $pop \
--window_size 50000 \
--n_cores 20 \
--output_prefix $filename \
--bypass_invariant_check 'yes'

done >> $log
