#!/bin/bash

# REMEMBER TO ACTIVATE PIXY ENV conda

# Define variables
files="/mnt/DATA/B_JUN_reseq/A_chr_split/B_VCFs/J_INVAR/Sco*gz" # modify to ScU*gz for sex chr
WorkDir="/mnt/DATA/B_JUN_reseq/B_PopIndexes/G_TajD/C_pop/A_auto"

for i in DOE DOW THU
do
    indv=$(cat pop_pop.txt | grep $i | awk '{print "--indv " $1}' | tr "\n" "\t")

    # Get list of chromosome files
    ls $files | parallel vcftools \
    --gzvcf {} \
    $indv \
    --TajimaD 50000 \
    --out {.}_lin_${i}_tajD.txt

done
