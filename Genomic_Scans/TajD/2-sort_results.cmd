#!/bin/bash

# REMEMBER TO ACTIVATE PIXY ENV conda

# Define variables
THU_TajD="/mnt/DATA/B_JUN_reseq/B_PopIndexes/G_TajD/C_pop/C_All/unsorted_DOE_DOW_THU/*_THU_*"
DOE_TajD="/mnt/DATA/B_JUN_reseq/B_PopIndexes/G_TajD/C_pop/C_All/unsorted_DOE_DOW_THU/*_DOE_*"
DOW_TajD="/mnt/DATA/B_JUN_reseq/B_PopIndexes/G_TajD/C_pop/C_All/unsorted_DOE_DOW_THU/*_DOW_*"

for i in $THU_TajD
do
	filename=$(basename "$i" .vcf_lin_THU_tajD.txt.Tajima.D)
	awk 'NR==1{$0=$0"\tpop"} NR>1{$0=$0"\tTHU"} 1' "$i" > "FRAGMENTS/${filename}_THU_TajD.txt"
done

for i in $DOE_TajD
do
        filename=$(basename "$i" .vcf_lin_DOE_tajD.txt.Tajima.D)
        awk 'NR==1{$0=$0"\tpop"} NR>1{$0=$0"\tDOE"} 1' "$i" > "FRAGMENTS/${filename}_DOE_TajD.txt"
done

for i in $DOW_TajD
do
        filename=$(basename "$i" .vcf_lin_DOW_tajD.txt.Tajima.D)
        awk 'NR==1{$0=$0"\tpop"} NR>1{$0=$0"\tDOW"} 1' "$i" > "FRAGMENTS/${filename}_DOW_TajD.txt"
done

head -n 1 FRAGMENTS/ScoVZU6_1043__HRSCAF___1065_DOW_TajD.txt > remaining_POPs_50kb_TajD.txt

for i in $(ls FRAGMENTS/*_TajD.txt)
do
tail -n +2 $i >> linages_50kb_TajD.txt
done

awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="window_pos_1"} NR>1{$(NF+1)=$2+1} 1' linages_50kb_TajD.txt | sed 's/CHROM/chromosome/g'> tmp

mv tmp linages_50kb_TajD.txt
    
