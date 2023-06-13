#!/bin/bash

# Extract the region of interest:
bcftools view -r ScoVZU6_963__HRSCAF___984:20500001-20550000 JUN_227_VAR.vcf.gz -o JUN_227_ScoVZU6_963__HRSCAF___984_20500001-20550000.vcf.gz -Oz

# Intersect it with the gene annotation
bedtools intersect -a Jhye_annotation.gff3 -b JUN_227_ScoVZU6_963__HRSCAF___984_20500001-20550000.vcf.gz > JUN_227_ScoVZU6_963__HRSCAF___984_20500001-20550000.gff3

# Extract gene names
grep -oP '(?<=ref-gene=).*?(?=;)' JUN_227_ScoVZU6_963__HRSCAF___984_20500001-20550000.gff3 | sort | uniq > genes_peak_chr_1.txt

# Check if the gene format of genes_peak_chr_1.txt is consistent with ensembl IDS 
# DONE
# NOW MOVE TO R
