#!/bin/bash

# Create a vcf file for the outgroup individual only
#vcftools --gzvcf test.vcf.gz --indv PHA1 --recode --recode-INFO-all --out PHA1
bcftools view -s PHA2 --threads 18 -Oz -o PHA2.vcf.gz PhaPDC_68_bial_maf.05.vcf.gz

# Create a table of SNP positions and alleles using bcftools
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' PHA2.vcf.gz  > file.tab

# Create a table file with ancestral allele information. 
# Here we select the reference allele (REF) or the variant allele (ALT) and put it in the 5th column. 
# If the site is not covered by our outgroup reads, i.e. missing '.', then we record it as missing.
awk '{OFS="\t";if($5=="0/0"){print $1,$2,$3,$4,$3} \
        if($5=="0/1"){print $1,$2,$3,$4,$4} \
        if($5=="./."){print $1,$2,$3,$4,$5}}' file.tab > file_aa.tab

# Compress and index the table file and the original vcf file
bgzip file_aa.tab
tabix -s1 -b2 -e2 file_aa.tab.gz
#bgzip file_BI_SNPS.recode.vcf

# Create an INFO file line for the new vcf file
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt

# Using bcftools to annotate the vcf file with the ancestral allele information
bcftools annotate -a file_aa.tab.gz \
 -c CHROM,POS,REF,ALT,INFO/AA -h hdr.txt -Oz \
 -o Polarized_PhaPDC_68_bial_maf.05.vcf.gz PhaPDC_68_bial_maf.05.vcf.gz

#  Check that it has worked. There should be an info field AA
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' Polarized_PhaPDC_68_bial_maf.05.vcf.gz | less

# Count how many sites have ancestral allele information

bcftools view -e 'INFO/AA=="."' Polarized_PhaPDC_68_bial_maf.05.vcf.gz -H | wc -l
