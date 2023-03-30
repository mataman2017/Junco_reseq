#!/bin/bash

## Phasing chromosome VCF's with shapeit, parallel run. 


parallel \
shapeit --input-vcf /DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/{}.recode.vcf.gz \
	-O /DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/{}_phased \
	--window 0.5 -T 20 --force "&&" \
shapeit -convert --input-haps /DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/{}_phased \
	--output-vcf /DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/{}_phased.vcf "&&" \
bgzip /DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/{}_phased.vcf "&&" \
bcftools index /DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/{}_phased.vcf.gz \
:::: /mnt/6TB/1-Juncos/4-PopIndexes/2-Rehh/A_cmd/236/A_list/chr_list.txt

