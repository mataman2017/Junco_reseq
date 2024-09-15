#!/bin/bash

# Select 100Kb region with ASIP inside
bcftools view -r ScoVZU6_1043__HRSCAF___1065:13200000-13300000 -Oz -o JUN_227_ASIP_joint.vcf.gz /mnt/DATA/B_JUN_reseq/A_chr_split/B_VCFs/J_Joint/ScoVZU6_1043__HRSCAF___1065.vcf.gz     

#Calculate depth
vcftools --gzvcf JUN_227_ASIP_joint.vcf.gz --depth --out vcf_depth.txt   

#paste pop
cat vcf_depth.txt | awk '{print $1}'| cut -b 1-3 > c
paste vcf_depth.txt c > vcf_depth_pop.txt

# keep 2 top depth per pop
awk 'NR>1 {print $4, $3, $1}' vcf_depth_pop.txt | sort -k1,1 -k2,2nr | awk '!seen[$1]++{count[$1]=0} count[$1]<2 {print; count[$1]++}' | sort -k1,1 -k2,2nr

# Convert VCF2fasta
vcf2phylip.py -i JUN_227_ASIP_joint_2ind.vcf.gz -p -f -m 1

# Convert fasta2HMM 
./Fasta2HMMFeature -i /mnt/DATA/B_JUN_reseq/I_phylogeny/G_local_phy/A_ASIP/JUN_227_ASIP_joint_2ind.min1.fasta -nosame -o /mnt/DATA/B_JUN_reseq/I_phylogeny/G_local_phy/A_ASIP/JUN_227_ASIP_joint_2ind.min1.hmm

# Run Saguaro
./Saguaro -f /mnt/DATA/B_JUN_reseq/I_phylogeny/G_local_phy/A_ASIP/JUN_227_ASIP_joint_2ind.min1.hmm -o /mnt/DATA/B_JUN_reseq/I_phylogeny/G_local_phy/A_ASIP/ASIP
