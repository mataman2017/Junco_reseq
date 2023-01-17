
OUT=paldorcan_nohyb_8_biall_maf005_hwe00001_mac1_q30_dp5_50_miss075_thin100.recode

vcftools --gzvcf paldorcan_nohyb_109_mac1.recode.vcf.gz \
  --indv junco_bm42_dorsalis_flags \
  --indv junco_ric4_caniceps_color \
  --indv junco_bm29_caniceps_utah \
  --indv junco_gre5_dorsalis_white \
  --indv junco_sac9_dorsalis_lincoln \
  --indv junco_sac10_dorsalis_lincoln \
  --indv junco_bm76_palliatus_mexic \
  --indv junco_bm71_palliatus_mexic \
  --remove-indels \
  --max-missing 0.75 \
  --minQ 30 \
  --minDP 5 --maxDP 50 \
  --min-alleles 2 --max-alleles 2 \
  --maf 0.05 \
  --hwe 0.0001 \
  --thin 100 \
  --recode --recode-INFO-all --stdout  | \
  bgzip -c > ${OUT}.vcf.gz
  
vcf2phylip.py -i ${OUT}.vcf.gz \
  --output-prefix $OUT  
 

  
  
