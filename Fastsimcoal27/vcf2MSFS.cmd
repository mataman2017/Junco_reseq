
### First run easySFS with the --preview option to check the suitability of the number of samples to use in each population. 

./easySFS.py \
	-i /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/A_BCN/C_paldorcan/paldorcan_mac1_hwe00001_bial_maf005_miss08_neutral.LDpruned.vcf.gz \
	-p /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/A_BCN/C_paldorcan/pop.file \
	-a -f --preview

### Then run easySFS
./easySFS.py \
	-i /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/A_BCN/C_paldorcan/paldorcan_mac1_hwe00001_bial_maf005_miss08_neutral.LDpruned.vcf.gz \
	-p /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/A_BCN/C_paldorcan/pop.file \
	-a -f --proj 34,24,36,40
