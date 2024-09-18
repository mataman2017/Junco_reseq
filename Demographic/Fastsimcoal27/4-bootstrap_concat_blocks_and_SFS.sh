#!/bin/bash
#SBATCH --job-name=1_concat_blocks_SFS
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 06:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap/OE/%x-N%N-J%j-%t.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap/OE/%x-N%N-J%j-%t.out
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=javier.sala@mncn.csic.es #Dirección a la que se envía

# ACTIVATE THE EASYSFS ENV PRIOR TO TUN THIS SCRIPT

PREFIX=PALDORCAN_69_bial_miss1_maf05_hwe0001.LDpruned.neutral
INFILE=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/E_VCFs/A_paldorcan_69/${PREFIX}.vcf.gz
WD=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap
POP=/home/csic/bbe/jsg/B_JUN_reseq/E_demography/C_fsc_bootstraps/pop_62.txt

# Generate 50 files each with randomly concatenated blocks and compute the SFS for each:
for i in {1..50}
do
  cd $WD
  # Make a new folder for each bootstrapping iteration:
  mkdir bs$i
  cd bs$i

  # Add the header to our new bootstrapped vcf file
  cat ../header > $PREFIX.bs.$i.vcf
  # Randomly add 100 blocks
  for r in {1..100}
  do
    cat `shuf -n1 -e ../$PREFIX.sites.*` >> ${PREFIX}.bs.$i.vcf
  done
  # Compress the vcf file again
  bgzip ${PREFIX}.bs.$i.vcf

  # Make an SFS from the new bootstrapped file
  easySFS.py -i ${PREFIX}.bs.$i.vcf.gz -p $POP -a -f --proj 30,24,38,32

  # Copy the observed SFS file into this folder renaming it to match the .tpl prefix
  cp ./output/fastsimcoal2/${PREFIX}_MSFS.obs  ${PREFIX}.bs.${i}_MSFS.obs

  # Say that it is finished with iteration $i
  echo bs$i" ready"

  cd ..
done
