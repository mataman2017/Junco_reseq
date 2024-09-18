#!/bin/bash
#SBATCH --job-name=split_boot_fsc
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap/OE/%x-N%N-J%j-%t.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap/OE/%x-N%N-J%j-%t.out
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=javier.sala@mncn.csic.es #Dirección a la que se envía

PREFIX=PALDORCAN_69_bial_miss1_maf05_hwe0001.LDpruned.neutral
INFILE=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/E_VCFs/A_paldorcan_69/${PREFIX}.vcf.gz
WD=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/C_fsc_bootstrap

# Get all lines with genomic data
#zgrep -v "^#" $INFILE > $WD/$PREFIX.allSites

# Get the header
#zgrep "^#" $INFILE > $WD/header

# get 100 files with 4338 sites each (number 101 removed due to only 90 sites)
split -l 34247 $WD/$PREFIX.allSites $WD/$PREFIX.sites.
