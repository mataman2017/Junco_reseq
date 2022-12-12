#!/bin/bash

#SBATCH -J consoGVCF
##SBATCH --test-only
#SBATCH -N 1
##SBATCH -n 4
#SBATCH -c 4
##SBATCH --mem-per-cpu=5GB
#SBATCH --mem 230GB
#SBATCH -t 6-00:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/4-ConsoGVCF
#SBATCH -o /home/csic/bbe/esl/2-Envios/1-Muestras_Maria/3-OE-gvcf2conso/%x-%N-%j.out
#SBATCH -e /home/csic/bbe/esl/2-Envios/1-Muestras_Maria/3-OE-gvcf2conso/%x-%N-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enrique.saez@mncn.csic.es


module load cesga/2020 gatk/4.2.0.0

time gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx200g -Xms200g" GenomicsDBImport \
 --genomicsdb-workspace-path alljuncos240DB \
 -L /home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles/chr_30zSplit.bed \
 --sample-name-map /home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles/alljuncos240.sample_map \
 --tmp-dir $LUSTRE_SCRATCH \
 --genomicsdb-shared-posixfs-optimizations true \
 --create-output-variant-md5 true \
 --max-num-intervals-to-import-in-parallel 10
