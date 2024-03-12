#!/bin/bash
#SBATCH -J filtervcf
#SBATCH -n 29
#SBATCH --mem-per-cpu 8GB
#SBATCH -t 0-06:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/5-FiltVCFs/2-Filtered
#SBATCH -o /home/csic/bbe/esl/2-Envios/1-Muestras_Maria/4-OE-filtvcfs/%x-%N-%j.out
#SBATCH -e /home/csic/bbe/esl/2-Envios/1-Muestras_Maria/4-OE-filtvcfs/%x-%N-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enrique.saez@mncn.csic.es

OUTERR="/home/csic/bbe/esl/2-Envios/1-Muestras_Maria/4-OE-filtvcfs/OE-INVAR"
SRUN="srun --exclusive -N 1 -n 1 --mem-per-cpu=${SLURM_MEM_PER_CPU} -o $OUTERR/%x-%N-%j-%t-{}.out -e $OUTERR/%x-%N-%j-%t-{}.err"
INPATH="/mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/4-ConsoGVCF"

cat /home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles/chr_30zXX.txt | \
parallel --delay .2 -j $SLURM_NTASKS $SRUN time vcftools \
 --gzvcf $INPATH/{}.vcf.gz \
 --max-maf 0 \
 --remove-filtered-all \
 --recode-INFO-all \
 --recode \
 --out {}-invar \
'&&' bgzip {}-invar.recode.vcf '&&' tabix -p vcf {}-invar.recode.vcf.gz
