#!/bin/bash
#SBATCH -J JointVarInvar
#SBATCH -N 56
#SBATCH -n 56
#SBATCH -c 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -t 0-03:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/5-FiltVCFs/3-Joint
#SBATCH -o /home/csic/bbe/esl/2-Envios/1-Juncos/4-OE-filtvcfs/%x-%N-%j.out
#SBATCH -e /home/csic/bbe/esl/2-Envios/1-Juncos/4-OE-filtvcfs/%x-%N-%j.err

OUTERR="/home/csic/bbe/esl/2-Envios/1-Juncos/4-OE-filtvcfs/OE-Joint"
INPATH="/mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/5-FiltVCFs/2-Filtered/"
SRUN="srun --exclusive -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem-per-cpu=${SLURM_MEM_PER_CPU} -o $OUTERR/%x-%N-%j-%t-{}.out -e $OUTERR/%x-%N-%j-%t-{}.err"

parallel -j 56 $SRUN bcftools concat -a -O z -o {}.vcf.gz \
 $INPATH{}-var.vcf.gz \
 $INPATH{}-invar.recode.vcf.gz \
'&&' tabix -p vcf {}.vcf.gz  :::: /home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles/chr_30zSplit.txt
