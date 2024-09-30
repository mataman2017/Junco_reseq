#!/bin/bash
#SBATCH -J VCF2BSC
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=24GB
#SBATCH -t 0-20:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/7-Selection/0-VCF2BSC/PRUEBAS
#SBATCH -o /home/csic/bbe/esl/2-Envios/1-Juncos/6-OE-Selection/OE-VCF2BSC/%x-%a-%N-%j-%t-ScoVZU6_1043__HRSCAF___1065-7500001-14944543.out
#SBATCH -e /home/csic/bbe/esl/2-Envios/1-Juncos/6-OE-Selection/OE-VCF2BSC/%x-%a-%N-%j-%t-ScoVZU6_1043__HRSCAF___1065-7500001-14944543.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enrique.saez@mncn.csic.es

#--paths--
#INPATH="/mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/5-FiltVCFs/2-Filtered"
#OUTERR="/home/csic/bbe/esl/2-Envios/1-Juncos/6-OE-Selection/OE-VCF2BSC"
UTILS="/home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles"

#POP=$(sed -n "$SLURM_ARRAY_TASK_ID"p $UTILS/chr_30zSplit.txt)
#SRUN="srun --exclusive -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem-per-cpu=${SLURM_MEM_PER_CPU} -o $OUTERR/%x-%a-%N-%j-%t-{}.out -e $OUTERR/%x-%a-%N-%j-%t-{}.err"

perl /home/csic/bbe/esl/0-Soft_y_Utils/1-Soft/vcf2bayescan/vcf2bayescanESL.pl \
 -p $UTILS/pops.txt \
 -v ScoVZU6_1043__HRSCAF___1065-7500001-14944543-var.vcf \
 -o ScoVZU6_1043__HRSCAF___1065_2.bsc
