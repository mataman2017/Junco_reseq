#!/bin/bash
#SBATCH -J Bayescan1
#SBATCH -N 9
#SBATCH -n 9
#SBATCH -c 12
#SBATCH --mem=8GB
#SBATCH -t 4-00:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/7-Selection/1-Bayescan
#SBATCH -o /home/csic/bbe/esl/2-Envios/1-Juncos/6-OE-Selection/%x-%a-%N-%j.out
#SBATCH -e /home/csic/bbe/esl/2-Envios/1-Juncos/6-OE-Selection/%x-%a-%N-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enrique.saez@mncn.csic.es

module load gcccore/system bayescan/2.1

#--paths--
INPATH="/mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/7-Selection/0-VCF2BSC"
OUTERR="/home/csic/bbe/esl/2-Envios/1-Juncos/6-OE-Selection/OE-Bayescan"
UTILS="/home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles"

#POP=$(sed -n "$SLURM_ARRAY_TASK_ID"p $UTILS/chr_30zSplit.txt)
SRUN="srun --exclusive -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem=${SLURM_MEM_PER_NODE} -o $OUTERR/%x-%a-%N-%j-%t-{}.out -e $OUTERR/%x-%a-%N-%j-%t-{}.err"

parallel -j 9 $SRUN bayescan_2.1 \
 -snp $INPATH/{}.bsc \
 -threads 12 \
 -od Chr_under10M \
 -out_freq :::: $UTILS/chr_30zSplit1.txt
