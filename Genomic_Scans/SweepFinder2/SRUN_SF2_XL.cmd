#!/bin/bash
#SBATCH --job-name=SF2_SRUN_XL
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/0_log/SF2/236/XL/
#SBATCH -N 8
#SBATCH --ntasks-per-node 8
#SBATCH -n 60
#SBATCH -c 1
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/0_log/SF2/236/XL/%x-N%N-J%j-%t.err
#SBATCH --mail-type=all #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=javier.sala@mncn.csic.es #Dirección a la que se envía

## UN TOTAL DE 60 TRABAJOS

OUTERR=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/0_log/SF2/227/XL
IN=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/B_SF2_in/227/ALL_in/XL/
OUT=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/C_SF2_out/227/XL/
SRUN="srun --exclusive -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} -e $OUTERR/%x-N%N-J%j-%t-{1}.err -o /dev/null"

mkdir -p $OUTERR
mkdir -p $OUT

parallel --delay .2 --joblog parallel.log --colsep ' ' -j $SLURM_NTASKS $SRUN \
        SweepFinder2 -lu \
        ${IN}{1}.grid \
        ${IN}{1}.in \
        ${IN}{2}_SpectFile \
        ${OUT}{1}_2.out :::: /home/csic/bbe/jsg/B_JUN_reseq/B_Indexes/C_227/A_lists/SF2_XL.txt
