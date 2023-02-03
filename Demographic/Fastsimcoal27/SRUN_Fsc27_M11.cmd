#!/bin/bash
#SBATCH -J fsc_M11
#SBATCH -N 13
#SBATCH --ntasks-per-node 8
#SBATCH -n 100
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G
#SBATCH -t 6:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/B_fastsimcoal/B_PDC_62/M11/
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/B_fastsimcoal/B_PDC_62/M11/OE/%x-%a-%N-%j.out
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/B_fastsimcoal/B_PDC_62/M11/OE/%x-%a-%N-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javier.sala@mncn.csic.es

# Script to create directories for each of the fsc runs and copy the infiles required in each directory.
# Then, execute parallel runs of fastsimcoal in different directories through srun. 

module load fastsimcoal2/fsc2709

#--paths--
MODEL="M11"
INPATH="/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/F_demography/B_fastsimcoal/B_PDC_62/${MODEL}"
OUTERR="$INPATH/OE"
SRUN="srun --exclusive -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem-per-cpu=${SLURM_MEM_PER_CPU} -D $INPATH/run{}/ -o $OUTERR/%x-%a-%N-%j-%t-{}.out -e $OUTERR/%x-%a-%N-%j-%t-{}.err"

mkdir -p $INPATH/run{}/

for i in {1..100}
do
  scp  $INPATH/${MODEl}.tpl $INPATH/${MODEl}.est $INPATH/${MODEl}_MSFS.obs $INPATH/run${i}/
done

parallel --delay .2 --colsep ' ' -j $SLURM_NTASKS $SRUN \
        fsc2709 -t $MODEL.tpl -e $MODEL.est -n 150000 -M -u -m -L 50 -c 1 -B 1 --seed $RANDOM \
	::: {1..100}
