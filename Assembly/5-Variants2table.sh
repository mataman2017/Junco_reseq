#!/bin/bash
#SBATCH -J var2table
#SBATCH -a 1-31
#SBATCH --mem-per-cpu 4GB
#SBATCH -t 0-06:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/4-ConsoGVCF/
#SBATCH -o /home/csic/bbe/esl/2-Envios/1-Muestras_Maria/4-OE-filtvcfs/OE-var2table/%x-%a-%N-%j.out
#SBATCH -e /home/csic/bbe/esl/2-Envios/1-Muestras_Maria/4-OE-filtvcfs/OE-var2table/%x-%a-%N-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enrique.saez@mncn.csic.es

module load gatk

CHR=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/csic/bbe/esl/0-Soft_y_Utils/2-Utiles/chr_30z.txt)
OUT="/mnt/lustre/scratch/nlsas/home/csic/bbe/esl/1-Datos/5-FiltVCFs/1-Before"

gatk --java-options "-Xmx3g" VariantsToTable \
 -V $CHR.vcf.gz \
 -F POS -F TYPE -F MULTI-ALLELIC -F QUAL -F DP -F NSAMPLES -F NCALLED \
 -F VAR -F HOM-REF -F HET -F HOM-VAR \
 -F AF -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -F InbreedingCoeff \
 -O $OUT/$CHR.table
