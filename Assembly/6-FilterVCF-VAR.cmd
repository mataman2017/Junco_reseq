#!/bin/bash
#SBATCH -J paldorcan_filtervcf-VAR
#SBATCH -n 43
#SBATCH --mem-per-cpu 16GB
#SBATCH -t 06:00:00
#SBATCH -D /mnt/lustre/scratch/nlsas/home/csic/bbe/bmi/A_JAVI/A_filtered_VCFs/B_paldorcan/
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/bbe/bmi/A_JAVI/B_OUTERR/B_filtering_paldorcan/%x-%N-%j.out
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/bbe/bmi/A_JAVI/B_OUTERR/B_filtering_paldorcan/%x-%N-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javier.sala@mncn.csic.es

module load gatk/4.2.0.0

OUTERR="/mnt/lustre/scratch/nlsas/home/csic/bbe/bmi/A_JAVI/B_OUTERR/B_filtering_paldorcan"
SRUN="srun --exclusive -N 1 -n 1 --mem-per-cpu=${SLURM_MEM_PER_CPU} -o $OUTERR/%x-%N-%j-%t-{}.out -e $OUTERR/%x-%N-%j-%t-{}.err"
INPATH="/mnt/lustre/scratch/nlsas/home/csic/bbe/bmi/4-ConsoGVCF"
OUT="/mnt/lustre/scratch/nlsas/home/csic/bbe/bmi/A_JAVI/A_filtered_VCFs/B_paldorcan/"

cat /home/csic/bbe/bmi/0-Soft_y_Utils/2-Utiles/chr_30z.txt | \
parallel --delay .2 -j $SLURM_NTASKS $SRUN time gatk --java-options \"-Xmx10g\" SelectVariants \
 -V $INPATH/{}.vcf.gz \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 --max-nocall-fraction 0.2 \
 --exclude-sample-name /home/csic/bbe/bmi/0-Soft_y_Utils/2-Utiles/exclude_list_paldorcan.args \
 --exclude-filtered \
 --output {}-prefiltered1.vcf.gz \
 -OVI \
'&&' \
time gatk --java-options \"-Xmx10g\" VariantFiltration \
 -V {}-prefiltered1.vcf.gz \
 --output {}-prefiltered2.vcf.gz \
 --filter-name \"AF0.05\" -filter \"AF'<'0.05\" \
 --filter-name \"QD2\" -filter \"QD'<'2.0\" \
 --filter-name \"FS60\" -filter \"FS'>'60.0\" \
 --filter-name \"SOR3\" -filter \"SOR'>'3.0\" \
 --filter-name \"MQ40\" -filter \"MQ'<'40.0\" \
 --filter-name \"MQRS-12.5\" -filter \"MQRankSum'<'-12.5\" \
 --filter-name \"RPRS-8\" -filter \"ReadPosRankSum'<'-8.0\" \
 -OVI \
'&&' \
rm {}-prefiltered1.vcf.gz* \
'&&' \
time gatk --java-options \"-Xmx10g\" SelectVariants \
 -V {}-prefiltered2.vcf.gz \
 --exclude-filtered \
 --output {}_var.vcf.gz \
 -OVI \
'&&' \
rm {}-prefiltered2.vcf.gz*
