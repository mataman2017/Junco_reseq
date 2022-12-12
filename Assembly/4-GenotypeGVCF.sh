#!/bin/bash

OUT=/mnt/DATA/B_JUN_reseq/H_mitocondrial/D_geno/B_VCFs
INFILE=/mnt/DATA/B_JUN_reseq/H_mitocondrial/C_bam2gvcf/B_gvcf
GENOME=/mnt/DATA/0_REF_genome/Dryad/fasta/B_mitochondrion/mitochondrion.fasta
READS=$(find $INFILE -name "*.g.vcf.gz")

gatk GenotypeGVCFs \
 -R $GENOME \
 -V gendb://MITO_alljuncos240DB \
 -O $OUT/mitochondrial_JUN_240.vcf.gz \
 --genomicsdb-shared-posixfs-optimizations TRUE \
 -all-sites true
