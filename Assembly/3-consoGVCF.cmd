#!/bin/bash

OUT=/mnt/DATA/B_JUN_reseq/H_mitocondrial/D_geno/B_VCFs
INFILE=/mnt/DATA/B_JUN_reseq/H_mitocondrial/C_bam2gvcf/B_gvcf
GENOME=/mnt/DATA/0_REF_genome/Dryad/fasta/B_mitochondrion/mitochondrion.fasta
READS=$(find $INFILE -name "*.g.vcf.gz")

### loop to inflate the variant option
OPTION=""
for i in $READS
do
	OPTION="${OPTION} -V ${i}"
done
echo $OPTION


gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
 --genomicsdb-workspace-path MITO_alljuncos240DB \
 -L /mnt/DATA/B_JUN_reseq/H_mitocondrial/D_geno/A_cmd/chr_mitoch.bed \
 $OPTION \
 --create-output-variant-md5 true \
 --max-num-intervals-to-import-in-parallel 10

