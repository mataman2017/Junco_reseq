#!/bin/bash

IN=/mnt/DATA/B_JUN_reseq/H_mitocondrial/B_raw2bam/B_BAMs
OUT=/mnt/DATA/B_JUN_reseq/H_mitocondrial/C_bam2gvcf/B_gvcf
REF=/mnt/DATA/0_REF_genome/Dryad/fasta/B_mitochondrion/mitochondrion.fasta

for i in $(cat list8.txt)
do
	samtools index -@4 $IN/$i.bam
	gatk HaplotypeCaller -R $REF \
	-I $IN/$i.bam \
	--native-pair-hmm-threads 4 \
	-O $OUT/$i.g.vcf.gz \
	-ERC GVCF
done

for j in $(cat list8.txt)
do
	rm $IN/$j.bam $IN/$j.bam.bai
done
