#!/bin/bash
for j in CAN PAL MON
do
DIR=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/A_POP_VCFs/227/$j
OUT=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/B_SF2_in/227/$j
CHR=$(find $DIR -name "*vcf.gz" )
	for i in $CHR
	do
		NAME=$(basename "$i")   #ESTO SE QUEDA SOLO CON NOMBRE DEL ARCHIVO, NO CON $PATH
		NAME2=${NAME%.recode.vcf.gz}
		mkdir -p $OUT
		vcftools --gzvcf $i --counts2 --out $OUT/$NAME2
		tail -n+2 $OUT/$NAME2.frq.count | awk -v OFS="\t" '{print $2,$6,$4,"1"}' > $OUT/$NAME2.in
		echo -e 'position\tx\tn\tfolded' | cat - $OUT/$NAME2.in > temp_$NAME2 && mv temp_$NAME2 $OUT/$NAME2.in
		echo "done for ${NAME2} ${j}"
	done
done

for j in CAN PAL MON
do
DIR=/mnt/lustre/scratch/nlsas/home/csic/bbe/jsg/B_JUN_reseq/B_SF2_in/227/$j
CHR=$(find $DIR -name "*.in" )
	for i in $CHR
	do
		NAME=$(basename "$i")   #ESTO SE QUEDA SOLO CON NOMBRE DEL ARCHIVO, NO CON $PATH
		NAME2=${NAME%.in}
		awk 'FNR>1 {print $1}' $i > $DIR/$NAME2.grid
		echo "done for ${NAME2} ${j}"
	done
done
