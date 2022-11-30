#!/bin/bash

### Split the phased chromosome VCFs for each population. LAUNCH THIS SCRIPT USING PARALLEL (e.g. parallel $ bash parallel_split.cmd :::: CAN DOR PAL AIK MEA MON ORE # etc...)
PAR=$1
DIR=/mnt/DATA/B_JUN_reseq/A_chr_split/B_VCFs/D_JUN_236/B_VCFs/
OUT=/mnt/DATA/B_JUN_reseq/A_chr_split/D_CHR_POP_VCFs/B_236/${PAR}/        #NO CAMBIAR, DIR DONDE VAN LOS OUTFILES
CHR=$(find $DIR -maxdepth 1 -name "*.vcf.gz" )
LIST=/mnt/DATA/B_JUN_reseq/A_chr_split/A_cmd/A_list/names_${PAR}.txt 

for i in $CHR
do
	NAME=$(basename "$i")   #ESTO SE QUEDA SOLO CON NOMBRE DEL ARCHIVO, NO CON $PATH
	NAME2=${PAR}${NAME%_Sc*}
	mkdir -p ${OUT}
	cd ${OUT}
	bcftools view -S $LIST -O z -o ${NAME2}.vcf.gz ${i} &&
	bcftools index ${NAME2}.vcf.gz &&
	echo "File $NAME2 splitted"
done

time
 
