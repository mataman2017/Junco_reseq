#!/bin/bash

DIRGEN=/mnt/DATA/0_REF_genome/Dryad/fasta/B_mitochondrion/mitochondrion.fasta
DIR=/mnt/DATA/B_JUN_reseq/H_mitocondrial/A_DATA_sym
OUT=/mnt/DATA/B_JUN_reseq/H_mitocondrial/B_raw2bam/B_BAMs

for i in $(ls $DIR)
        do cd ${DIR}/${i}
        if [  $(ls -1 *fq.gz | wc -l) == 2  ]
        then
		echo "##Run de una sola lane##"
		LANE=$(ls -1 *1.fq.gz | sed -E 's/(^.+L[0-9])(_[1-2].+$)/\1/')
		echo "---------------------"
		echo "Run de secuenciacion: "$LANE
		FASTQ1=$(ls -1 *1.fq.gz)
		FASTQ2=$(ls -1 *2.fq.gz)
		echo "Archivos FASTQ de la muestra "$i": "$FASTQ1" y "$FASTQ2
		ID=$(zgrep -m 1 '@' ${LANE}_1.fq.gz | cut -d ":" -f 3)
		echo "ID de la muestra: "$ID
		optR=$(echo "@RG\tID:"${i}"_"$ID"\tSM:"$i"\tLB:"$i)
		echo "Opcion -R en bwa: "$optR
		echo "Filtrando reads..."
		time fastp -i $FASTQ1 -o ${OUT}/${LANE}_1-clean.fq.gz -I $FASTQ2 -O ${OUT}/${LANE}_2-clean.fq.gz \
			-g -w 10 -q 5 -u 50 -n 15 -l 150 \
			--min_trim_length 10 \
			--overlap_diff_limit 1 \
			--overlap_diff_percent_limit 10 \
			-j $LANE.json -h $LANE.html
		echo "Mapeando reads..."
		time bwa mem -t 10 -k 32 -M -R $optR $DIRGEN ${OUT}/${LANE}_1-clean.fq.gz ${OUT}/${LANE}_2-clean.fq.gz | \
			samtools fixmate -@10 -m - - | \
			samtools sort -@5 -T ${OUT}/tmp${i}${ID} - | \
			samtools markdup -@10 -r -T ${OUT}/tmp${i}${ID} - ${OUT}/$i.bam
		rm ${OUT}/${LANE}_1-clean.fq.gz ${OUT}/${LANE}_2-clean.fq.gz
		echo "-----TERMINADO-----"
	else
		echo "@@Run de dos lanes@@"
		for j in $(ls -1 *1.fq.gz | sort)
			do
			LANE=$(echo $j | sed -E 's/(^.+L[0-9])(_[1-2].+$)/\1/')
			echo "---------------------"
			echo "Run de secuenciacion: "$LANE
			FASTQ1=$(ls -1 ${LANE}_1.fq.gz)
			FASTQ2=$(ls -1 ${LANE}_2.fq.gz)
			echo "Archivos FASTQ de la muestra "$i": "$FASTQ1" y "$FASTQ2
			ID=$(zgrep -m 1 '@' ${LANE}_1.fq.gz | cut -d ":" -f 3)
			echo "ID de la muestra: "$ID
			optR=$(echo "@RG\tID:"${i}"_"$ID"\tSM:"$i"\tLB:"$i)
			echo "Opcion -R en bwa: "$optR
			echo "Filtrando reads..."
			time fastp -i $FASTQ1 -o ${OUT}/${LANE}_1-clean.fq.gz -I $FASTQ2 -O ${OUT}/${LANE}_2-clean.fq.gz \
				-g -w 10 -q 5 -u 50 -l 150 -n 15 \
				--min_trim_length 10 \
				--overlap_diff_limit 1 \
				--overlap_diff_percent_limit 10 \
				-j $LANE.json -h $LANE.html
			echo "Mapeando reads..."
			time bwa mem -t 10 -k 32 -M -R $optR $DIRGEN ${OUT}/${LANE}_1-clean.fq.gz ${OUT}/${LANE}_2-clean.fq.gz | \
				samtools fixmate -@10 -m - - | \
				samtools sort -@5 -T ${OUT}/tmp${i}${ID} - | \
				samtools markdup -@10 -r -T ${OUT}/tmp${i}$ID - ${OUT}/${i}${ID}.bam
			rm ${OUT}/${LANE}_1-clean.fq.gz ${OUT}/${LANE}_2-clean.fq.gz
			echo "Lane terminada"
		done
		echo "Juntando lanes de esta muestra..."
		BAMFILES=$(ls -1 ${OUT}/${i}???*.bam)
		samtools merge -@4 ${OUT}/${i}.bam $BAMFILES
		rm $BAMFILES
		echo "-----TERMINADO-----"
	fi
done
