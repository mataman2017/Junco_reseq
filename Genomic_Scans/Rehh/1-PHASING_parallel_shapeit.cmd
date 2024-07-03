REC=/mnt/DATA/B_JUN_reseq/N_genetic_map/E_GOOD_rmap/D_rmap_shapeit_format/rmap_
DIR=/mnt/DATA/B_JUN_reseq/A_chr_split/B_VCFs/E_JUN_227/chromosome
OUT=/mnt/DATA/B_JUN_reseq/A_chr_split/B_VCFs/E_JUN_227/C_polarize_correct
CHR=$(find $DIR -name "*.vcf.gz" )

SHAPEIT=/home/milab/lib/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit

for i in $CHR
do
	NAME=$(basename "$i")   #ESTO SE QUEDA SOLO CON NOMBRE DEL ARCHIVO, NO CON $PATH
	NAME2=${NAME%.vcf.gz}_phased
	NAME3=${NAME%.vcf.gz}
	mkdir -p $OUT
	$SHAPEIT --input-vcf $i \
		-O $OUT/$NAME2 \
		--window 0.5 \
		-T 10 \
		--force \
		--input-map ${REC}$NAME3.txt
	echo "Chromosome $NAME phased"
	$SHAPEIT -convert \
        	--input-haps $OUT/$NAME2 \
        	--output-vcf $OUT/$NAME2.vcf &&

	echo "Chromosome $NAME2 converted to .vcf"

	bgzip $OUT/$NAME2.vcf &&
	bcftools index $OUT/$NAME2.vcf.gz

	echo "Chromosome $NAME2 bgzipped"
done
