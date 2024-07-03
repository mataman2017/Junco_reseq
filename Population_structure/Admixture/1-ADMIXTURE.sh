#/bin/bash 

module load ADMIXTURE/1.3.0
export _JAVA_OPTIONS='-Xmx100g'

plink2 --vcf $FILE.vcf.gz --make-bed --out $FILE --allow-extra-chr --const-fid

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

for K in 2 3 4 5 6 7 8 9 10 
do

    admixture --cv -B200 $FILE.bed $K | tee log${K}.out

done

grep -h CV log*.out > bestK.txt

