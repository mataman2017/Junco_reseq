ldPruning.sh data.vcf.gz

plink2 --vcf data_ldpruned.vcf.gz --make-bed --out data --allow-extra-chr

plink2 --bfile data --out data --allow-extra-chr --pca
