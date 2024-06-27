# 1. filter SNPs

vcftools --gzvcf JUN.vcf.gz --keep group.txt --max-missing 1 --mac 1 --exclude-bed exclude.bed --recode --stdout | gzip -c > JUN.group.vcf.gz

# 2. construction of the site frequency spectrum
easySFS.py -i JUN.group.vcf.gz -p JUN.group.txt --proj={2*sample size} -a -o output_group --prefix group

# 3. run stairway_plot
java -cp stairway_plot_es Stairbuilder group.blueprint
sh group.blueprint.sh
