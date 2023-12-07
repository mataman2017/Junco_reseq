


# Extracts population information from a VCF file and creates a clusters file.
bcft query -l $file.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[1]}' > birds.clust

# Converts PED and MAP files to PLINK binary format.
plink --file $file --make-bed --out $file --allow-no-sex --allow-extra-chr 0

# Calculates allele frequencies and missingness within specified populations.
plink --bfile $file --freq --missing --within dogs.clust --out $file --allow-no-sex --allow-extra-chr 0

# Compresses a file containing stratified allele frequencies.
gzip $file".frq.strat"

# Converts PLINK frequency file to TreeMix format.
./plink2treemix.py $file".frq.strat.gz" $file".treemix.frq.gz"

# Decompresses the TreeMix frequency file.
gunzip $file".treemix.frq.gz"

# Decompresses the stratified allele frequencies file.
gunzip $file".frq.strat.gz"

# Extracts scaffold position information from a PLINK map file.
awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' $file".map" > $file".positions"

# Combines scaffold positions and TreeMix allele frequencies.
paste $file".positions" $file".treemix.frq" > $file".frequencies"

# Calculates frequency values and creates a new frequency file.
awk '{printf $0 for(i = 4; i <= NF; i++){split($i,values,",") if((values[1]+values[2])>0) freq=values[1]/(values[1]+values[2]) else freq=0 printf freq"\t" } printf "\n"}' $file".frequencies" > $file".frequencies2"

# Renames the newly created frequency file.
mv $file".frequencies2" $file".frequencies"

# Adjusts scaffold positions in the frequency file.
awk 'BEGIN{scaffold="";pos=0;newpos=0} {if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}' $file".frequencies" > $file".frequencies.newpos"

# Compresses the TreeMix frequency file.
gzip $file".treemix.frq"
