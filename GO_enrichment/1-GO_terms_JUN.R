##topGO 4 Species
##https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
##GET ZEBRA FINCH GENES ID#
#grep -oP '(?<=zfinch_gene-).*?(?=;)' PIN_caididates_upper.txt | sort -u > PIN_genes_ZEFIID.txt
##REMOVE in EACH LINE EVEYTHING AFTER FIRST COMA
#sed -i.bak 's/,.*$//' genes_annotation.gff
rm(ls=list())
source("https://bioconductor.org/biocLite.R")
BiocManager::install("biocLite")
source("https://bioconductor.org/biocLite.R")

setwd("/mnt/DATA/B_JUN_reseq/L_GO_enrichment/B_inversion_chr5")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")

library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)

# create GO db for genes to be used using biomaRt - please note that this takes a while
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='tguttata_gene_ensembl', host="https://www.ensembl.org")

# Read all junco genes
exp_dataJUN= read.table('Jhye_gene_names.txt', header=TRUE)
bg_genesJUN=as.character(exp_dataJUN[,1])

# Convert junco genes to Go terms
go_idsJUN= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genesJUN_1, mart=db,useCache = FALSE)

# Read candidate genes
INV<- read.table("genes_list_edited.txt", header=T)
INV=as.character(INV[,1])

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GOJUN=unstack(go_idsJUN[,c(1,2)])

# remove any candidate genes without GO annotation
keepJUN = INV %in% go_idsJUN[,2]
keepJUN = which(keepJUN==TRUE)
INV=INV[keepJUN]

# make named factor showing which genes are of interest
geneListINV=factor(as.integer(bg_genesJUN %in% INV))
names(geneListINV)= bg_genesJUN

myGOdataINV <- new("topGOdata", description="My project", ontology="BP", allGenes=geneListINV,  annot = annFUN.gene2GO, gene2GO = gene_2_GOJUN)
myGOdataINV

#The 'ontology' argument can be 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component).
#The 'description' argument has a short description of your project.
#The 'allGenes' argument specifies all the genes in the 'gene universe', and which of those are  your genes of interest.
#The 'annot' argument tells topGO how to map genes to GO annotations. 'annot'=annFUN.gene2GO means that the user provides gene-to-GO annotations, and we specify here that they are in object 'geneID2GO' (see above for how this was created).
#An optional argument is 'nodeSize': this is used to prune the GO hierarchy, eg. nodesize=10 prunes the GO hierarchy, to remove terms which have less than 10 annotated genes

#The list of genes of interest can be accessed using the method sigGenes():

sgINV <- sigGenes(myGOdataINV)
str(sgINV)
numSigGenes(myGOdataINV)

##Performing enrichment tests
resultFisherINV <- runTest(myGOdataINV, algorithm="weight01", statistic="fisher")
resultFisher_classicINV <- runTest(myGOdataINV, algorithm="classic", statistic="fisher")


#Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account,
#so each GO term is tested independently. (Usually in fact we want to take the GO hierarchy into account,
#so would use algorithm='weight01' for example.)

##Results summary
resultFisherINV
resultFisher_classicINV

#We can list the top ten significant results found:
allGOINV=usedGO(myGOdataINV)
#allResJUN_3 <- GenTable(myGOdataJUN_3, classicFisher = resultFisherJUN_3, orderBy = "weightFisher", numChar=100, ranksOf = "classicFisher",  topNodes=length(allGOJUN_3))
allResINV <- GenTable(myGOdataINV, classicFisher = resultFisherINV, orderBy = "weightFisher", numChar=100, ranksOf = "classicFisher",  topNodes=10)
allResINV

