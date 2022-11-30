#!/usr/bin/env Rscript
# Calculate Xpehh and iHS in this case, for the populations DOR and CAN
# 14/03/22
rm(list = ls())

setwd("/mnt/6TB/1-Juncos/4-PopIndexes/2-Rehh/A_cmd/236")

library(rehh)
library(tidyverse)

#######################  ONE BY ONE CHR ##############

args = commandArgs(trailingOnly=TRUE) 

### Set the variables. Change the name of the pops. 
P1<-"DOR"
P2<-"CAN"

# Create the OUT dirs if they do not exist
O_P1_name<-paste("/DATA/2-WorkData/juncos_reseq/4-PopIndexes/2-Rehh/236/B_OUT/",P1,"/",sep="")
O_P2_name<-paste("/DATA/2-WorkData/juncos_reseq/4-PopIndexes/2-Rehh/236/B_OUT/",P2,"/",sep="")
O_name<-paste("/mnt/6TB/1-Juncos/4-PopIndexes/2-Rehh/C_outfiles/",P1,P2,"/",sep="") #make sure that the directory of destination has its name ordered this way
dir.create(O_P1_name,recursive=TRUE)
dir.create(O_P2_name,recursive=TRUE)
dir.create(O_name,recursive=TRUE)

# Location of the infiles PLEASE MAKE SURE THE PATH AND THE NAME ARE CORRECT
I_P1<-paste("/DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/",P1,"/",P1,"_",args[1],"-var_phased.vcf.gz",sep="") 
I_P2<-paste("/DATA/2-WorkData/juncos_reseq/2-FilteredBiallelicSNPs/236/phased/",P2,"/",P2,"_",args[1],"-var_phased.vcf.gz",sep="")
# Destination and name of the outfile
O_P1<-paste(O_P1_name,P1,"_",args[1],"_IHS.txt",sep="")
O_P2<-paste(O_P2_name,P2,"_",args[1],"_IHS.txt",sep="")
O<-paste(O_name,P1,P2,"_",args[1],"_XPEHH.txt",sep="")

O_P1_scan<-paste("/DATA/2-WorkData/juncos_reseq/4-PopIndexes/2-Rehh/236/B_OUT/",P1,"/",P1,"_",args[1],"_scan_hh.txt",sep="")
O_P2_scan<-paste("/DATA/2-WorkData/juncos_reseq/4-PopIndexes/2-Rehh/236/B_OUT/",P2,"/",P2,"_",args[1],"_scan_hh.txt",sep="")

cat("All variables for", args[1], "of", P1,"have been assigned")
### READ DATA

P1_hh <- data2haplohh(hap_file = I_P1, remove_multiple_markers=TRUE, polarize_vcf = FALSE)
P2_hh <- data2haplohh(hap_file = I_P2, remove_multiple_markers=TRUE, polarize_vcf = FALSE)

cat("Data for", args[1], "of", P1,"has been read")

### SCAN (long part) This computes Ihs per population. 
P1_scan <- scan_hh(P1_hh, polarized = FALSE)
P2_scan <- scan_hh(P2_hh, polarized = FALSE)

cat("SCAN, for", args[1], "of", P1,"has been performed")
### write table from SCAN so that this just needs to be done once
write.table(P1_scan, O_P1_scan)
#write.table(P2_scan, O_P2_scan)

cat("All Scans saved for", args[1], "of", P1) # I save the scans so that I can compute XPEHH really fast on different population comparisons.

## Read the saved Scans in case they have already been computed
# P1_scan<-read.table(O_P1_scan)
# P2_scan<-read.table(O_P2_scan)

### PERFORM ihs
P1_ihs <- ihh2ihs(P1_scan, freqbin = 1)
P2_ihs <- ihh2ihs(P2_scan, freqbin = 1)

### write out Ihs POP1
write.table(P1_ihs, O_P1)

### write out Ihs POP2
write.table(P2_ihs, O_P2) 

cat("IHS calculated and saved for", args[1], "of", P1)

### perform xp-ehh
XPEHH <- ies2xpehh(P1_scan, P2_scan, popname1 = P1, popname2=P2, include_freq = T) 

## write out can-dor xpEHH
write.table(XPEHH, O)

#cat("XPEHH calculated and saved for", args[1], "of", P1,P2)
cat("END OF THE JOB FOR", args[1], "OF", P1)
