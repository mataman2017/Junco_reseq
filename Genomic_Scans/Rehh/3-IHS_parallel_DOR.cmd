#!/bin/bash

#Script to launch an R Rehh script for all the Chr at a time for the population dorsalis (DOR)
### set the variables
LIST=/DATA/2-WorkData/juncos_reseq/4-PopIndexes/2-Rehh/236/A_cmd/A_list/chr.list
#### command
cat $LIST | parallel Rscript R/Rehh_IHS_DOR.R   # This should be changed to whatever population you want to run the analysis
