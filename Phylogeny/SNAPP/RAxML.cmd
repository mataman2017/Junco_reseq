#!/bin/bash

# Using a filtered and subseted VCF file to calculate a tree to use as input for SNAPP

VCF=paldorcan_nohyb_4_biall_maf005_hwe00001_mac1_q30_dp5_50_miss075.recode.vcf.gz
NAME=${VCF%.recod*}

### CONVERT TO RAXML-COMPATIBLE PHYLIP ###

vcf2phylip.py -i $VCF --output-prefix $NAME -o junco_bm76_palliatus_mexic

ascbias.py -p ${NAME}.min4.phy -o {NAME}_var.min4.phy #

### RUN RAXML ###

raxml -s ${NAME}_var.min4.phy \
	-d \
	-p 12345 \
	-x 12345 \
	-f a -# 100 \
	-m ASC_GTRGAMMA \
	--asc-corr=lewis \
	-o junco_bm76_palliatus_mexic \
	-n paldorcan_single \
	-T 8
