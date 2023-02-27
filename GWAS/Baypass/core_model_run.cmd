#!/bin/bash/

GFILE=/mnt/DATA/B_JUN_reseq/J_GWAS/A_baypass/A_PDC_61/baypass_infile.geno

baypass -gfile $GFILE -outprefix PDC_61 -nthreads 20
