#/bin/bash

## Genome-environment association analysis using BayPass
## Javi Sala 2024

# Generate input file from VCF
zcat JUN_205_DY.vcf.gz | perl vcf2baypass.pl pop.txt JUN_DY

# Run BayPass under the core model mode to generate covariance matrix
/home/milab/bin/baypass -npop 14 -gfile JUN_DY.txt -outprefix all -nthreads 20

# Run BayPass under the standard covariate model using importance sampling (IS) estimator
i_baypass -npop 20 -gfile GSD_RAD.HA412.good_snps2.baypass.txt -efile soil.baypass.txt -scalecov -omegafile all_mat_omega.out -outprefix all_soilGEA -nthreads 20

# Produce POD samples with 1,000 SNPs and run BayPass
Rscript pod.R
i_baypass -npop 20 -gfile G.GSD_RAD.HA412.good_snps2.baypass.pod -efile soil.baypass.txt -scalecov -omegafile all_mat_omega.out -outprefix pod_soilGEA -nthreads 20
