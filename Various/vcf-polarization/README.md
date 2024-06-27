# VCF Polarization Pipeline

This pipeline polarizes a VCF file by switching REF and ALT alleles based on a specified individual's genotype.

## Prerequisites

- Python 3.x
- bcftools


## Usage

Prior to running the polarization script, annotate and modify your VCF file using `bcftools`:
```bash
bcftools annotate -x INFO,FORMAT,^INFO/FORMAT/GT JUN_227_VAR_nomiss.vcf.gz | bgzip -c > JUN_227_VAR_nomiss_GT.vcf.gz


This step removes unnecessary INFO and FORMAT fields, leaving only the genotype information needed for polarization.

Usage
To polarize a VCF file, use the following command:

python polarizeVCFbyOutgroup.py -vcf JUN_227_VAR_nomiss_GT.vcf.gz -out output.vcf.gz -ind 0 -add -keep -gz

Options:
-vcf: Specify the input VCF file (post bcftools annotate).
-out: Specify the output file name.
-ind: Index (number) of the individual to polarize.
-keep: Keep sites with undefined ancestral state.
-add: Add ancestral allele INFO field.
-gz: Compress output using gzip and index with bcftools.


### Notes:

- Ensure that `bcftools` and `bgzip` are installed and accessible in your environment.
- The `-x INFO,FORMAT,^INFO/FORMAT/GT` flag removes all INFO and FORMAT fields except for the genotype (GT) field, which is necessary for the polarization script.
- After running `bcftools annotate`, the output VCF file (`JUN_227_VAR_nomiss_GT.vcf.gz`) should be used as input (`-vcf`) for the polarization script.

By adding these instructions, users will understand the necessary preprocessing step before using your polarization pipeline. This ensures clarity and completeness in your GitHub repository documentation.

