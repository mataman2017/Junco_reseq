import sys
import argparse
import gzip
import re
import subprocess

"""
usage: python polarizeVCFbyOutgroup [-h] [-vcf VCF] [-out OUT] [-ind IND] [-keep] [-add] [-gz]

Switch REF and ALT allele of a vcf file if the specified individual is homozygous ALT.

This script was modified by Javi Sala-Garcia (javier.sala@mncn.csic.es),
from Museo Nacional de Ciencias Naturales (MNCN-CSIC) at JUNE 2023,
from a script provided by Krisian Ullrich's personal GitHub.

Prior to running the script, first remove info from the VCF:
bcftools annotate -x INFO,FORMAT,^INFO/FORMAT/GT DATA.vcf.gz | bgzip -c > DATA_GT.vcf.gz

-vcf can be gzipped
-out complete name of the outfile
-ind IT IS NOT THE NAME OF THE SAMPLE BUT THE POSITIONAL NUMBER OF THE SAMPLE IN THE VCF
-keep
-add it adds the field AA, indicating the ancestral variant
-gz to compress outfile

...
"""

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in sorted(rep_dict,key=len,reverse=True)]), flags=re.DOTALL)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

def parse_lines(fin, fou, ind, keep, add):
    switchcount = 0
    removecount = 0
    outmissing = 0
    totalcount = 0
    for line in fin:
        line = line.rstrip('\r\n')
        if line.startswith('#'):
            if add and line.startswith('##INFO'):
                fou.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n')
            fou.write(line + '\n')
        else:
            totalcount += 1
            linesplit = line.split('\t')
            format_fields = linesplit[8].split(':')
            genotypes = [fields.split(':') for fields in linesplit[9:]]
            ind_genotype = genotypes[ind][0]

            if ind_genotype == './.' or ind_genotype == '.|.':
                outmissing += 1
            elif ind_genotype == '0/0' or ind_genotype == '0|0':
                if add:
                    linesplit[7] = 'AA=' + linesplit[3] + ';' + linesplit[7]
                fou.write('\t'.join(linesplit) + '\n')
            elif ind_genotype == '0/1' or ind_genotype == '0|1' or ind_genotype == '1/0' or ind_genotype == '1|0':
                removecount += 1
                if keep:
                    fou.write(line + '\n')
            elif ind_genotype == '1/1' or ind_genotype == '1|1':
                switchcount += 1
                REF = linesplit[3]
                ALT = linesplit[4]
                linesplit[3] = ALT
                linesplit[4] = REF
                changed_genotypes = [multiple_replace(gt[0], {'0':'1', '1':'0'}) + ':' + ':'.join(gt[1:]) if gt[1:] else multiple_replace(gt[0], {'0':'1', '1':'0'}) for gt in genotypes]
                if add:
                    linesplit[7] = 'AA=' + linesplit[3] + ';' + linesplit[7]
                fou.write('\t'.join(linesplit[:8] + format_fields) + '\t' + '\t'.join(changed_genotypes) + '\n')

    print('Parsed ' + str(totalcount) + ' sites.')
    print('Removed ' + str(outmissing) + ' sites due to missing allele info in switch individual.')
    if keep:
        print('Kept ' + str(removecount) + ' sites with undefined ancestral state.')
    else:
        print('Removed ' + str(removecount) + ' sites with undefined ancestral state.')
    print('Switched REF and ALT allele for ' + str(switchcount) + ' sites.')


def main():
    parser = argparse.ArgumentParser(prog='polarizeVCFbyOutgroup', description='Switch REF and ALT allele of a vcf file if the specified individual is homozygous ALT.')
    parser.add_argument('-vcf', help='specify vcf input file')
    parser.add_argument('-out', help='specify output file')
    parser.add_argument('-ind', type=int, help='specify index of individual to polarize')
    parser.add_argument('-keep', action='store_true', help='keep sites with undefined ancestral state')
    parser.add_argument('-add', action='store_true', help='add ancestral allele INFO field')
    parser.add_argument('-gz', action='store_true', help='compress output using gzip and index with bcftools')
    args = parser.parse_args()

    vcf_file = args.vcf
    output_file = args.out
    ind_index = args.ind
    keep_undefined = args.keep
    add_ancestral = args.add
    compress_output = args.gz

    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as fin:
        with open(output_file, 'w') as fou:
            parse_lines(fin, fou, ind_index, keep_undefined, add_ancestral)

    if compress_output:
        subprocess.run(['bgzip', output_file])
        subprocess.run(['bcftools', 'index', output_file + '.gz'])

if __name__ == '__main__':
    main()
