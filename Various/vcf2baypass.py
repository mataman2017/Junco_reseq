# Script to convert vcf or gzvcf into the baypass format. 
# Launching: $ python vcf_to_genotype.py input.vcf population.txt

import sys
import gzip

# Function to read VCF file
def read_vcf_file(vcf_file):
    if vcf_file.endswith('.gz'):
        with gzip.open(vcf_file, 'rt') as f:
            lines = [line.strip() for line in f if not line.startswith('#')]
    else:
        with open(vcf_file, 'r') as f:
            lines = [line.strip() for line in f if not line.startswith('#')]
    genotypes = []
    for line in lines:
        fields = line.split('\t')
        genotypes.append([int(x) for x in fields[9:]])
    return genotypes

# Function to write genotyping data file
def write_genotyping_data_file(genotypes, pop_assignments, outfile):
    npop = len(pop_assignments)
    nsnp = len(genotypes)
    with open(outfile, 'w') as f:
        for i in range(nsnp):
            counts = []
            for j in range(npop):
                pop = pop_assignments[j]
                if pop == '.':
                    counts.append('0 0')
                else:
                    gt = genotypes[i][j]
                    if gt == 0:
                        counts.append('2 0')
                    elif gt == 1:
                        counts.append('1 1')
                    elif gt == 2:
                        counts.append('0 2')
                    else:
                        counts.append('0 0')
            f.write(' '.join(counts) + '\n')

# Main program
if __name__ == '__main__':
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print(f'Usage: {sys.argv[0]} <vcf_file> <pop_assignments_file> <outfile>')
        sys.exit(1)
    vcf_file = sys.argv[1]
    pop_assignments_file = sys.argv[2]
    outfile = sys.argv[3]

    # Read population assignments
    with open(pop_assignments_file, 'r') as f:
        pop_assignments = [line.strip() for line in f]

    # Read genotypes from VCF file
    genotypes = read_vcf_file(vcf_file)

    # Write genotyping data file
    write_genotyping_data_file(genotypes, pop_assignments, outfile)
