import argparse
from Bio import SeqIO

# List of sample names you want to keep
samples_to_keep = ["AIK4", "AIK5", "AIK7", "ALT1", "ALT2", "BAI1", "BAI2", "CAN1", "CAN2", "CAN3", "CAR12", "CAR13", "CAR2", "DOR13", "DOR15", "DOR18", "DOR32", "DOR33", "DOR34", "FUL1", "FUL2", "HYE11", "HYE13", "HYE9", "INS5", "INS8", "INS9", "MEA1", "MEA13", "MEA17", "MON13", "MON19", "MON9", "ORE11", "ORE12", "ORE9", "PAL11", "PAL16", "PAL6", "PHA13", "PHA2", "PHA6", "PIN4", "PIN7", "PIN8", "THU4", "THU5", "THU7", "TOW10", "TOW3", "TOW5", "VUL1", "VUL2"]

def main(input_file, output_file, samples_to_keep):
    # Create a dictionary to store sequences for each sample
    sample_sequences = {}

    # Iterate through the input FASTA file and filter sequences
    for record in SeqIO.parse(input_file, "fasta"):
        sample_name = record.id.split("|")[0]  # Adjust this based on your FASTA format

        if sample_name in samples_to_keep:
            if sample_name not in sample_sequences:
                sample_sequences[sample_name] = []
            sample_sequences[sample_name].append(str(record.seq))

    # Write the filtered sequences to the output file
    with open(output_file, "w") as output_handle:
        for sample_name, sequences in sample_sequences.items():
            combined_sequence = "".join(sequences)
            output_handle.write(f">{sample_name}\n{combined_sequence}\n")

    print(f"Filtering completed. Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and concatenate sequences from a FASTA file.")
    parser.add_argument("input_file", help="Input FASTA file to process")
    parser.add_argument("output_file", help="Output FASTA file to save concatenated sequences")

    args = parser.parse_args()

    main(args.input_file, args.output_file, samples_to_keep)
