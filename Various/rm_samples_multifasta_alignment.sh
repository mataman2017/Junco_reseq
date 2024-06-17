#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <multifasta_alignment_file> <samples_to_remove_file> <output_alignment_file>"
    exit 1
fi

# Assign input arguments to variables
alignment_file="$1"
samples_to_remove_file="$2"
output_alignment_file="$3"

# Check if input files exist
if [ ! -f "$alignment_file" ]; then
    echo "Error: Multifasta alignment file '$alignment_file' not found."
    exit 1
fi

if [ ! -f "$samples_to_remove_file" ]; then
    echo "Error: Samples to remove file '$samples_to_remove_file' not found."
    exit 1
fi

# Create a temporary file to store the filtered alignment
temp_alignment_file=$(mktemp)

# Read the samples to remove into an array
mapfile -t samples_to_remove < "$samples_to_remove_file"

# Flag to indicate whether a sequence should be removed or not
remove_sequence=false

# Loop through the alignment file
while IFS= read -r line; do
    # If the line is a header line
    if [[ "$line" == ">"* ]]; then
        # Extract sample name from the header line
        sample=$(echo "$line" | sed 's/^>//')
        # Check if the sample should be removed
        if [[ " ${samples_to_remove[@]} " =~ " $sample " ]]; then
            remove_sequence=true
        else
            remove_sequence=false
            echo "$line" >> "$temp_alignment_file"
        fi
    elif ! $remove_sequence; then
        # Write sequence data to the temporary alignment file if not flagged for removal
        echo "$line" >> "$temp_alignment_file"
    fi
done < "$alignment_file"

# Move the filtered alignment to the specified output file path
mv "$temp_alignment_file" "$output_alignment_file"

echo "Samples removed successfully. Output saved to: $output_alignment_file"
