#!/bin/bash

# Define the input and output directories
input_dir="./out"
output_dir="./out_results"

# Check if the output directory exists, create it if not
if [ ! -d "$output_dir" ]; then
  mkdir "$output_dir"
fi

# Loop through all PDB files in the input directory
for pdb_file in "$input_dir"/*.pdb; do
  # Check if the file exists
  if [ -f "$pdb_file" ]; then
    # Get the base name of the PDB file without the extension
    base_name=$(basename "$pdb_file" .pdb)

    # Print a message to show the progress
    echo "Processing $pdb_file..."

    # Run prank on the current PDB file and save the results in the output directory
    ./prank predict -f "$pdb_file" -o "$output_dir/$base_name"
    
    # Print a message to show that the processing has finished for the current file
    echo "Finished processing $pdb_file. Results saved in $output_dir/$base_name."
    echo "--------------------------------------------------------"
  fi
done

# Print a message to show that the script has finished running
echo "All PDB files processed. Results saved in $output_dir."
