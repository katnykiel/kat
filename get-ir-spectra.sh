#!/bin/bash

# List of target directories
directories=("ta4c3" "ta4c3o" "ta4c3o2") # ignore ta2co for now

# Loop through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: $dir"
    
    # Change to the directory
    cd $dir || continue
    
    echo "Getting ir spectra in $dir"

    phonopy-ir --linewidth=16.5 --spectrum-range="0 1200"

    # Return to original directory
    cd ..
done