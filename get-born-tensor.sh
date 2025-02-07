#!/bin/bash

# List of target directories
directories=("ta4c3" "ta4c3o" "ta4c3o2") # ignore ta2co for now

# Loop through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: $dir"
    
    # Change to the directory
    cd $dir || continue
    
    echo "Getting born in $dir"

    cd born || continue
    phonopy-vasp-born > BORN

    cp BORN ..
    cd ..
    
    # Return to original directory
    cd ..
done