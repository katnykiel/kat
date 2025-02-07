#!/bin/bash

# List of target directories
directories=("ta4c3" "ta4c3o" "ta4c3o2") # ignore ta2co for now

# Loop through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: $dir"
    
    # Change to the directory
    cd $dir || continue
    
    # Count number of POSCAR files and extract maximum number
    max_jobs=$(ls POSCAR-* 2>/dev/null | sed 's/POSCAR-//' | sort -n | tail -n1)
    
    # Loop through numbers 1 to max_jobs
    for i in $(seq $max_jobs); do
        folder=$(printf "%03d" $i)
        cp "$(printf '%03d' $i)/vasprun.xml" "vasprun-$(printf '%03d' $i).xml"
    done

    # Run phonopy on all vasprun files
    phonopy -f vasprun-*.xml
    wait

    phonopy --dim="3 3 1" --fc-symmetry --mesh="1 1 1" --eigenvectors
    wait
    
    # Return to original directory
    cd ..
done