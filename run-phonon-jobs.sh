#!/bin/bash

# List of target directories
directories=("ta4c3" "ta4c3o" "ta4c3o2") # ignore ta2co for now

# Loop through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: $dir"
    
    # Change to the directory
    cd $dir || continue
    
    # Run phonopy command
    echo "Running phonopy in $dir"
    phonopy -d --dim="3 3 1" && wait
    
    # Count number of POSCAR files and extract maximum number
    max_jobs=$(ls POSCAR-* 2>/dev/null | sed 's/POSCAR-//' | sort -n | tail -n1)
    
    # Loop through numbers 1 to max_jobs
    for i in $(seq -f "%03g" 1 $max_jobs); do
        echo "Processing job $i in $dir"
        
        # Create directory
        mkdir -p $i
        
        # Copy POSCAR file with corresponding number
        cp POSCAR-$i $i/POSCAR
        cp POTCAR $i/POTCAR
        # Copy other required files
        for file in KPOINTS INCAR fw.script; do
            cp ../$file $i/
        done
        
        # Change to subdirectory, submit job, and return
        (cd $i && sbatch fw.script)
    done
    
    # Return to original directory
    cd ..
done