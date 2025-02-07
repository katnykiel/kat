#!/bin/bash

# List of target directories
directories=("ta4c3" "ta4c3o" "ta4c3o2") # ignore ta2co for now

# Loop through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: $dir"
    
    # Change to the directory
    cd $dir || continue
    
    echo "Running born in $dir"
    
    mkdir -p born_supercell
    cp ../INCAR_born born_supercell/INCAR
    cp ../KPOINTS born_supercell/
    cp ../fw_born.script born_supercell/
    cp POTCAR born_supercell/
    cp SPOSCAR born_supercell/POSCAR

    cd born_supercell

    # run vasp job
    sbatch fw_born.script
    cd ..
    
    # Return to original directory
    cd ..
done