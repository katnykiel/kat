import subprocess
from ase.io import read, write
from fireworks import LaunchPad
import glob
import os
from pymatgen.core import Structure
from pymatgen.io.vasp import Xdatcar
import re
from pymatgen.phonon.dos import PhononDos, CompletePhononDos
from pymatgen.phonon.plotter import PhononDosPlotter
from pymatgen.core import Element

lpad = LaunchPad.auto_load()

def get_frames(fw_id):
    """
    Given an atomate2 fw_id, write the frames to a file
    """

    # get the Firework document from the LaunchPad and extract file path
    fw = lpad.get_fw_by_id(fw_id)
    file_path = fw.launches[-1].launch_dir

    # look for a .traj file at this location, return if none is found
    traj_files = glob.glob(f"{file_path}/*.traj")
    if not traj_files:
        return
    trajectory = read(traj_files[0], index=':')

    # write each frame to a LAMMPS dump file
    for i, atoms in enumerate(trajectory):
        # Construct the filename for each frame
        filename = f"frame_{i}.lammps"
        # Write the current frame to a LAMMPS dump file
        path = os.path.join(file_path, filename)
        write(path, atoms, format='lammps-data') # type: ignore
    
    return file_path

def get_image(struct):
    """
    Given a structure, write the image to a file
    """
    # Write the structure to a LAMMPS dump file
    
def run_dos_pt(dir):
    """
    Given a VASP MD directory, run the VDOS calculation
    """

    # Convert XDATCAR to traj.xyz in the given directory
    xdatcar_path = os.path.join(dir, 'XDATCAR')
    traj = read(xdatcar_path, index=':')
    traj_xyz_path = os.path.join(dir, 'traj.xyz')
    write(traj_xyz_path, traj)

    # Extract information from the XDATCAR file using pymatgen
    xdatcar = Xdatcar(xdatcar_path)
    structure = xdatcar.structures[0]  # Use the first structure for cell parameters and masses
    
    # Extract cell parameters
    cell = ' '.join(' '.join(map(str, vec)) for vec in structure.lattice.matrix)
    
    # Extract temperature (assuming it's in the comment line)
    temperature = 300  # Default value, modify if temperature is specified elsewhere
    
    # Write input file
    input_file_path = os.path.join(dir, 'input')
    with open(input_file_path, 'w') as f:
        f.write(f"points = {len(traj)}\n")
        f.write(f"tau = {len(traj) * 0.001}\n")  # Assuming 1 fs timestep
        f.write(f"cell = {cell}\n")
        f.write(f"temperature = {temperature}\n")
        f.write("format = xyz\n")
        f.write("estimate_velocities = .true.\n")
        f.write("hs_formalism = lin\n")
        f.write("f_opt = .false.\n")
    
    # Write masses file using pymatgen
    masses_file_path = os.path.join(dir, 'masses')
    with open(masses_file_path, 'w') as f:
        for element in structure.composition.elements:
            f.write(f"{element.symbol} {element.atomic_mass}\n".replace(" amu", ""))
    
    # Write groups file
    groups_file_path = os.path.join(dir, 'groups')
    with open(groups_file_path, 'w') as f:
        num_sites = len(traj[0])  # Number of sites in the first frame
        f.write(f"{num_sites} {num_sites}\n")
        for i in range(1, num_sites + 1):
            f.write(f"1 1\n{i}\n")
    
    # Write supergroups file
    supergroups_file_path = os.path.join(dir, 'supergroups')
    with open(supergroups_file_path, 'w') as f:
        f.write(f"1-{num_sites}\n")

    # Change directory, run the dospt command, and change back
    subprocess.run(["DoSPT"], cwd = dir)

    complete_phonon_dos = read_dospt_output(structure, os.path.join(dir, 'dos'))

    return complete_phonon_dos

def read_dospt_output(structure, file_path):
    """
    Read the dospt output file and extract the frequency and density data for each group.
    Also read the total VDOS from the start of the file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    data = {}
    total_vdos = []
    current_group = None
    reading_total_vdos = False

    for line in lines:
        if line.startswith(" # Total density of states for all the groups involved"):
            reading_total_vdos = True
            continue
        if line.startswith(" # Density of states (in ps) for group"):
            current_group = int(re.search(r"group\s+(\d+)", line).group(1))
            data[current_group] = []
            reading_total_vdos = False
        elif reading_total_vdos:
            parts = line.split()
            if len(parts) >= 2:
                frequency = float(parts[0])
                density = float(parts[1])
                total_vdos.append((frequency, density))
        elif current_group is not None:
            parts = line.split()
            if len(parts) >= 2:
                frequency = float(parts[0])
                density = float(parts[1])
                data[current_group].append((frequency, density))

    # Create a PhononDos object from the total VDOS
    frequencies, densities = zip(*total_vdos)
    total_phonon_dos = PhononDos(frequencies, densities)

    # Create a CompletePhononDos object from the group data
    group_dos = {structure.sites[i-1]: [density for _, density in dos] for i, dos in data.items()}
    complete_phonon_dos = CompletePhononDos(structure, total_phonon_dos, group_dos)

    return complete_phonon_dos