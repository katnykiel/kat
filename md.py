from ase.io import read, write
from fireworks import LaunchPad
import glob
import os

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

get_frames(65257)