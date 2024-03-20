from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatgen.core.composition import Composition
from fireworks import LaunchPad
import json
import pandas as pd
from monty.json import MontyEncoder
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.computed_entries import ComputedEntry
from mp_api.client import MPRester
from collections import Counter
import os
import numpy as np
import plotly.graph_objects as go
from pymatgen.io.vasp import Outcar

# load the launchpad for querying
launchpad = LaunchPad.auto_load()

def get_convegence_plots(doc, energies = True, forces = False, energy_type = "e_fr_energy"):
    """
    Given an atomate2 document, returns the convergence plot as a plotly figure.
    """
    # get the convergence data
    run_data = doc["output"]["calcs_reversed"][0]["output"]

    # Extract the energies, lining up electronic and ionic steps
    i_energies = [i_step[energy_type] for i_step in run_data['ionic_steps']]
    e_steps_sets = [i_step['electronic_steps'] for i_step in run_data['ionic_steps']]
    e_energies_sets = [[e_step[energy_type] for e_step in e_steps] for e_steps in e_steps_sets]
    e_count = [len(e) for e in e_energies_sets]
    i_steps = [sum(e_count[0:n+1])-1 for n in range(len(e_count))]
    e_energies = [item for sublist in e_energies_sets for item in sublist]
    e_steps = np.linspace(1,len(e_energies),len(e_energies))

    # Create and save a plotly figure with the energy at each ionic and electronic step
    fig_energies = go.Figure()
    fig_energies.add_trace(go.Scatter(x = e_steps, y = e_energies, name = 'electronic'))
    fig_energies.add_trace(go.Scatter(x = i_steps, y = i_energies, name = 'ionic'))

    scaling_factor = 1.001
    fig_energies.update_layout(
        # title = f'Energy Convergence: {doc["output"]["formula_pretty"]} {doc["output"]["task_type"]}',
        xaxis_title = 'electronic steps',
        yaxis_title = 'energy (eV)',
        yaxis_range = [min(i_energies)*scaling_factor,max(i_energies)/scaling_factor],
        legend=dict(x=.97, y=.97, xanchor='right', yanchor='top'),
        template="simple_white",
        font=dict(size=18),
        width = 800,
        height = 600,
        margin=dict(l=50, r=50, t=50, b=50)
        
        )

    fig_energies.update_xaxes(mirror=True, showgrid=True)
    fig_energies.update_yaxes(mirror=True, showgrid=True)
    

    return fig_energies

def check_convergence(fw_id):
    """
    Given a firework ID, create a plot of the convergence data and return the plotly figure.
    """

    # Get the current run directory, and any past runs from archived_launches
    fw = launchpad.get_fw_by_id(fw_id)
    run_dir = fw.launches[-1].launch_dir
    archived_launches = fw.archived_launches
    archived_run_dirs = [launch.launch_dir for launch in archived_launches if len(archived_launches) > 0]
    run_dirs = archived_run_dirs + [run_dir]

    energies = []
    for run_dir in run_dirs:
        # load a task doc from the run directory using pymatgen outcar
        # if OUTCAR.gz exists, unzip it
        if os.path.exists(run_dir + "/OUTCAR.gz"):
            os.system(f"gunzip {run_dir}/OUTCAR.gz")
        outcar = Outcar(run_dir + "/OUTCAR")
        outcar.read_pattern({"energy": r"energy\s*\(sigma->0\)\s*=\s*(-?\d+\.\d+)"})
        energies.append([float(e[-1]) for e in outcar.data["energy"]])

    # flatten energies
    energies = [item for sublist in energies for item in sublist]
    # Get the steps for each energy
    e_steps = np.linspace(1,len(energies),len(energies))
    fig_energies = go.Figure()
    fig_energies.add_trace(go.Scatter(x = e_steps, y = energies, name = 'electronic'))

    fig_energies.update_layout(
    # title = f'Energy Convergence: {doc["output"]["formula_pretty"]} {doc["output"]["task_type"]}',
    xaxis_title = 'electronic steps',
    yaxis_title = 'energy (eV)',
    # yaxis_range = [min(i_energies)*scaling_factor,max(i_energies)/scaling_factor],
    template="simple_white",
    font=dict(size=18),
    margin=dict(l=50, r=50, t=50, b=50)
        )

    fig_energies.update_xaxes(mirror=True, showgrid=True)
    fig_energies.update_yaxes(mirror=True, showgrid=True)

    fig_energies.show()
    
    return fig_energies

def convert_old_structure_to_new(structure):
    """
    Convert a pymatgen structure to a dictionary with the new format.
    """

    # write structure to a POSCAR file
    structure.to("POSCAR", fmt="POSCAR")

    # read the POSCAR file and convert to a dictionary
    structure = Structure.from_file("POSCAR")

    # remove the POSCAR file
    os.remove("POSCAR")
    return structure

def continue_firework(fw_id):
    """Given a firework ID of a job that hit the waltime, update the firework spec with the last output and rerun the firework.

    Args:
        fw_id (int): id of firework which hit the walltime
    """    

    # Get the last geometric output from the run directory
    fw = launchpad.get_fw_by_id(fw_id)
    struct = Structure.from_file(fw.launches[-1].launch_dir + "/CONTCAR")

    # Update the firework spec in the job store with the new structure
    launchpad.update_spec([fw_id], {"_tasks.0.job.function_args": [struct.as_dict()]})

    # Re-run the firework # TODO: check why this sometimes runs in the same directory
    launchpad.rerun_fw(fw_id)

    return

def continue_workflows(pattern = "continue_workflow"):
    """
    Continue workflows that match the given pattern.

    Args:
        pattern (str): The pattern to match against workflow names. Defaults to "continue_workflow".

    Returns:
        None
    """

    # query for all workflows that match the pattern
    pattern = f".*{pattern}.*"
    workflows = launchpad.workflows.find({"name": {"$regex": pattern}})

    for workflow in workflows:
        
        # remove workflows where the status of all fireworks is completed
        if Counter(workflow["fw_states"].values())["COMPLETED"] == len(workflow["fw_states"]):
            continue

        # get all the fizzled fireworks
        fizzled_fw = [int(id) for id, state in workflow["fw_states"].items() if state == "FIZZLED"]

        for id in fizzled_fw:
            try:
                continue_firework(id)
                print(f"Continued firework {id}")
            except:
                print(f"Could not continue firework {id}")

def save_docs_to_df(docs, file_name="docs.json"):
    """
    Save a list of documents to a JSON file.

    Args:
        docs (list): The list of atomate2 documents to be saved.
        file_name (str, optional): The name of the output JSON file. Defaults to "docs.json".
    """
    # turn the results into a dataframe
    df = pd.DataFrame(docs)

    # serialize the results using MontyEncoder
    serialized_results = df.to_dict(orient="records")

    # save results to .json file
    with open(file_name, "w") as f:
        json.dump(serialized_results, f, cls=MontyEncoder)

def get_convex_hulls(df, show_fig = False, write_results = False):
    """
    Given a dictionary `dataframe` containing a list of atomate2 docs, generates a set of convex hull diagrams using pymatgen.
    """

    # Get the list of chemical systems from the doc df
    df["chemical_system"] = df["output"].apply(
    lambda x: Composition.from_dict(x["composition"]).chemical_system)

    # Get sub_dfs grouped by unique_chemsys
    sub_dfs = [df[df["chemical_system"] == chemsys] for chemsys in df["chemical_system"].unique()]

    results = {'structures':[],'E_above_hull':[], "stable_structures":[]}

    figs = []
    for sub_df in sub_dfs:

        # Get the elements in chemical_system
        elements = sub_df["chemical_system"].iloc[0].split("-")

        # Sort the elements so that C/N are last
        elements.sort(key=lambda x: x not in ["C", "N"])

        # Get the list of computed entries from the doc df
        entries =sub_df["output"].apply(lambda x: x["entry"])
        computed_entries = [ComputedEntry.from_dict(e) for e in entries]

        # Query MP for the additional entries
        pmg_mapi_key = os.environ.get("PMG_MAPI_KEY")
        with MPRester(pmg_mapi_key) as mpr:

            # Obtain ComputedStructureEntry objects
            mp_entries = mpr.get_entries_in_chemsys(elements=elements) 
        
        # Create a PhaseDiagram object from the entries
        pd = PhaseDiagram(computed_entries + mp_entries)
        stable_entries = pd.stable_entries

        # get the energy above hull from phase diagram
        E_above_hull = [pd.get_e_above_hull(entry) for entry in computed_entries] # type: ignore

        # get the structures
        structures = sub_df["output"].apply(lambda x: x["structure"]).tolist()

        # add the energy above hull and structures to results dict
        results["E_above_hull"].append(E_above_hull)
        results["structures"].append(structures)
        # results["stable_structures"].append([entry.structure for entry in stable_entries if "structure" in entry.as_dict().keys()])

        for i, entry in enumerate(computed_entries):
            uuid = sub_df.iloc[i]["uuid"]
            df.loc[df["uuid"] == uuid, "energy_above_hull"] = E_above_hull[i]

        if len(elements) > 4:
            print(f"Too many elements to plot: {elements}")
            continue
        
        # Create a PDPlotter object from the PhaseDiagram object
        plotter = PDPlotter(pd)

        # Generate the convex hull diagram
        fig = plotter.get_plot(highlight_entries=computed_entries) # type: ignore

        # write fig to file
        fig.write_image(f"convex_hull_{'-'.join(elements)}.png", scale=3)# type: ignore
        
        if show_fig:
            fig.show()# type: ignore

        figs.append(fig)

    if write_results:
        # serialize the results using MontyEncoder
        serialized_results = df.to_dict(orient="records")

        # save results to .json file
        with open("convex_hull_docs.json", "w") as f:
            json.dump(serialized_results, f, cls=MontyEncoder)

        # write tailored results to a .json file
        with open("convex_hull_results.json", "w") as f:
            json.dump(results, f, cls=MontyEncoder)
    return figs