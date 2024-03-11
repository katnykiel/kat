from pymatgen.core.structure import Structure

def round_partial_occupancies(struct):
    """
    Round the occupancies of the sites in the structure to the nearest value that is an integer multiple of 1/num_sites
    """

    # Create a dictionary to store the mapping from occupancy_dict to the number of sites
    occupancy_mapping = {}

    # For each site in the structure, get the occupancy
    for site in struct.sites:
        occupancy_dict = site.species.as_dict()

        # Convert the occupancy_dict to a hashable type
        occupancy_key = tuple(occupancy_dict.items())

        # Check if the occupancy_dict already exists in the mapping dictionary
        if occupancy_key in occupancy_mapping:
            occupancy_mapping[occupancy_key] += 1
        else:
            occupancy_mapping[occupancy_key] = 1

    # For each site in the structure, get the occupancy
    for site in struct.sites:
        occupancy_dict = site.species.as_dict()

        # Convert the occupancy_dict to a hashable type
        occupancy_key = tuple(occupancy_dict.items())

        # Get the number of sites with the same occupancy_dict
        total_sites = occupancy_mapping[occupancy_key]
        
        # For each value in the dictionary, num_sites*value is not an integer, round it to the nearest value that is
        for key in occupancy_dict:
            occupancy_dict[key] = round(total_sites*occupancy_dict[key])/total_sites
        site.species = occupancy_dict
    return struct