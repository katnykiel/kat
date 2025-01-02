from pymatgen.core.structure import Structure

def round_partial_occupancies(struct: Structure) -> Structure:
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

    site_updates = []
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

        # check that the sum of the occupancies is equal to 1, and if not adjust the occupancy of the last element
        if sum(occupancy_dict.values()) != 1:
            occupancy_dict[list(occupancy_dict.keys())[-1]] += 1 - sum(occupancy_dict.values())
        site_updates.append(occupancy_dict)

    # Update the occupancies of the sites in the structure
    for site, updated_occupancy in zip(struct.sites, site_updates):
        site.species = updated_occupancy

    return struct