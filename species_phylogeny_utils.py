import parse_midas_data

# Returns a copy of the list of species sorted phylogenetically
def sort_phylogenetically(species):

    
    # Load the species tree from the midas directory
    from Bio import Phylo
    newick_filename = parse_midas_data.midas_directory+"/species_tree.newick"
    tree = Phylo.read(newick_filename, 'newick')
    ordered_species_ids = [term.name for term in tree.get_terminals()]

    # the ordered species IDs are just the numbers at the end of the full species name
    # so we need to convert the allowed list
    id_name_map = {}
    for species_name in species:
        items = species_name.split("_")
        id = items[-1]
        id_name_map[id] = species_name
        
    # walk through the sorted species_id list and output the ones that we want
    ordered_species = []
    for species_id in ordered_species_ids:
        if species_id in id_name_map:
            ordered_species.append(id_name_map[species_id])
            
    return ordered_species


if __name__=='__main__':
    # Test that it works
    good_species = parse_midas_data.parse_good_species_list()
    ordered_good_species = sort_phylogenetically(good_species)
    
    good_species_set = set(good_species)
    ordered_good_species_set = set(ordered_good_species)
    
    print ordered_good_species
    print good_species_set - ordered_good_species_set
    print ordered_good_species_set - good_species_set
    
    