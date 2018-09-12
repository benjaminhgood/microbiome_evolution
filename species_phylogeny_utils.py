import parse_midas_data

def get_genus_name(species_name):
    return species_name.split("_")[0]
    
def get_taxonomy_map():

    species_genome_map = {}
    genome_species_map = {}
    
    file = open("%sspecies_info.txt" % parse_midas_data.midas_directory,"r")
    file.readline() # header
    for line in file:
        items = line.split()
        species_name = items[0].strip()
        genome = items[1].strip()
        species_genome_map[species_name] = genome
        genome_species_map[genome] = species_name
    file.close()
    
    species_taxonomy_map = {}
    
    file = open("%sgenome_taxonomy.txt" % parse_midas_data.midas_directory,"r")
    file.readline() # header
    for line in file:
        items =  line.split("\t")
        
        genome_id = items[0].strip()
        
        kingdom = items[3].strip()
        phylum = items[4].strip()
        class_name = items[5].strip()
        order = items[6].strip()
        family = items[7].strip()
        genus = items[8].strip()
        
        if genome_id in genome_species_map:
        
            species_name = genome_species_map[genome_id]
            species_taxonomy_map[species_name] = (kingdom,phylum,class_name,order,family,genus)
        
    return species_taxonomy_map
    
   
def sort_phylogenetically(species_list, first_entry="", second_sorting_attribute=[]):

    species_taxonomy_map = get_taxonomy_map()
    
    kingdom_order_map = {}
    order_kingdom_map = []
    
    phylum_order_map = {}
    order_phylum_map = []
    
    class_order_map = {}
    order_class_map = []
            
    order_name_order_map = {}
    order_order_name_map = []
    
    family_order_map = {}
    order_family_map = []
    
    genus_order_map = {}
    order_genus_map = []
    
    if first_entry!="":
        
        kingdom,phylum,class_name,order,family,genus = species_taxonomy_map[first_entry]
        
        kingdom_order_map[kingdom] = 0
        order_kingdom_map.append(kingdom)
        
        phylum_order_map[phylum] = 0
        order_phylum_map.append(phylum)
        
        class_order_map[class_name] = 0
        order_class_map.append(class_name)
        
        order_name_order_map[order] = 0
        order_order_name_map.append(order)
        
        family_order_map[family] = 0
        order_family_map.append(family) 
        
        genus_order_map[genus] = 0
        order_genus_map.append(genus) 
        
    # Now walk through the species
    # to get sorting order
    for species_name in species_list:
        
        kingdom,phylum,class_name,order,family,genus = species_taxonomy_map[species_name]
        
        if kingdom not in kingdom_order_map:
            kingdom_order_map[kingdom] = len(order_kingdom_map)
            order_kingdom_map.append(kingdom)
        
        if phylum not in phylum_order_map:
            phylum_order_map[phylum] = len(order_phylum_map)
            order_phylum_map.append(phylum)
        
        if class_name not in class_order_map:
            class_order_map[class_name] = len(order_class_map)
            order_class_map.append(class_name)
        
        if order not in order_name_order_map:
            order_name_order_map[order] = len(order_order_name_map)
            order_order_name_map.append(order)
        
        if family not in family_order_map:
            family_order_map[family] = len(order_family_map)
            order_family_map.append(family) 
        
        if genus not in genus_order_map:
            genus_order_map[genus] = len(order_genus_map)
            order_genus_map.append(genus)     
    
    
    # Now actually do the sort
    order_list = []
    for species_name in species_list:
        
        kingdom,phylum,class_name,order,family,genus = species_taxonomy_map[species_name]     
        order_list.append( (kingdom_order_map[kingdom], phylum_order_map[phylum], class_order_map[class_name], order_name_order_map[order], family_order_map[family], genus_order_map[genus]) )
    
    if len(second_sorting_attribute)==0:
        second_sorting_attribute = list(species_list)
    
    first_entry_list = [1-(species==first_entry) for species in species_list]
        
    # sort in descending order of sample size
    # Sort by num haploids    
    sorted_order_list, sorted_first_entry_list, sorted_second_sorting_attribute, sorted_species_list = zip(*sorted(zip(order_list, first_entry_list, second_sorting_attribute, species_list)))

    return sorted_species_list


# Returns a copy of the list of species sorted phylogenetically
def old_sort_phylogenetically(species):

    
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
    
    