import parse_midas_data

good_species_list = parse_midas_data.parse_good_species_list()

num_genomes=0
for species_name in good_species_list: 
        genome_ids=parse_midas_data.get_ref_genome_ids(species_name)
        for genome_id in genome_ids:
            print genome_id
            num_genomes +=1
