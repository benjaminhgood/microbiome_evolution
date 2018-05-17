import config
import gzip

def get_pangenome_map(species_name):
    
    gene_info_filename = '%span_genomes/%s/gene_info.txt.gz' % (config.midas_directory, species_name)
    file = gzip.open(gene_info_filename, 'r')
    file.readline() # header
    
    pangenome_map = {}
    
    for line in file:
        items = line.split("\t")
        gene_id = items[0].strip()
        genome_id = items[1].strip()
        centroid_99 = items[2].strip()
        centroid_95 = items[3].strip()
        
        if genome_id not in pangenome_map:
            pangenome_map[genome_id] = {}
        
        pangenome_map[genome_id][gene_id] = (centroid_99, centroid_95)
        
    file.close()
    return pangenome_map
    
def get_number_of_genomes(species_name):
    
    return len(get_pangenome_map(species_name))
    
if __name__=='__main__':
    
    import parse_midas_data
    good_species_list = parse_midas_data.parse_good_species_list()
    for species_name in good_species_list:
        print species_name
        print get_number_of_genomes(species_name)