import numpy
import parse_midas_data

def parse_sample_metadata_map():
    
    isolate_metadata_map = {}
    
    # load simulations
    file = open(parse_midas_data.scripts_directory+"isolates_genome_list.txt","r")
    file.readline() # 
    for line in file:
        items = line.strip().split("\t")
        subject_id = items[0] 
        sample_id = subject_id 
        accession_id=subject_id
        country = "isolate"
        continent = "isolate"
        order = 1
        
        isolate_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)
        
    file = open(parse_midas_data.scripts_directory+"mixture_labels.txt","r")
    file.readline() # header
    for line in file:
        items = line.strip().split("\t")
        subject_id = items[0] # this is one of two of the 90/10 mixtures
        sample_id = items[1] # This is the exact simulation
        accession_id=sample_id # same as sample
        country = "mixture"
        continent = "mixture"
        order = 1
        isolate_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)
        
    return isolate_metadata_map

    
def filter_sample_metadata_map(sample_metadata_map, field, field_value):
    
    if field=="country":
        field_idx = 3
    elif field=="continent":
        field_idx = 4
    elif field=="order":
        field_idx = 5
    else:
        return sample_metadata_map
    
    filtered_sample_metadata_map = {}
    for sample in sample_metadata_map:
        if sample_metadata_map[sample][field_idx]==field_value:
            filtered_sample_metadata_map[sample] = sample_metadata_map[sample]
            
    return sample_metadata_map
     


###############################################################################
#
# returns map from sample to (subject_id, temporal_order)
#
###############################################################################
def parse_sample_order_map(sample_metadata_map = {}): 
    
    import config
    
    if len(sample_metadata_map)==0:
        # Load it 
        sample_metadata_map = parse_sample_metadata_map()
    
    sample_order_map = {}
    for sample in sample_metadata_map:
        sample_order_map[sample] = (sample_metadata_map[0], sample_metadata_map[5])

    return sample_order_map



###############################################################################
#
# Loads metadata for HMP samples 
# Returns map from subject -> map of samples -> set of accession IDs
#
###############################################################################
def parse_subject_sample_map(sample_metadata_map = {}): 
    
    import config
    
    if len(sample_metadata_map)==0:
        # Load it 
        sample_metadata_map = parse_sample_metadata_map()

    
    subject_sample_map = {}
    for sample_id in sample_metadata_map:
        subject_id, dummy, accession_id, country, continent, order = sample_metadata_map[sample_id]
     
        if subject_id not in subject_sample_map:
            subject_sample_map[subject_id] = {}
            
        if sample_id not in subject_sample_map[subject_id]:
            subject_sample_map[subject_id][sample_id] = set()
            
        subject_sample_map[subject_id][sample_id].add(accession_id)
    
 
    return subject_sample_map 
   
#####
#
# Loads country metadata for samples
#
#####
def parse_sample_country_map(sample_metadata_map = {}): 
    
    import config
    
    if len(sample_metadata_map)==0:
        # Load it 
        sample_metadata_map = parse_sample_metadata_map()
        if sample_id not in sample_country_map:
            sample_country_map[sample_id] = country
    
    return sample_country_map


def list_of_isolates_and_mixtures():
    
    isolate_metadata_map = parse_isolate_metadata_map()
    
    isolates=[]
    mixtures=[]
    
    for sample_id in isolate_metadata_map:
        subject_id, dummy, accession_id, country, continent, order = sample_metadata_map[sample_id]
        if country=='isolate':
            isolates.append(sample_id)
        elif country=='mixture':
            mixtures.append(sample_id)
             
    return isolates, mixtures


###############################################################################
#
# Prunes sample list to only include samples from allowed countries
# Returns len(sampe_list) boolean array with element=False if sample was pruned  
#
###############################################################################
def calculate_country_samples(sample_country_map, sample_list=[], allowed_countries=set([])):

    if len(sample_list)==0:
        sample_list = list(sorted(sample_country_map.keys()))
        
    allowed_idxs = []
    for sample in sample_list:
        
        desired_sample = sample
            
        if (len(allowed_countries))==0 or (sample_country_map[desired_sample] in allowed_countries):
            allowed_idxs.append(True)
        else:
            allowed_idxs.append(False)
            
    allowed_idxs = numpy.array(allowed_idxs)
    return allowed_idxs


