import parse_midas_data
import numpy


###############################################################################
#
# Loads metadata for HMP samples 
# Returns map from subject -> map of samples -> set of accession IDs
#
###############################################################################
def parse_sample_order_map(): 

    sample_order_map = {}
    
    # First load HMP metadata
    file = open(parse_midas_data.scripts_directory+"HMP_ids_order.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        order = long(items[5].strip())
        
        sample_order_map[sample_id] = (subject_id, order)
        
    file.close()

    # Then load Kuleshov data 
    file = open(parse_midas_data.scripts_directory+"kuleshov_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        
        sample_order_map[sample_id] = (subject_id, 1)
    file.close()
    
    # Then load Qin data
    file = open(parse_midas_data.scripts_directory+"qin_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        sample_id = accession_id # Nandita used accession id as MIDAS header for this dataset
        country = items[3].strip()
        continent = items[4].strip()
        
        sample_order_map[sample_id] = (subject_id, 1)
        
    file.close()
    
    return sample_order_map

###############################################################################
#
# Loads metadata for HMP samples 
# Returns map from subject -> map of samples -> set of accession IDs
#
###############################################################################
def parse_subject_sample_map(): 

    subject_sample_map = {}
    
    # First load HMP metadata
    file = open(parse_midas_data.scripts_directory+"HMP_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        
        if subject_id not in subject_sample_map:
            subject_sample_map[subject_id] = {}
            
        if sample_id not in subject_sample_map[subject_id]:
            subject_sample_map[subject_id][sample_id] = set()
            
        subject_sample_map[subject_id][sample_id].add(accession_id)
    file.close()
    
    # Then load Kuleshov data 
    file = open(parse_midas_data.scripts_directory+"kuleshov_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        
        if subject_id not in subject_sample_map:
            subject_sample_map[subject_id] = {}
            
        if sample_id not in subject_sample_map[subject_id]:
            subject_sample_map[subject_id][sample_id] = set()
            
        subject_sample_map[subject_id][sample_id].add(accession_id)
    file.close()
    
    # Then load Qin data
    file = open(parse_midas_data.scripts_directory+"qin_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        sample_id = accession_id # Nandita used accession id as MIDAS header for this dataset
        country = items[3].strip()
        continent = items[4].strip()
        
        if subject_id not in subject_sample_map:
            subject_sample_map[subject_id] = {}
            
        if sample_id not in subject_sample_map[subject_id]:
            subject_sample_map[subject_id][sample_id] = set()
            
        subject_sample_map[subject_id][sample_id].add(accession_id)
    file.close()
      
    
    # Repeat for other data
    # Nothing else so far
     
    return subject_sample_map 
   

#####
#
# Loads country metadata for samples
#
#####
def parse_sample_country_map(): 

    sample_country_map = {}
    
    # First load HMP metadata
    file = open(parse_midas_data.scripts_directory+"HMP_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
            
        if sample_id not in sample_country_map:
            sample_country_map[sample_id] = country
            
    file.close()
    
    # Then load Kuleshov data 
    file = open(parse_midas_data.scripts_directory+"kuleshov_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        
        if sample_id not in sample_country_map:
            sample_country_map[sample_id] = country
    
    file.close()
    
    # Then load Qin data
    file = open(parse_midas_data.scripts_directory+"qin_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        sample_id = accession_id # Nandita used accession id as MIDAS header for this dataset
        country = items[3].strip()
        continent = items[4].strip()
        
        if sample_id not in sample_country_map:
            sample_country_map[sample_id] = country
    
    file.close()
      
    
    # Repeat for other data
    # Nothing else so far
     
    return sample_country_map

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


def parse_subject_phenotype_map():
    
    column_idx = 5

    subject_phenotype_map = {}
    
    # First load HMP metadata
    file = open(parse_midas_data.scripts_directory+"HMP_phenotypes.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        city = items[column_idx].strip()
        
        if city=='city1':
            phenotype=0
        elif city=='city2':
            phenotype=1
        else:
            phenotype=-1
            
        subject_phenotype_map[subject_id] = phenotype
            
    file.close()
    
    
    return subject_phenotype_map
    
def parse_sample_phenotype_map():

    subject_phenotype_map = parse_subject_phenotype_map()
    subject_sample_map = parse_subject_sample_map()
    sample_subject_map = parse_midas_data.calculate_sample_subject_map(subject_sample_map)
    sample_phenotype_map = {}
    for sample in sample_subject_map.keys():
        subject = sample_subject_map[sample]
        if subject in subject_phenotype_map:
            sample_phenotype_map[sample] = subject_phenotype_map[subject]
        else:
            sample_phenotype_map[sample] = -1    
    
    return sample_phenotype_map

