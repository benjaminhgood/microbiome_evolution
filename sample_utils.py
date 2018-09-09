import numpy
from config import *

###############################################################################
#
# Methods for parsing sample metadata
#
###############################################################################

def parse_merged_sample_names(items):
    samples = []
    for item in items:
        sample = item.strip()
        if sample.endswith('c'):
            sample = sample[:-1]
        samples.append(sample)
        
    samples = numpy.array(samples)
    return samples


def calculate_sample_subject_map(subject_sample_map):
    
    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject
            
    return sample_subject_map
    

###############################################################################
#
# Creates a map of indexes from one list of samples (sample_list_from)
# to another list of samples (sample_list_to). The from list must be a 
# strict subset of the to list. 
#
###############################################################################
def calculate_sample_idx_map(sample_list_from, sample_list_to):
    
    sample_list_to = list(sample_list_to)
    sample_map = {}
    for i in xrange(0,len(sample_list_from)):
        sample_map[i] = sample_list_to.index(sample_list_from[i])
    
    return sample_map

def apply_sample_index_map_to_indices(sample_idx_map, idxs):
    new_idxs = (numpy.array([sample_idx_map[i] for i in idxs[0]]), numpy.array([sample_idx_map[i] for i in idxs[1]]))
    return new_idxs

def sample_name_lookup(sample_name, samples):
    
    for sample in samples:
    
        if sample.startswith(sample_name):
            return sample
            
    return ""

###############################################################################
#
# Prunes sample list to remove multiple timepoints from same subject
# Returns len(sampe_list) boolean array with element=False if sample was pruned  
#
###############################################################################
def calculate_unique_samples(subject_sample_map, sample_list=[]):

    if len(sample_list)==0:
        sample_list = list(sorted(flatten_samples(subject_sample_map).keys()))
    
    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject
    
    subject_idx_map = {}
        
    for i in xrange(0,len(sample_list)):
        sample = sample_list[i]
        if sample.endswith('c'):
            sample = sample[:-1]
        subject = sample_subject_map[sample]
        if not subject in subject_idx_map:
            subject_idx_map[subject] = i
            
    unique_idxs = numpy.zeros(len(sample_list),dtype=numpy.bool_)
    for i in subject_idx_map.values():
        unique_idxs[i]=True
    
    return unique_idxs
    
### 
#
# Returns a vector of true false values 
#
###
def calculate_samples_in_different_subjects(subject_sample_map, sample_list, focal_sample):

    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject

    in_different_subject = []
    for sample in sample_list:
        if sample_subject_map[sample] == sample_subject_map[focal_sample]:
            in_different_subject.append(False)
        else:
            in_different_subject.append(True)
                        
    return numpy.array(in_different_subject)
    
    
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
        
        if sample.endswith('c'):
            desired_sample = sample[:-1]
        else:
            desired_sample = sample
            
        if (len(allowed_countries))==0 or (sample_country_map[desired_sample] in allowed_countries):
            allowed_idxs.append(True)
        else:
            allowed_idxs.append(False)
            
    allowed_idxs = numpy.array(allowed_idxs)
    return allowed_idxs
    

###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with idx1 and idx2. All pairs are included 
# only once. 
#
###############################################################################
def calculate_subject_pairs(subject_sample_map, sample_list=[]):

    if len(sample_list)==0:
        sample_list = list(sorted(flatten_samples(subject_sample_map).keys()))
    
    new_sample_list = []
    for sample in sample_list:
        if sample.endswith('c'):
            new_sample_list.append(sample[:-1])
        else: 
            new_sample_list.append(sample)
    
    sample_list = new_sample_list
    
    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject
    
    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []
        
    for i in xrange(0,len(sample_list)):
        same_sample_idx_lower.append(i)
        same_sample_idx_upper.append(i)
        for j in xrange(0,i):
            if sample_subject_map[sample_list[i]]==sample_subject_map[sample_list[j]]:
                same_subject_idx_lower.append(i)
                same_subject_idx_upper.append(j)
            else: 
                diff_subject_idx_lower.append(i)
                diff_subject_idx_upper.append(j)
    
    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))
    
    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))
    
    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))
    
    return same_sample_idxs, same_subject_idxs, diff_subject_idxs

###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with idx1 and idx2. All pairs are included 
# only once. 
#
###############################################################################
def calculate_old_ordered_subject_pairs(sample_order_map, sample_list=[]):

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        for j in xrange(i,len(sample_list)):
            # loop over all pairs of samples
            
            if i==j:
                same_sample_idx_lower.append(i)
                same_sample_idx_upper.append(j)
            else:
                
                subject1, order1 = sample_order_map[sample_list[i]]
                subject2, order2 = sample_order_map[sample_list[j]]
                
                if subject1==subject2:
                    # same subject!
                    if order2-order1==1:
                        # consecutive samples
                        same_subject_idx_lower.append(i)
                        same_subject_idx_upper.append(j)
                    elif order1-order2==1:
                        # consecutive samples
                        same_subject_idx_lower.append(j)
                        same_subject_idx_upper.append(i)
                    else:
                        # do not add
                        pass
                    
                else:
                    # different subjects!
                    # Only take first one (to prevent multiple comparisons)
                    if order1==1 and order2==1:
                        diff_subject_idx_lower.append(i)
                        diff_subject_idx_upper.append(j)
        
    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))
    
    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))
    
    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))
    
    return same_sample_idxs, same_subject_idxs, diff_subject_idxs
    
###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with idx1 and idx2. All pairs are included 
# only once. 
#
###############################################################################
def calculate_ordered_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive'):

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    diff_subject_pair_map = {}
    same_subject_pair_map = {}
    # in this list, store lower_idx, upper_idx, interval_length

    # same sample pairs.. trivial
    same_sample_idx_lower = numpy.arange(0,len(sample_list))
    same_sample_idx_upper = numpy.arange(0,len(sample_list))

    # reconstruct "timeseries" for each subject
    subject_order_idx_map = {}
    for i in xrange(0,len(sample_list)):
        subject, order = sample_order_map[sample_list[i]]    
        
        if subject not in subject_order_idx_map:
            subject_order_idx_map[subject] = {}
        
        subject_order_idx_map[subject][order] = i
    
    # create index pairs within subjects
    for subject in subject_order_idx_map:
        
        sorted_orders = list(sorted(subject_order_idx_map[subject].keys()))
        
        # if at least two samples
        if len(sorted_orders)>1.5:
            
            if within_host_type=='longest':
                same_subject_idx_lower.append( subject_order_idx_map[subject][sorted_orders[0]] )
                same_subject_idx_upper.append( subject_order_idx_map[subject][sorted_orders[-1]] )
            elif within_host_type=='consecutive':
                for order_idx in xrange(1,len(sorted_orders)):
                    same_subject_idx_lower.append( subject_order_idx_map[subject][sorted_orders[order_idx-1]] )
                    same_subject_idx_upper.append( subject_order_idx_map[subject][sorted_orders[order_idx]] )
            elif within_host_type=='nonconsecutive':
                for order_idx_i in xrange(0,len(sorted_orders)):
                    for order_idx_j in xrange(order_idx_i+1,len(sorted_orders)):
                        same_subject_idx_lower.append( subject_order_idx_map[subject][sorted_orders[order_idx_i]] )
                        same_subject_idx_upper.append( subject_order_idx_map[subject][sorted_orders[order_idx_j]] )
               
    # now create index pairs in different subjects
    sorted_subjects = sorted(subject_order_idx_map.keys())
    
    for subject_i_idx in xrange(0,len(sorted_subjects)):
        subject_i = sorted_subjects[subject_i_idx]
        i = subject_order_idx_map[subject_i][min(subject_order_idx_map[subject_i].keys())]
        
        for subject_j_idx in xrange(subject_i_idx+1,len(sorted_subjects)):
            subject_j = sorted_subjects[subject_j_idx]
            j = subject_order_idx_map[subject_j][min(subject_order_idx_map[subject_j].keys())]
            
            diff_subject_idx_lower.append(i)
            diff_subject_idx_upper.append(j)
                 
    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))
    
    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))
    
    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))
    
    return same_sample_idxs, same_subject_idxs, diff_subject_idxs
        
    
###############################################################################
#
# For a given list of samples, calculates which belong to different subjects
# which belong to different timepoints in same subject, and which are the same
# timepoint.
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with idx1 and idx2. 
# The same pair can be included multiple times if there are multiple timepoints
#
###############################################################################
def calculate_nonconsecutive_ordered_subject_pairs(sample_order_map, sample_list=[]):

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        for j in xrange(i,len(sample_list)):
            # loop over all pairs of samples
            
            if i==j:
                same_sample_idx_lower.append(i)
                same_sample_idx_upper.append(j)
            else:
                
                subject1, order1 = sample_order_map[sample_list[i]]
                subject2, order2 = sample_order_map[sample_list[j]]
                
                if subject1==subject2:
                    # same subject!
                    if order2-order1>0.5:
                        # consecutive samples
                        same_subject_idx_lower.append(i)
                        same_subject_idx_upper.append(j)
                    elif order1-order2>0.5:
                        # consecutive samples
                        same_subject_idx_lower.append(j)
                        same_subject_idx_upper.append(i)
                    else:
                        # do not add
                        pass
                    
                else:
                    # different subjects!
                    # Only take first one (to prevent multiple comparisons)
                    if order1==1 and order2==1:
                        diff_subject_idx_lower.append(i)
                        diff_subject_idx_upper.append(j)
        
    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))
    
    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))
    
    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))
    
    return same_sample_idxs, same_subject_idxs, diff_subject_idxs
    

###############################################################################
# For a given list of samples, calculates which belong to different timepoints in same subject
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with idx1 and idx2. All pairs are included 
# only once. 
#
###############################################################################
def calculate_ordered_subject_triplets(sample_order_map, sample_list=[]):

    same_subject_idxs = []
    
    for i in xrange(0,len(sample_list)):
        
        subject1, order1 = sample_order_map[sample_list[i]]
        if order1 != 1:
            continue
                
        for j in xrange(0,len(sample_list)):
            
            subject2, order2 = sample_order_map[sample_list[j]]
            
            if subject2 != subject1:
                continue
                
            if order2 != 2:
                continue
                
            for k in xrange(0,len(sample_list)):
            
                subject3, order3 = sample_order_map[sample_list[k]]
            
                if subject3 != subject1:
                    continue
                    
                if order3 != 3:
                    continue
                    
                # if you get here, a triplet! 
                same_subject_idxs.append((i,j,k))    
            
    return same_subject_idxs



###############################################################################
#
# Calculates the subset that are sampled three times
#
###############################################################################
def calculate_triple_samples(sample_order_map, sample_list=[]):

    sample_idx_map = {}
    
    
    for i in xrange(0,len(sample_list)):
        subject, order = sample_order_map[sample_list[i]]
        
        if subject not in sample_idx_map:
            sample_idx_map[subject] = {}
            
        sample_idx_map[subject][order] = i
        
    triple_samples = []
    for subject in sample_idx_map.keys():
        if len(sample_idx_map[subject].keys()) > 2:
            triple_samples.append(numpy.array([sample_idx_map[subject][order] for order in sorted(sample_idx_map[subject].keys())]))
            
    return triple_samples


###############################################################################
#
# Returns matrix where rows are subjects and columns are hosts 
# where A_ih = 1 if sample i is in host h 
#
###############################################################################
def calculate_sample_subject_matrix(samples):

    sample_idx_map = {samples[i]:i for i in xrange(0,len(samples))}

    subject_sample_map = parse_subject_sample_map()
    subjects = subject_sample_map.keys()
    
    sample_subject_matrix = numpy.zeros((len(samples),len(subjects)),dtype=numpy.bool)
    
    for subject_idx in xrange(0,len(subjects)):
        for sample in subject_sample_map[subjects[subject_idx]]:
            if sample in sample_idx_map:
                sample_subject_matrix[sample_idx_map[sample], subject_idx] = True
    
    return sample_subject_matrix, subjects
    
    
    
    

    same_sample_idx_lower = []
    same_sample_idx_upper = []
    same_subject_idx_lower = []
    same_subject_idx_upper = []
    diff_subject_idx_lower = []
    diff_subject_idx_upper = []

    for i in xrange(0,len(sample_list)):
        for j in xrange(i,len(sample_list)):
            # loop over all pairs of samples
            
            if i==j:
                same_sample_idx_lower.append(i)
                same_sample_idx_upper.append(j)
            else:
                
                subject1, order1 = sample_order_map[sample_list[i]]
                subject2, order2 = sample_order_map[sample_list[j]]
                
                if subject1==subject2:
                    # same subject!
                    if order2-order1==1:
                        # consecutive samples
                        same_subject_idx_lower.append(i)
                        same_subject_idx_upper.append(j)
                    elif order1-order2==1:
                        # consecutive samples
                        same_subject_idx_lower.append(j)
                        same_subject_idx_upper.append(i)
                    else:
                        # do not add
                        pass
                    
                else:
                    # different subjects!
                    # Only take first one (to prevent multiple comparisons)
                    if order1==1 and order2==1:
                        diff_subject_idx_lower.append(i)
                        diff_subject_idx_upper.append(j)
        
    same_sample_idxs = (numpy.array(same_sample_idx_lower,dtype=numpy.int32), numpy.array(same_sample_idx_upper,dtype=numpy.int32))
    
    same_subject_idxs = (numpy.array(same_subject_idx_lower,dtype=numpy.int32), numpy.array(same_subject_idx_upper,dtype=numpy.int32))
    
    diff_subject_idxs = (numpy.array(diff_subject_idx_lower,dtype=numpy.int32), numpy.array(diff_subject_idx_upper,dtype=numpy.int32))
    
    return same_sample_idxs, same_subject_idxs, diff_subject_idxs

###############################################################################
#
# Returns a flat map of all the replicate sets for
# the samples in subject_sample_map, indexed by sample key        
#
###############################################################################
def flatten_samples(subject_sample_map):
    
    grouping_replicate_map = {}
    for subject in sorted(subject_sample_map.keys()):
        for sample in sorted(subject_sample_map[subject].keys()):
            grouping_replicate_map[sample] = subject_sample_map[subject][sample]
    
    return grouping_replicate_map


###############################################################################
#
# Returns a flat map of the merged replicate sets for each subject, 
# indexed by subject key 
#   
###############################################################################    
def flatten_subjects(subject_sample_map):
    
    grouping_replicate_map = {}
    for subject in sorted(subject_sample_map.keys()):
        merged_replicates = set()
        for sample in subject_sample_map[subject].keys():
            merged_replicates.update(subject_sample_map[subject][sample])
        grouping_replicate_map[subject] = merged_replicates
        
    return grouping_replicate_map


###############################################################################
#
# groupings = ordered list of nonoverlapping sets of sample names
# samples = ordered list of samples
#
# returns: list whose i-th element contains a numpy array of idxs
#          of the items in samples that are present in the ith grouping
#   
###############################################################################       
def calculate_grouping_idxs(groupings, samples):
    
    grouping_idxs = []
    for i in xrange(0,len(groupings)):
    
        idxs = []
        for j in xrange(0,len(samples)):
            if samples[j] in groupings[i]:
                idxs.append(j)
        idxs = numpy.array(idxs,dtype=numpy.int32)
        #print idxs
        grouping_idxs.append(idxs)
    
    return grouping_idxs
