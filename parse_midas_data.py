import numpy
import sys
import bz2
import os.path 
import stats_utils

data_directory = os.path.expanduser("~/ben_nandita_hmp_data/")
analysis_directory = os.path.expanduser("~/ben_nandita_hmp_analysis/")

default_directory_prefix =  data_directory

###############################################################################
#
# Returns: the string to append at the front of the filename for a given 
#          combination type.
#
# Valid combination types are None, "sample", "subject"
#
###############################################################################
def get_snp_prefix_from_combination_type(combination_type=None):
    if combination_type==None:
        snp_prefix = ""
    else:
        snp_prefix = combination_type + "_"
    
    return snp_prefix

###############################################################################
#
# Loads metadata for HMP samples 
# Returns map from subject -> map of samples -> set of accession IDs
#
###############################################################################
def parse_subject_sample_map(filename=os.path.expanduser("~/projectBenNandita/HMP_ids.txt")): # NOTE THAT THE PATH HAS CHANGED. Also, use os.path.expanduser to read what the absolute path of '~' is
    file = open(filename,"r")
    file.readline() # header
    
    
    subject_sample_map = {}
    
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
        
    return subject_sample_map 

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
        subject = sample_subject_map[sample_list[i]]
        if not subject in subject_idx_map:
            subject_idx_map[subject] = i
            
    unique_idxs = numpy.zeros(len(sample_list),dtype=numpy.bool_)
    for i in subject_idx_map.values():
        unique_idxs[i]=True
    
    return unique_idxs

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
    

def parse_gene_coverages(desired_species_name, combination_type=None, directory_prefix=default_directory_prefix):

    snp_prefix = get_snp_prefix_from_combination_type(combination_type)
    
    coverage_file = bz2.BZ2File("%ssnps/%s/%sgene_coverage.txt.bz2" % (default_directory_prefix, desired_species_name, snp_prefix))

    line = coverage_file.readline() # header
    items = line.split()
    samples = items[1:]
    gene_coverages = {}
    
    for line in coverage_file:
        items = line.split()
        gene_name = items[0]
        depths = numpy.array([float(item) for item in items[1:]])
        
        gene_coverages[gene_name] = depths
        
    return gene_coverages, samples

def parse_coverage_distribution(desired_species_name, combination_type=None, directory_prefix=default_directory_prefix):

    snp_prefix = get_snp_prefix_from_combination_type(combination_type)
    
    coverage_distribution_file = bz2.BZ2File("%ssnps/%s/%scoverage_distribution.txt.bz2" % (default_directory_prefix, desired_species_name, snp_prefix))

    line = coverage_distribution_file.readline() # header
    samples = []
    sample_coverage_histograms = []
    for line in coverage_distribution_file:
        items = line.split()
        sample_coverage_histogram = {}
        for item in items[1:]:
            subitems = item.split(",")
            sample_coverage_histogram[float(subitems[0])] = float(subitems[1])
        sample_coverage_histograms.append(sample_coverage_histogram)
        samples.append(items[0])
        
    return sample_coverage_histograms, samples
    

def parse_marker_gene_coverages(desired_species_name, combination_type=None, directory_prefix=default_directory_prefix):
    
    snp_prefix = get_snp_prefix_from_combination_type(combination_type)
    
    marker_file = bz2.BZ2File("%ssnps/%s/%smarker_coverage.txt.bz2" % (default_directory_prefix, desired_species_name, snp_prefix))
    
    line = marker_file.readline() # header
    samples = line.split()[1:]
    species = []
    species_coverage_matrix = []
    
    for line in marker_file:
        items = line.split()
        species_name = items[0]
        coverages = numpy.array([float(item) for item in items[1:]])
        species.append(species_name)
        species_coverage_matrix.append(coverages)
    
    marker_file.close()    
    
    species_coverage_matrix = numpy.array(species_coverage_matrix)
    return species_coverage_matrix, samples, species

    
def parse_global_marker_gene_coverages(combination_type=None, depth_threshold=10, directory_prefix=default_directory_prefix):

    if combination_type==None:
        coverage_file_prefix = ""
    else:
        coverage_file_prefix = combination_type + "_"
    

    file = bz2.BZ2File("%sspecies/%scoverage.txt.bz2" %  (directory_prefix, coverage_file_prefix),"r")
    line = file.readline() # header
    samples = line.split()[1:]
    species = []
    species_coverage_matrix = []
    for line in file:
        items = line.split()
        species_name = items[0]
        #print items
        coverages = numpy.array([float(item) for item in items[1:]])
        
        if coverages.sum() > depth_threshold:
            species.append(species_name)
            species_coverage_matrix.append(coverages)
    
    file.close()    
    species, species_coverage_matrix = zip(*sorted(zip(species, species_coverage_matrix), key=lambda pair: pair[1].sum(), reverse=True))
    
    species_coverage_matrix = numpy.array(species_coverage_matrix)
    return species_coverage_matrix, samples, species


def calculate_relative_depth_threshold_map(species_coverage_vector, samples, avg_depth_threshold=20, lower_factor=0.5, upper_factor=2):
    
    # returns map of sample name: coverage threshold
    # essentially filtering out samples whose marker depth coverage
    # does not exceed the average coverage threshold
    
    depth_threshold_map = {}
    for i in xrange(0,len(samples)):
        
        if species_coverage_vector[i]<avg_depth_threshold:    
            lower_depth_threshold=1000000001
            upper_depth_threshold=1000000001
        else:
            lower_depth_threshold = species_coverage_vector[i]*lower_factor
            upper_depth_threshold = species_coverage_vector[i]*upper_factor
    
        depth_threshold_map[samples[i]] = (lower_depth_threshold, upper_depth_threshold)
        
    return depth_threshold_map


def calculate_absolute_depth_threshold_map(species_coverage_vector, samples, avg_depth_threshold=20, site_depth_threshold=15):
    
    # returns map of sample name: coverage threshold
    # essentially filtering out samples whose marker depth coverage
    # does not exceed the average coverage threshold
    
    depth_threshold_map = {}
    for i in xrange(0,len(samples)):
        
        if species_coverage_vector[i]<avg_depth_threshold:    
            lower_depth_threshold=1000000001
        else:
            lower_depth_threshold=site_depth_threshold
    
        upper_depth_threshold = 1000000001
        depth_threshold_map[samples[i]] = (lower_depth_threshold, upper_depth_threshold)
        
    return depth_threshold_map
    
    
def combine_replicates(species_name, grouping_replicate_map, output_prefix, directory_prefix=default_directory_prefix,debug=False):
    # writes new files after combining technical replicates according to groupings in the grouping -> {replicates} map
      
    ref_freq_file = bz2.BZ2File("%ssnps/%s/snps_ref_freq.txt.bz2" % (default_directory_prefix, species_name),"r")
    depth_file = bz2.BZ2File("%ssnps/%s/snps_depth.txt.bz2" % (default_directory_prefix, species_name),"r")
    
    combined_ref_freq_file = bz2.BZ2File("%ssnps/%s/%s_snps_ref_freq.txt.bz2" % (default_directory_prefix, species_name, output_prefix),"w")
    combined_depth_file = bz2.BZ2File("%ssnps/%s/%s_snps_depth.txt.bz2" % (default_directory_prefix, species_name, output_prefix),"w")
    
    # get header lines from each file
    depth_line = depth_file.readline()
    ref_freq_line = ref_freq_file.readline()
    
    # get list of samples
    depth_items = depth_line.split()
    samples = depth_items[1:]
    
    # get list of samples
    ref_freq_items = ref_freq_line.split()

    # calculate new groupings and their associated sample idxs
    grouping_ids = []
    replicate_groupings = []
    for grouping_id in sorted(grouping_replicate_map.keys()):
        
        grouping_ids.append(grouping_id)
        replicate_groupings.append(grouping_replicate_map[grouping_id])
        
    grouping_idxs = calculate_grouping_idxs(replicate_groupings, samples)

    # write header lines for both output files
    combined_depth_file.write("\t".join([depth_items[0]]+grouping_ids))
    combined_ref_freq_file.write("\t".join([ref_freq_items[0]]+grouping_ids))
    
    
    num_sites_processed = 0
    while True:
            
        # load next lines
        depth_line = depth_file.readline()
        ref_freq_line = ref_freq_file.readline()
        
        # quit if file has ended
        if depth_line=="":
            break
        
        depth_items = depth_line.split()
        ref_freq_items = ref_freq_line.split()

        # now parse allele count info
        depths = numpy.array([float(item) for item in depth_items[1:]])
        ref_freqs = numpy.array([float(item) for item in ref_freq_items[1:]]) 
        refs = ref_freqs*depths    
        
        combined_depths = numpy.array([depths[grouping_idxs[i]].sum() for i in xrange(0,len(grouping_idxs))])
        combined_refs = numpy.array([refs[grouping_idxs[i]].sum() for i in xrange(0,len(grouping_idxs))])
        combined_ref_freqs = (combined_refs+(combined_depths==0))*1.0/(combined_depths+(combined_depths==0))
        
        # write output lines for both output files
        combined_depth_file.write("\n")
        combined_depth_file.write("\t".join([depth_items[0]]+["%g" % d for d in combined_depths]))
        combined_ref_freq_file.write("\n")
        combined_ref_freq_file.write("\t".join([ref_freq_items[0]]+["%g" % f for f in combined_ref_freqs]))
    
        num_sites_processed+=1
        if num_sites_processed%10000==0:
            sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
            if debug:
                break
    
    ref_freq_file.close()
    depth_file.close()
    combined_ref_freq_file.close()
    combined_depth_file.close()
    
    # Done!
    # no return value

def combine_species_marker_gene_replicates(species_name, grouping_replicate_map, output_prefix, directory_prefix=default_directory_prefix):
    # writes new files after combining technical replicates according to groupings in the grouping -> {replicates} map
    
    coverage_file = bz2.BZ2File("%ssnps/%s/marker_coverage.txt.bz2" % (default_directory_prefix, species_name),"r")
    
    combined_coverage_file = bz2.BZ2File("%ssnps/%s/%s_marker_coverage.txt.bz2" % (default_directory_prefix, species_name, output_prefix),"w")
    
    # get header line
    line = coverage_file.readline()
    
    # get list of samples
    items = line.split()
    samples = items[1:]
    
    # calculate new groupings and their associated sample idxs
    grouping_ids = []
    replicate_groupings = []
    for grouping_id in sorted(grouping_replicate_map.keys()):
        
        grouping_ids.append(grouping_id)
        replicate_groupings.append(grouping_replicate_map[grouping_id])
        
    grouping_idxs = calculate_grouping_idxs(replicate_groupings, samples)

    # write header lines for output file
    combined_coverage_file.write("\t".join([items[0]]+grouping_ids))
    
    num_sites_processed = 0
    while True:
            
        # load next lines
        line = coverage_file.readline()
        
        # quit if file has ended
        if line=="":
            break
        
        items = line.split()
        
        # now parse allele count info
        depths = numpy.array([float(item) for item in items[1:]])
        
        combined_depths = numpy.array([depths[grouping_idxs[i]].sum() for i in xrange(0,len(grouping_idxs))])
        
        # write output lines for both output files
        combined_coverage_file.write("\n")
        combined_coverage_file.write("\t".join([items[0]]+["%g" % d for d in combined_depths]))
        
        
    coverage_file.close()
    combined_coverage_file.close()
    
    # Done!
    # no return value

    
def combine_marker_gene_replicates(grouping_replicate_map, output_prefix, directory_prefix=default_directory_prefix):
    # writes new files after combining technical replicates according to groupings in the grouping -> {replicates} map
      
    coverage_file = bz2.BZ2File("%sspecies/coverage.txt.bz2" % (default_directory_prefix),"r")
    combined_coverage_file = bz2.BZ2File("%sspecies/%s_coverage.txt.bz2" % (default_directory_prefix, output_prefix),"w")
    
    # get header line
    line = coverage_file.readline()
    
    # get list of samples
    items = line.split()
    samples = items[1:]
    
    # calculate new groupings and their associated sample idxs
    grouping_ids = []
    replicate_groupings = []
    for grouping_id in sorted(grouping_replicate_map.keys()):
        
        grouping_ids.append(grouping_id)
        replicate_groupings.append(grouping_replicate_map[grouping_id])
        
    grouping_idxs = calculate_grouping_idxs(replicate_groupings, samples)

    # write header lines for output file
    combined_coverage_file.write("\t".join([items[0]]+grouping_ids))
    
    num_sites_processed = 0
    while True:
            
        # load next lines
        line = coverage_file.readline()
        
        # quit if file has ended
        if line=="":
            break
        
        items = line.split()
        
        # now parse allele count info
        depths = numpy.array([float(item) for item in items[1:]])
        
        combined_depths = numpy.array([depths[grouping_idxs[i]].sum() for i in xrange(0,len(grouping_idxs))])
        
        # write output lines for both output files
        combined_coverage_file.write("\n")
        combined_coverage_file.write("\t".join([items[0]]+["%g" % d for d in combined_depths]))
        
        
    coverage_file.close()
    combined_coverage_file.close()
    
    # Done!
    # no return value


###############################################################################
#
# Loads list of SNPs and counts of target sites from annotated SNPs file
#
# returns (lots of things, see below)
#
###############################################################################
def parse_snps(species_name, combination_type=None, debug=False):
    
    if combination_type==None:
        snp_prefix = ""
    else:
        snp_prefix = combination_type + "_"
    
    
    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File("%ssnps/%s/%sannotated_snps.txt.bz2" % (default_directory_prefix, species_name, snp_prefix),"r")
    
    line = snp_file.readline() # header
    items = line.split()
    samples = items[1:]
    # Only going to look at 1D and 4D sites
    allowed_variant_types = set(['1D','4D'])
    
    allele_counts_map = {}
    
    # map from gene_name -> var_type -> (list of locations, matrix of allele counts)
    passed_sites_map = {}
    # map from gene_name -> var_type -> (location, sample x sample matrix of whether both samples can be called at that site)
    
    
    num_sites_processed = 0
    for line in snp_file:
        
        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        variant_type = info_items[3]
        pvalue = float(info_items[4])
        
        # make sure it is either a 1D or 4D site
        # (shouldn't be needed anymore)
        if not variant_type in allowed_variant_types:
            continue
        
        # Load alt and depth counts
        alts = []
        depths = []
        for item in items[1:]:
            subitems = item.split(",")
            alts.append(float(subitems[0]))
            depths.append(float(subitems[1]))
        alts = numpy.array(alts)
        depths = numpy.array(depths)
        refs = depths-alts

        passed_sites = (depths>0)*1.0
        if gene_name not in passed_sites_map:
            passed_sites_map[gene_name] = {v: {'location': (chromosome,location), 'sites': numpy.zeros((len(samples), len(samples)))} for v in allowed_variant_types}
            
            allele_counts_map[gene_name] = {v: {'locations':[], 'alleles':[]} for v in allowed_variant_types}
        
        
        passed_sites_map[gene_name][variant_type]['sites'] += passed_sites[:,None]*passed_sites[None,:]
        
        # zero out non-passed sites
        # (shouldn't be needed anymore)    
        refs = refs*passed_sites
        alts = alts*passed_sites
        depths = depths*passed_sites
        
        # calculate whether SNP has passed
        
        # Criteria used in Schloissnig et al (Nature, 2013)
        #total_alts = alts.sum()
        #total_depths = depths.sum()
        #pooled_freq = total_alts/((total_depths+(total_depths==0))
        #snp_passed = (freq>0.01) and (total_alts>=4) and ((total_depths-total_alts)>=4)
        
        # new version
        #alt_threshold = numpy.ceil(depths*0.05)+0.5 #at least one read above 5%.
        #alts = alts*((alts>alt_threshold))
        #snp_passed = (alts.sum()>0) and (pvalue<0.05)
        
        # "pop gen" version
        alt_lower_threshold = numpy.ceil(depths*0.05)+0.5 #at least one read above 5%.
        alts = alts*((alts>alt_lower_threshold))
        
        alt_upper_threshold = alt_lower_threshold
        snp_passed = ((alts>alt_upper_threshold).sum()>0) and (pvalue<0.05)
        
        # consensus approximation
        #alt_upper_threshold = depths*0.95
        #snp_passed = ((alts>alt_upper_threshold).sum()>0)
        
        
        #print alts.sum()>0, pvalue, (pvalue < 5e-02), snp_passed
        
        allele_counts = numpy.transpose(numpy.array([alts,refs]))
        
        allele_counts_map[gene_name][variant_type]['locations'].append((chromosome, location))
        allele_counts_map[gene_name][variant_type]['alleles'].append(allele_counts)
        
        num_sites_processed+=1
        if num_sites_processed%10000==0:
            sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
            if debug:
                break
    
    snp_file.close()

    for gene_name in passed_sites_map.keys():
        for variant_type in passed_sites_map[gene_name].keys():
            
            allele_counts_map[gene_name][variant_type]['alleles'] = numpy.array(allele_counts_map[gene_name][variant_type]['alleles'])

    return samples, allele_counts_map, passed_sites_map
    

  

###############################################################################
#
# Reads midas output and prints to stdout in a format 
# suitable for further downstream processing
#
###############################################################################
def pipe_snps(species_name, combination_type=None, avg_depth_threshold=20, directory_prefix=default_directory_prefix,debug=False):
    
    
    # Load marker gene coverage data
    species_coverage_matrix, sample_list, species_list = parse_marker_gene_coverages(species_name, combination_type, directory_prefix=directory_prefix)
    marker_gene_coverages = species_coverage_matrix[species_list.index(species_name),:]
    # Load genomic coverage distributions
    sample_coverage_histograms, samples = parse_coverage_distribution(species_name, combination_type=combination_type)

    median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])

    #depth_vector =  marker_gene_coverages # we don't use this anymore
    depth_vector = median_coverages
    
    depth_threshold_map = calculate_relative_depth_threshold_map(depth_vector, sample_list, avg_depth_threshold=avg_depth_threshold) #, site_depth_threshold)
    
    if combination_type==None:
        snp_prefix = ""
    else:
        snp_prefix = combination_type + "_"
    
    # Open MIDAS output files
    ref_freq_file = bz2.BZ2File("%ssnps/%s/%ssnps_ref_freq.txt.bz2" % (default_directory_prefix, species_name, snp_prefix),"r")
    depth_file = bz2.BZ2File("%ssnps/%s/%ssnps_depth.txt.bz2" % (default_directory_prefix, species_name, snp_prefix),"r")
    alt_allele_file = bz2.BZ2File("%ssnps/%s/snps_alt_allele.txt.bz2" % (default_directory_prefix,species_name),"r")
    info_file = bz2.BZ2File("%ssnps/%s/snps_info.txt.bz2" % (default_directory_prefix, species_name),"r")
    marker_file = bz2.BZ2File("%ssnps/%s/%smarker_coverage.txt.bz2" % (default_directory_prefix, species_name, snp_prefix))
    
    # get header lines from each file
    depth_line = depth_file.readline()
    ref_freq_line = ref_freq_file.readline()
    alt_line = alt_allele_file.readline()
    info_line = info_file.readline()
    marker_line = marker_file.readline()
    
    # get list of samples
    depth_items = depth_line.split()
    samples = numpy.array(depth_items[1:])
    
    # get marker gene coverages
    
    # create depth threshold vector from depth threshold map
    lower_depth_threshold_vector = []
    upper_depth_threshold_vector = []
    for sample in samples:
        lower_depth_threshold_vector.append(depth_threshold_map[sample][0])
        upper_depth_threshold_vector.append(depth_threshold_map[sample][1])
        
    lower_depth_threshold_vector = numpy.array(lower_depth_threshold_vector)
    upper_depth_threshold_vector = numpy.array(upper_depth_threshold_vector)
    
    # Figure out which samples passed our avg_depth_threshold
    passed_samples = (lower_depth_threshold_vector<1e09)
    total_passed_samples = passed_samples.sum()
    
    # Let's focus on those from now on
    samples = list(samples[passed_samples])
    lower_depth_threshold_vector = lower_depth_threshold_vector[passed_samples]
    upper_depth_threshold_vector = upper_depth_threshold_vector[passed_samples]
    
    #print lower_depth_threshold_vector
    
    # print header
    print_str = "\t".join(["site_id"]+samples)
    print print_str
    
    # Only going to look at 1D and 4D sites
    allowed_variant_types = set(['1D','4D'])
    
    allele_counts_syn = [] # alt and reference allele counts at 4D synonymous sites with snps
    locations_syn = [] # genomic location of 4D synonymous sites with snps
    genes_syn = [] # gene name of 4D synonymous sites with snps
    passed_sites_syn = numpy.zeros(len(samples))*1.0
    
    allele_counts_non = [] # alt and reference allele counts at 1D nonsynonymous sites with snps
    locations_non = [] # genomic location of 1D nonsynonymous sites
    genes_non = [] # gene name of 1D nonsynonymous sites with snps
    passed_sites_non = numpy.zeros_like(passed_sites_syn)
    
    num_sites_processed = 0
    while True:
            
        # load next lines
        depth_line = depth_file.readline()
        ref_freq_line = ref_freq_file.readline()
        alt_line = alt_allele_file.readline()
        info_line = info_file.readline()
        
        # quit if file has ended
        if depth_line=="":
            break
        
        # parse site info
        info_items = info_line.split()
        variant_type = info_items[5]
        
        # make sure it is either a 1D or 4D site
        if not variant_type in allowed_variant_types:
            continue
    
        # continue parsing site info
        gene_name = info_items[6]
        site_id_items = info_items[0].split("|")
        contig = site_id_items[0]
        location = site_id_items[1]
        new_site_id_str = "|".join([contig, location, gene_name, variant_type])
        
        
    
        # now parse allele count info
        depths = numpy.array([float(item) for item in depth_line.split()[1:]])[passed_samples]
        ref_freqs = numpy.array([float(item) for item in ref_freq_line.split()[1:]])[passed_samples]
        refs = numpy.round(ref_freqs*depths)   
        alts = depths-refs
        
        passed_sites = (depths>=lower_depth_threshold_vector)*1.0
        passed_sites *= (depths<=upper_depth_threshold_vector)
        
        #print passed_sites.sum(), total_passed_samples, passed_sites.sum()/total_passed_samples
        
        # make sure the site is prevalent in enough samples to count as "core"
        if (passed_sites).sum()*1.0/total_passed_samples < 0.5:
            continue
            #passed_sites *= 0
            
        refs = refs*passed_sites
        alts = alts*passed_sites
        depths = depths*passed_sites
        
        total_alts = alts.sum()
        total_refs = depths.sum()
        total_depths = total_alts+total_refs
        
        
        # polarize SNP based on consensus in entire dataset
        if total_alts>total_refs:
            alts,refs = refs,alts
            total_alts, total_refs = total_refs, total_alts
        
        # print string
        read_strs = ["%g,%g" % (A,A+R) for A,R in zip(alts, refs)]
        print_str = "\t".join([new_site_id_str]+read_strs)
        
        print print_str
        #print total_alts
        
        num_sites_processed+=1
        if num_sites_processed%10000==0:
            sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
            if debug:
                break
    
    ref_freq_file.close()
    depth_file.close()
    alt_allele_file.close()
    info_file.close()
    
    # returns nothing
    

if __name__=='__main__':

    pass
    
