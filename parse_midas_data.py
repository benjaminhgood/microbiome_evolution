import numpy
import sys
import bz2
import gzip
import os.path 
import stats_utils
from math import floor, ceil
import gene_diversity_utils

import config
import sample_utils

###############################################################################
#
# Set up default source and output directories
#
###############################################################################

data_directory = config.data_directory
analysis_directory = config.analysis_directory
scripts_directory = config.scripts_directory
patric_directory = config.patric_directory
midas_directory = config.midas_directory

# We use this one to debug because it was the first one we looked at
debug_species_name = config.debug_species_name

from sample_utils import *

###############################################################################
#
# Methods for parsing species metadata
#
###############################################################################

#############
#
# Returns a list of all species that MIDAS called SNPS for
#
#############
def parse_species_list():
    
    species_names = []
    
    file = open(data_directory+"snps/species_snps.txt","r")
    for line in file:
        species_names.append(line.strip())
    file.close()
    
    return species_names

#############
#
# Returns a list of all species that MIDAS called SNPS for
# sorted in order of decreasing total sequencing depth
#
#############
def parse_depth_sorted_species_list():
    species_coverage_matrix, samples, species = parse_global_marker_gene_coverages()
    return species


#############
#
# Returns a list of all species that MIDAS called SNPS for
# that passed a certain depth / prevalence requirement,
# again sorted in order of decreasing total sequencing depth
#
#############
def parse_good_species_list(min_marker_coverage=config.good_species_min_coverage, min_prevalence=config.good_species_min_prevalence):
    good_species_list = []
    
    species_coverage_matrix, samples, species = parse_global_marker_gene_coverages()
    for i in xrange(0,len(species)):
        
        species_coverages = species_coverage_matrix[i,:]
        if (species_coverages>=min_marker_coverage).sum() >= min_prevalence:
            good_species_list.append(species[i])
    
    return good_species_list

    
    
    
###############################################################################
#
# Methods for parsing coverage of different species across samples
#
###############################################################################

###############################################################################
#
# Loads marker gene coverages produced by MIDAS
# for all species in which SNPs were called 
#
# Returns: species-by-sample matrix of marker gene coverages,
#          with species sorted in descending order of total coverage;
#          ordered list of sample ids; ordered list of species names;
#
###############################################################################
def parse_global_marker_gene_coverages():

    desired_species_names = set(parse_species_list())

    file = bz2.BZ2File("%sspecies/coverage.txt.bz2" %  (data_directory),"r")
    line = file.readline() # header
    samples = line.split()[1:]
    species = []
    species_coverage_matrix = []
    for line in file:
        items = line.split()
        species_name = items[0]
        #print items
        coverages = numpy.array([float(item) for item in items[1:]])
        
        if species_name in desired_species_names:
            species.append(species_name)
            species_coverage_matrix.append(coverages)
    
    file.close()    
    species, species_coverage_matrix = zip(*sorted(zip(species, species_coverage_matrix), key=lambda pair: pair[1].sum(), reverse=True))
    
    species_coverage_matrix = numpy.array(species_coverage_matrix)
    return species_coverage_matrix, samples, species

def parse_species_marker_gene_coverages(desired_species_name):
    
    species_coverage_matrix, samples, species = parse_species_marker_gene_coverages()
    
    for species_idx in xrange(0,len(species)):
        species_name = species[species_idx]
        if species_name==desired_species_name:
            return species_coverage_matrix[species_idx,:], samples

    return None
    
def parse_sample_coverage_map(desired_species_name):
    import stats_utils
    
    # Load genomic coverage distributions
    sample_coverage_histograms, samples = parse_coverage_distribution(desired_species_name)
    median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    return sample_coverage_map

def parse_gene_coverages(desired_species_name):

    coverage_file = bz2.BZ2File("%ssnps/%s/gene_coverage.txt.bz2" % (data_directory, desired_species_name))

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

def parse_marker_gene_coverage_distribution(desired_species_name):
    
    coverage_distribution_file = bz2.BZ2File("%ssnps/%s/marker_coverage_distribution.txt.bz2" % (data_directory, desired_species_name))

    line = coverage_distribution_file.readline() # header
    marker_gene_coverages = {}
    for line in coverage_distribution_file:
        items = line.split("\t")
        subitems = items[0].split(",")
        sample = subitems[0].strip()
        gene_name = subitems[1].strip()
        
        if sample not in marker_gene_coverages:
            marker_gene_coverages[sample] = {}
            
        locations = []
        coverages = []
        for item in items[1:]:
            subitems = item.split(",")
            locations.append( long(subitems[0]) )
            coverages.append( float(subitems[1]) )
        
        locations = numpy.array(locations)
        coverages = numpy.array(coverages)
        
        marker_gene_coverages[sample][gene_name] = (locations,coverages)
        
    return marker_gene_coverages

def parse_coverage_distribution(desired_species_name,prevalence_filter=True,remove_c=True):

    if prevalence_filter:
        full_str = ""
    else:
        full_str = "full_"
        

    coverage_distribution_file = bz2.BZ2File("%ssnps/%s/%scoverage_distribution.txt.bz2" % (data_directory, desired_species_name, full_str))

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
    
    if remove_c == True:
        samples = parse_merged_sample_names(samples)    
    return sample_coverage_histograms, samples
    
## 
# 
# Loads species-specific marker gene coverage
#
##
def parse_marker_gene_coverages(desired_species_name):
    
    marker_file = bz2.BZ2File("%ssnps/%s/marker_coverage.txt.bz2" % (data_directory, desired_species_name))
    
    line = marker_file.readline() # header
    samples = parse_merged_sample_names(line.split()[1:])
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


################
#
# Methods for determining whether samples or sites pass certain depth requirements
#
################
  

def calculate_relative_depth_threshold_map(sample_coverage_histograms, samples, min_nonzero_median_coverage=5, lower_factor=0.5, upper_factor=2):
    
    # returns map of sample name: coverage threshold
    # essentially filtering out samples whose marker depth coverage
    # does not exceed the average coverage threshold
    
    depth_threshold_map = {}
    for i in xrange(0,len(samples)):
        
        # Check if coverage distribution meets certain requirements
        is_bad_coverage_distribution = False
        
        # First check if passes median coverage requirement
        nonzero_median_coverage = stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histograms[i])
        if round(nonzero_median_coverage) < min_nonzero_median_coverage:
            is_bad_coverage_distribution=True
    
        # Passed median coverage requirement
        # Now check whether a significant number of sites fall between lower and upper factor. 
        lower_depth_threshold = floor(nonzero_median_coverage*lower_factor)-0.5 # why is 0.5 being added/subtracted? NRG
        upper_depth_threshold = ceil(nonzero_median_coverage*upper_factor)+0.5
    
        depths, depth_CDF = stats_utils.calculate_CDF_from_histogram(sample_coverage_histograms[i])
        # remove zeros
        if depths[0]<0.5:
            depth_CDF -= depth_CDF[0]
            depth_CDF /= depth_CDF[-1]
        
        fraction_in_good_range = depth_CDF[(depths>lower_depth_threshold)*(depths<upper_depth_threshold)].sum()
    
        if fraction_in_good_range < 0.6: #where does 0.6 come from? NRG
            is_bad_coverage_distribution=True
            
        if is_bad_coverage_distribution:
            lower_depth_threshold = 1000000001
            upper_depth_threshold = 1000000001
        
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
  

###############################################################################
#
# Reads midas output and prints to stdout in a format 
# suitable for further downstream processing
#
# In the process, filters sites that fail to meet the depth requirements
#
###############################################################################
def pipe_snps(species_name, min_nonzero_median_coverage=config.pipe_snps_min_nonzero_median_coverage, lower_factor=config.pipe_snps_lower_depth_factor, upper_factor=config.pipe_snps_upper_depth_factor, min_samples=config.pipe_snps_min_samples, debug=False):


# lower_factor = 0.3 is the default to be consistent with MIDAS gene presence criterion
# upper factor = 3 is the default for (logarithmic) symmetry 
# min_samples=4 is the default because then the site is guaranteed to be present in 
# at least 2 independent people. 
# NRG: Why is a site guaranteed to be in at least 2 independent people?
# BG: In our cohort, the maximum number of samples per person is 3. 
#     If there are 4 samples, then they must be spread across at least 2 people. 
    
    # Load genomic coverage distributions
    sample_coverage_histograms, sample_list = parse_coverage_distribution(species_name, remove_c=False)
    # depth threshold map returns the lower and upper depth values that are 0.3*median and 3*median depth in the data. 
    depth_threshold_map = calculate_relative_depth_threshold_map(sample_coverage_histograms, sample_list, min_nonzero_median_coverage, lower_factor, upper_factor)
    
   
    # Open MIDAS output files
    ref_freq_file = bz2.BZ2File("%ssnps/%s/snps_ref_freq.txt.bz2" % (data_directory, species_name),"r")
    depth_file = bz2.BZ2File("%ssnps/%s/snps_depth.txt.bz2" % (data_directory, species_name),"r")
    alt_allele_file = bz2.BZ2File("%ssnps/%s/snps_alt_allele.txt.bz2" % (data_directory, species_name),"r")
    info_file = bz2.BZ2File("%ssnps/%s/snps_info.txt.bz2" % (data_directory, species_name),"r")
    marker_file = bz2.BZ2File("%ssnps/%s/marker_coverage.txt.bz2" % (data_directory, species_name))
    
    # get header lines from each file
    depth_line = depth_file.readline()
    ref_freq_line = ref_freq_file.readline()
    alt_line = alt_allele_file.readline()
    info_line = info_file.readline()
    marker_line = marker_file.readline()
    
    # get list of samples
    depth_items = depth_line.split()
    samples = numpy.array(depth_items[1:])
    
    # BHG (06/24/17) removed this so that all "raw" data have "c"s. 
    # All functions that write something keep them. 
    # All loading functions strip them. 
    #samples= parse_merged_sample_names(samples) # NRG (06/06/17): I added this so that the keys in dictionary are compatible. 

    # samples
    prevalence_threshold = min([min_samples*1.0/len(samples), 0.5])
    
    # create depth threshold vector from depth threshold map
    lower_depth_threshold_vector = []
    upper_depth_threshold_vector = []
    for sample in samples:
        lower_depth_threshold_vector.append(depth_threshold_map[sample][0])
        upper_depth_threshold_vector.append(depth_threshold_map[sample][1])
        
    lower_depth_threshold_vector = numpy.array(lower_depth_threshold_vector)
    upper_depth_threshold_vector = numpy.array(upper_depth_threshold_vector)
    
    # Figure out which samples passed our avg_depth_threshold
    passed_samples = (lower_depth_threshold_vector<1e09) #1e09 comes from the calculate_relative_depth_threshold_map definition above, which is a code for a bad sample. A bad sample has median depth less than 5 or greater than 0.6 fraction of the genome is outside the acceptable range of good depths. 
    total_passed_samples = passed_samples.sum()
    
    # Let's focus on those from now on
    samples = list(samples[passed_samples])
    lower_depth_threshold_vector = lower_depth_threshold_vector[passed_samples]
    upper_depth_threshold_vector = upper_depth_threshold_vector[passed_samples]
    
    #print lower_depth_threshold_vector
    
    # print header
    print_str = "\t".join(["site_id"]+samples)
    print print_str
    
    # Only going to look at 1D, 2D, 3D, and 4D sites
    # (we will restrict to 1D and 4D downstream)
    allowed_variant_types = set(['1D','2D','3D','4D'])
    
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
        info_items = info_line.split("\t")
        variant_type = info_items[5]
        
        # make sure it is either a 1D or 4D site
        if not variant_type in allowed_variant_types:
            continue
    
        # continue parsing site info
        gene_name = info_items[6]
        site_id_items = info_items[0].split("|")
        # NRG: added this if condition to deal with extra 'accn' in db swap. 
        if site_id_items[0]=='accn':
            contig = site_id_items[1]
            location = site_id_items[2]
        else:
            contig = site_id_items[0] 
            location = site_id_items[1]
        
    
        # now parse allele count info
        depths = numpy.array([float(item) for item in depth_line.split()[1:]])
        ref_freqs = numpy.array([float(item) for item in ref_freq_line.split()[1:]])
        
        depths = depths[passed_samples]
        ref_freqs = ref_freqs[passed_samples]
        
        #if (ref_freqs==1.0).all():
        #    sys.stderr.write("Non-polymorphic site!\n")
        
        refs = numpy.round(ref_freqs*depths)   
        alts = depths-refs
        
        passed_sites = (depths>=lower_depth_threshold_vector)*1.0
        passed_sites *= (depths<=upper_depth_threshold_vector)
        
        #print passed_sites.sum(), total_passed_samples, passed_sites.sum()/total_passed_samples
        
        # make sure the site is prevalent in enough samples to count as "core"
        if (passed_sites).sum()*1.0/total_passed_samples < prevalence_threshold:
            continue
            #passed_sites *= 0
            
        refs = refs*passed_sites
        alts = alts*passed_sites
        depths = depths*passed_sites
        
        total_alts = alts.sum()
        total_refs = depths.sum()
        total_depths = total_alts+total_refs
        
        # BG: 05/18: moving the polarization part to another part of the pipeline
        # so that we can use HMP polarization with other datasets. 
        # at the moment, still saving polarization state.
        
        polarization = "R"    
        new_site_id_str = "|".join([contig, location, gene_name, variant_type, polarization])
        
        
        # print string
        read_strs = ["%g,%g" % (A,A+R) for A,R in zip(alts, refs)]
        print_str = "\t".join([new_site_id_str]+read_strs)
        
        print print_str
        #print total_alts
        
        num_sites_processed+=1
        if num_sites_processed%10000==0:
            #sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
            if debug:
                break
    
    ref_freq_file.close()
    depth_file.close()
    alt_allele_file.close()
    info_file.close()
    
    # returns nothing

###############################################################################
#
# Loads list of SNPs and counts of target sites from annotated SNPs file
#
# returns (lots of things, see below)
#
###############################################################################
def parse_snps(species_name, debug=False, allowed_samples=[], allowed_genes=[], allowed_variant_types=['1D','2D','3D','4D'], initial_line_number=0, chunk_size=1000000000):
    
    import calculate_snp_prevalences
    # Load population freqs (for polarization purposes)    
    population_freqs = calculate_snp_prevalences.parse_population_freqs(species_name)
   
    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps.txt.bz2" % (data_directory, species_name),"r")
    
    line = snp_file.readline() # header
    items = line.split()[1:]    
    samples = parse_merged_sample_names(items)
    
    if len(allowed_samples)==0:
        allowed_sample_set = set(samples)
    else:
        allowed_sample_set = (set(allowed_samples) & set(samples))
        
    allowed_genes = set(allowed_genes)
    allowed_variant_types = set(allowed_variant_types)
    
    # This is a hack because there were some mistaken repeats in an old data file
    # should be able to remove later
    seen_samples = set()
    desired_sample_idxs = []
    for sample in allowed_samples:
        if (sample in allowed_sample_set) and (sample not in seen_samples):
            desired_sample_idxs.append( numpy.nonzero(samples==sample)[0][0] )
        else:
            pass
        
        seen_samples.add(sample)
        
    desired_sample_idxs = numpy.array(desired_sample_idxs)    
        
    desired_samples = samples[desired_sample_idxs]
    
    #print len(samples), len(desired_sample_idxs), len(allowed_samples), len(desired_samples), len(allowed_sample_set)

    
    # map from gene_name -> var_type -> (list of locations, matrix of allele counts)
    allele_counts_map = {}
    # map from gene_name -> var_type -> (location, sample x sample matrix of whether both samples can be called at that site)
    passed_sites_map = {}
    
    num_sites_processed = 0
    line_number = -1
    final_line_number = -1
    previous_gene_name = ""
    gene_name = ""
    for line in snp_file:
        
        line_number += 1
        previous_gene_name = gene_name
                
        if line_number < initial_line_number:
            continue
        
        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        variant_type = info_items[3]
        
        if len(info_items) > 5: # for backwards compatability
            polarization = info_items[4]
            pvalue = float(info_items[5])
        else: 
            polarization="R" # not correct, but avoids a crash
            pvalue = float(info_items[4])
            
        
        if num_sites_processed >= chunk_size and gene_name!=previous_gene_name:
            # We are done for now!
            final_line_number = line_number
            break
        
                
        if not variant_type in allowed_variant_types:
            continue
            
        if len(allowed_genes)>0 and (not gene_name in allowed_genes):
            continue
        
        # Load alt and depth counts
        # Load alt and depth counts
        alts = []
        depths = []
        
        for idx in desired_sample_idxs:    
            item = items[1+idx]
            subitems = item.split(",")
            alts.append(float(subitems[0]))
            depths.append(float(subitems[1]))
        alts = numpy.array(alts)
        depths = numpy.array(depths)
        
        
        # polarize
        if (chromosome, location) in population_freqs:
            population_freq = population_freqs[(chromosome, location)]
        else:
            population_freq = 0
        
        # polarize SFS according to population freq
        if population_freq>0.5:
            alts = depths-alts
            polarization = 'A'
     
        passed_sites = (depths>0)*1.0
        if gene_name not in passed_sites_map:
            passed_sites_map[gene_name] = {v: {'location': (chromosome,location), 'sites': numpy.zeros((len(desired_samples), len(desired_samples)))} for v in allowed_variant_types}
            
            allele_counts_map[gene_name] = {v: {'locations':[], 'alleles':[]} for v in allowed_variant_types}
        
        passed_sites_map[gene_name][variant_type]['sites'] += passed_sites[:,None]*passed_sites[None,:]
        
        # zero out non-passed sites
        # (shouldn't be needed anymore)    
        alts = alts*passed_sites
        depths = depths*passed_sites
        
        # calculate whether SNP has passed
        alt_threshold = numpy.ceil(depths*config.parse_snps_min_freq)+0.5 #at least one read above 5%.
        snp_passed = ((alts>alt_threshold).sum()>0) and (pvalue<0.05)
        
        # Criteria used in Schloissnig et al (Nature, 2013)
        #total_alts = alts.sum()
        #total_depths = depths.sum()
        #pooled_freq = total_alts/((total_depths+(total_depths==0))
        #snp_passed = (freq>0.01) and (total_alts>=4) and ((total_depths-total_alts)>=4)
        
        if snp_passed:
            allele_counts = numpy.transpose(numpy.array([alts,depths-alts]))
        
            allele_counts_map[gene_name][variant_type]['locations'].append((chromosome, location))
            allele_counts_map[gene_name][variant_type]['alleles'].append(allele_counts)
        
            num_sites_processed+=1
        
            if num_sites_processed>0 and num_sites_processed%1000==0:
                sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
                if debug:
                    break
         
    snp_file.close()

    #print line_number, final_line_number, num_sites_processed

    for gene_name in passed_sites_map.keys():
        for variant_type in passed_sites_map[gene_name].keys():
            
            allele_counts_map[gene_name][variant_type]['alleles'] = numpy.array(allele_counts_map[gene_name][variant_type]['alleles'])

    return desired_samples, allele_counts_map, passed_sites_map, final_line_number

###############################################################################
#
# Calculates within-sample synonymous pi directly from annotated_snps.txt.bz2. 
# Ugly hack (since it does not encourage code re-use and puts pop-gen logic
# in the basic parsing scripts) but we need it so that we can call parse_snps 
# on subsets of samples later on to improve performance
#
# returns vector of samples, vector of pi_s (raw counts), vector of opportunities
#
###############################################################################
def parse_within_sample_pi_new(species_name, allowed_genes=set([]), allowed_variant_types=set(['4D']), debug=False):
    
    samples, sfs_map = parse_within_sample_sfs(species_name, allowed_variant_types)
    
    total_pi = []
    total_opportunities = []
    for sample in samples: 
        p,n = diversity_utils.calculate_pi_from_sfs_map(sfs_map[sample])
        total_pi.append(p)
        total_opportunities.append(n)
        
    total_pi = numpy.array(total_pi)*1.0
    total_opportunities = numpy.array(total_opportunities)*1.0
    
    return samples, total_pi, total_opportunities
    


###############################################################################
#
# Calculates within-sample sfs directly from annotated_snps.txt.bz2. 
# Ugly hack (since it does not encourage code re-use and puts pop-gen logic
# in the basic parsing scripts) but we need it so that we can call parse_snps 
# on subsets of haploid samples later on to improve performance
#
# returns vector of samples, vector of sfs maps
#
# (
#
###############################################################################
def parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])):
    
    # First write (filtered) genome-wide coverage distribution
    sfs_file = bz2.BZ2File("%ssnps/%s/within_sample_sfs.txt.bz2" % (data_directory, species_name),"r")
    sfs_file.readline() # header
    
    sfs_map = {}
    samples = []
    for line in sfs_file:
        items = line.split("\t")
        sample = parse_merged_sample_names([items[0].strip()])[0]
        
        variant_type = items[1].strip()
        sfs_items = items[2:]
        
        
        if variant_type not in allowed_variant_types:
            continue
        
        if sample not in sfs_map:
            sfs_map[sample] = {}
            samples.append(sample)
            
        
            
        for sfs_item in sfs_items:
            subitems = sfs_item.split(",")
            D = long(subitems[0])
            A = long(subitems[1])
            n = long(subitems[2])
            reverse_n = float(subitems[3])
            
            if D<0.5:
                continue
    
            if (A,D) not in sfs_map[sample]:
                sfs_map[sample][(D,A)] = [0,0.0]
                
            sfs_map[sample][(D,A)][0] += n
            sfs_map[sample][(D,A)][1] += reverse_n
            
            
    return numpy.array(samples), sfs_map

def pangenome_data_exists(species_name):
   gene_reads_filename =  "%sgenes/%s/genes_reads.txt.bz2" % (data_directory, species_name)
   return os.path.isfile(gene_reads_filename)

###############################################################################
#
# Loads MIDAS's pangenome coverage data for a given species 
#
# returns (lots of things, see below)
#
###############################################################################
def parse_pangenome_data(species_name, allowed_samples = [], allowed_genes=[], convert_centroid_names=True, disallowed_genes=[]):
    
    if not pangenome_data_exists(species_name):
        return [], [], [], [], [], []
        
    # Open post-processed MIDAS output
    # Raw read counts
    gene_reads_file =  bz2.BZ2File("%sgenes/%s/genes_reads.txt.bz2" % (data_directory, species_name),"r")
    # Depth (read counts / length?)
    gene_depth_file =  bz2.BZ2File("%sgenes/%s/genes_depth.txt.bz2" % (data_directory, species_name),"r")
    # Presence absence calls
    gene_presabs_file =  bz2.BZ2File("%sgenes/%s/genes_presabs.txt.bz2" % (data_directory, species_name),"r")
    
    # First read through gene_summary_file to get marker gene coverages
    # Gene summary file
    gene_summary_file = file("%sgenes/%s/genes_summary.txt" % (data_directory, species_name),"r")
    marker_coverage_map = {}
    gene_summary_file.readline() # header
    marker_coverage_samples = []
    marker_coverages = []
    for summary_line in gene_summary_file:
        items = summary_line.split()
        sample = items[0].strip()
        marker_coverage = float(items[5])
        marker_coverage_samples.append(sample)
        marker_coverages.append(marker_coverage)

    gene_summary_file.close()

    marker_coverage_samples = parse_merged_sample_names(marker_coverage_samples)

    marker_coverage_map = {sample: marker_coverage for sample,marker_coverage in zip(marker_coverage_samples, marker_coverages)}
    
    # Now read through remaining files
    reads_line = gene_reads_file.readline() # header
    depth_line = gene_depth_file.readline() # header
    presabs_line = gene_presabs_file.readline() # header
    items = presabs_line.split()
    
    samples = parse_merged_sample_names(items[1:])
    
    # ordered vector of marker coverages (guaranteed to be in same order as samples)
    marker_coverages = numpy.array([marker_coverage_map[sample] for sample in samples])
    
    if len(allowed_samples)==0:
        allowed_samples = set(samples)
    else:
        allowed_samples = (set(allowed_samples) & set(samples))
        
    desired_sample_idxs = numpy.array([sample in allowed_samples for sample in samples])
    desired_samples = samples[desired_sample_idxs]
    
    marker_coverages = marker_coverages[desired_sample_idxs]
    gene_presence_matrix = []
    gene_reads_matrix = []
    gene_depth_matrix = []
    gene_names = []
    
    num_genes_processed = 0
    
    reads_line = gene_reads_file.readline() # header
    depth_line = gene_depth_file.readline() # header
    presabs_line = gene_presabs_file.readline() # header
    
    while reads_line!="":
        
        # Loop through!
        
        items = presabs_line.split()
        gene_name = items[0]
        gene_presences = numpy.array([float(item) for item in items[1:]])[desired_sample_idxs]
        
        if True: #gene_presences.sum() > 0.5:
        
            gene_reads = numpy.array([float(item) for item in reads_line.split()[1:]])[desired_sample_idxs]
            gene_depths = numpy.array([float(item) for item in depth_line.split()[1:]])[desired_sample_idxs]
            
            # Note to self: not uniform across samples!
            #gene_lengths = gene_reads/(gene_depths+(gene_reads<0.5))
            #print gene_lengths
            
            # gene is present in at least one individual! 
            gene_presence_matrix.append(gene_presences)
            gene_depth_matrix.append(gene_depths)
            gene_reads_matrix.append(gene_reads)
            gene_names.append(gene_name)    
        
        num_genes_processed+=1
        
        reads_line = gene_reads_file.readline() # header
        depth_line = gene_depth_file.readline() # header
        presabs_line = gene_presabs_file.readline() # header
    
        
    gene_reads_file.close()
    gene_depth_file.close()
    gene_presabs_file.close()
    gene_presence_matrix = numpy.array(gene_presence_matrix)
    gene_depth_matrix = numpy.array(gene_depth_matrix)
    gene_reads_matrix = numpy.array(gene_reads_matrix)

    if convert_centroid_names:
        new_gene_names = []
        centroid_gene_map = load_centroid_gene_map(species_name)
        for gene_name in gene_names:
            new_gene_names.append(centroid_gene_map[gene_name])
    else:
        new_gene_names=gene_names
    
    
    new_gene_names = numpy.array(new_gene_names)
        
    # Now weed out disallowed genes if provided
    disallowed_genes=set(disallowed_genes)
    allowed_gene_idxs = []
    for gene_idx in xrange(0,len(new_gene_names)):
        
        if new_gene_names[gene_idx] in disallowed_genes:
            # don't include
            pass
        else:
            allowed_gene_idxs.append(gene_idx)
    allowed_gene_idxs = numpy.array(allowed_gene_idxs)
    
    new_gene_names = new_gene_names[allowed_gene_idxs]
    gene_presence_matrix = gene_presence_matrix[allowed_gene_idxs,:]
    gene_depth_matrix = gene_depth_matrix[allowed_gene_idxs,:]
    gene_reads_matrix = gene_reads_matrix[allowed_gene_idxs,:]
    
    return desired_samples, new_gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix

###############################################################################
#
# A set of (lightweight) methods for loading gene_ids
# that satisfy various criteria (e.g. in reference genome
# part of metaphlan2 set, etc.
#
###############################################################################


####
#
# The gene_ids in the pangenome list are the centroids of gene clusters.
# Sometimes the gene in the reference genome is not chosen as the centroid.
# This function creates a map between pangenome_centroids and genes in 
# reference genome (if it exists)
#
###
def load_centroid_gene_map(desired_species_name):
    
    # First load reference genes
    reference_genes = load_reference_genes(desired_species_name)
    
    gene_info_file = gzip.open("%span_genomes/%s/gene_info.txt.gz" % (midas_directory, desired_species_name), 'r')
    
    gene_info_file.readline() # header
    
    centroid_gene_map = {}
    
    for line in gene_info_file:
        
        items = line.split("\t") 
        gene_id = items[0].strip()
        centroid_id = items[3].strip()
        
        if centroid_id not in centroid_gene_map:
            centroid_gene_map[centroid_id] = centroid_id
            
        if (gene_id in reference_genes) and (centroid_id not in reference_genes):
            centroid_gene_map[centroid_id] = gene_id
            
        
    gene_info_file.close()
    
    return centroid_gene_map


####
#
# Returns a map from gene name to (reference corrected centroids)
#
###
def load_gene_centroid_map(desired_species_name):
    
    gene_info_file = gzip.open("%span_genomes/%s/gene_info.txt.gz" % (midas_directory, desired_species_name), 'r')
    
    gene_info_file.readline() # header
    
    gene_centroid_map = {}
    
    for line in gene_info_file:
        
        items = line.split("\t") 
        gene_id = items[0].strip()
        centroid_id = items[3].strip()
        
        gene_centroid_map[gene_id] = centroid_id    
        
    gene_info_file.close()
    
    centroid_new_centroid_map = load_centroid_gene_map(desired_species_name)
    
    for gene_id in gene_centroid_map.keys():
        new_centroid = centroid_new_centroid_map[gene_centroid_map[gene_id]]
        gene_centroid_map[gene_id] = new_centroid
        
    return gene_centroid_map


####
#
# In this definition we return as a key the centroid and as value a list of ALL genes that are within the 95% clster.
#
###
def load_complete_centroid_gene_map(desired_species_name):
    
    gene_info_file = gzip.open("%span_genomes/%s/gene_info.txt.gz" % (midas_directory, desired_species_name), 'r')
    
    gene_info_file.readline() # header
    
    complete_centroid_gene_map = {}
    
    for line in gene_info_file:
        
        items = line.split("\t") 
        gene_id = items[0].strip()
        centroid_id = items[3].strip()
        
        if centroid_id not in complete_centroid_gene_map:
            complete_centroid_gene_map[centroid_id] = [gene_id]
        else:
            complete_centroid_gene_map[centroid_id].append(gene_id)
            
        
    gene_info_file.close()
    
    return complete_centroid_gene_map



###############################################################################
#
# Loads the MIDAS's pangenome (after clustering at X% identity)
# 
# returns a complete set of genes (irrespective of prevalence)
#
###############################################################################
def load_pangenome_genes(species_name):

    gene_names=[]
    # Open post-processed MIDAS output
    # Presence absence calls
    gene_presabs_file =  bz2.BZ2File("%sgenes/%s/genes_presabs.txt.bz2" % (data_directory, species_name),"r")
    presabs_line = gene_presabs_file.readline() # header
    
    while presabs_line!="":
        
        # Loop through!        
        items = presabs_line.split()
        gene_name = items[0]
        if gene_name !='gene_id':
            gene_names.append(gene_name)            
        presabs_line = gene_presabs_file.readline() # header

    centroid_gene_map = load_centroid_gene_map(species_name)

    new_species_names = []
    for gene_name in gene_names:
        new_species_names.append(centroid_gene_map[gene_name])
        
    return set(gene_names), set(new_species_names)


    
###############################################################################
#
# Loads list of genes in the reference genome used by MIDAS for a given species
#
###############################################################################
def load_reference_genes(desired_species_name):

    
    features_file = gzip.open("%srep_genomes/%s/genome.features.gz" % (midas_directory, desired_species_name), 'r')
    
    features_file.readline() # header
    reference_genes = []
    for line in features_file:
        items = line.split()
        gene_name = items[0].strip()
        reference_genes.append(gene_name)
    features_file.close()    
    
    return set(reference_genes)
    
################################################################################
#
# Loads metaphlan2 genes (that are present in the reference genome)
# returns a list of metaphlan2 genes
#
################################################################################
def load_metaphlan2_genes(desired_species_name):
    gene_file = open("%smetaphlan2_genes/%s_metaphlan2_genes_mapped.txt" % (data_directory, desired_species_name), 'r')
    
    reference_genes = load_reference_genes(desired_species_name)
    
    metaphlan2_genes=[]
    for line in gene_file:
        gene_name = line.strip()
        
        if gene_name in reference_genes:
            metaphlan2_genes.append(gene_name)
        else:
            pass
            #print gene_name

    gene_file.close()    
    
    return set(metaphlan2_genes)
    
###############################################################################
#
# Loads list of MIDAS marker genes for a given species
# (by default, load only those that are present in reference genome)
#
###############################################################################
def load_marker_genes(desired_species_name, require_in_reference_genome=True):   

    # Chosen markers (From table S7 of Nayfach et al (Genome Res 2016)
    marker_ids = set(['B000032', 'B000039', 'B000041', 'B000062', 'B000063', 'B000065', 'B000071', 'B000079', 'B000080', 'B000081', 'B000082', 'B000086', 'B000096', 'B000103', 'B000114'])

    marker_gene_map = {marker_id : [] for marker_id in marker_ids}

    reference_genes = set(load_reference_genes(desired_species_name))
    
    marker_gene_file = open("%smarker_genes/phyeco.map" % (midas_directory), 'r')
    
    marker_gene_file.readline() # header
    
    marker_genes = []
    #print "Markers absent from the reference genome"
    for line in marker_gene_file:
        
        items = line.split("\t")
        gene_name = items[0].strip()
        species_name = items[3].strip()
        marker_id = items[4].strip()
        
        if desired_species_name == species_name:
            if marker_id in marker_ids:
                if (not require_in_reference_genome) or (gene_name in reference_genes):
                    marker_genes.append(gene_name)
                    marker_gene_map[marker_id].append(gene_name)
                    
                else:
                    pass
                    #print gene_name, marker_id
            else:
                pass
    
    #print marker_gene_map
                
    #print 'Two markers that are present in the reference genome:'
    #print marker_genes[0]
    #print marker_genes[1]
 
    return set(marker_genes)
 
 
def parse_pangenome_species():
    pangenome_species = []
    
    species_list = parse_depth_sorted_species_list()
    
    for species_name in species_list:
    
        gene_directory = '%sgenes/%s' % (data_directory, species_name)
        if os.path.isdir(gene_directory):
            pangenome_species.append(species_name)
    pangenome_species.append('new_species')       
    return pangenome_species
     
###############################################################################
#
# Loads a subset of "core" genes. 
# *Deprecated: Use core_gene_utils.parse_core_genes instead
#
###############################################################################   
def load_core_genes(desired_species_name, min_copynum=0.3, min_prevalence=0.9, min_marker_coverage=20, unique_individuals=True):

    import core_gene_utils
    return core_gene_utils.parse_core_genes(desired_species_name)
    
    # Old
    #return load_core_genes_from_pangenome(desired_species_name, min_copynum, min_prevalence, min_marker_coverage, unique_individuals)

###############################################################################
#
# Loads a subset of "core" genes using copynum information in the genes/ folder 
#
###############################################################################   
def load_core_genes_from_pangenome(desired_species_name, min_copynum=0.3, min_prevalence=0.9, min_marker_coverage=20, unique_individuals=True):

    # Load subject and sample metadata
    subject_sample_map = parse_subject_sample_map()
    
    # Load reference genes
    reference_genes = load_reference_genes(desired_species_name)
    
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_pangenome_data(desired_species_name)
    
    gene_names = numpy.array(gene_names)
   
    reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])

    if unique_individuals:
        sample_idxs = (calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_marker_coverage)
    else:
        sample_idxs = marker_coverages>=min_marker_coverage
    
    prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,sample_idxs], marker_coverages[sample_idxs], min_copynum)

    core_gene_idxs = reference_gene_idxs*(prevalences>=min_prevalence)

    return set(gene_names[core_gene_idxs])
    

########################################################################################
#
# Loads time data for HMP samples 
# Returns map from subject_id -> visno -> [[sample_id, study_day_1], [sample_id, study_day_2], etc]
#
#######################################################################################
def parse_subject_sample_time_map(filename=os.path.expanduser("~/ben_nandita_hmp_data/HMP_ids_time.txt")): 
    file = open(filename,"r")
    file.readline() # header
    
    
    subject_sample_time_map = {}
    
    for line in file:
        items = line.split("\t")
        subject_id= items[0].strip()
        sample_id = items[1].strip()
        visno     = int(items[5].strip())
        study_day = int(items[6].strip())

        if subject_id not in subject_sample_time_map:
            subject_sample_time_map[subject_id] = {}
        if visno not in subject_sample_time_map[subject_id].keys():
            subject_sample_time_map[subject_id][visno]=[]
        
        found=False
        for i in range(0, len(subject_sample_time_map[subject_id][visno])):
            if sample_id == subject_sample_time_map[subject_id][visno][i][0]:
                found =True 
        if found==False:
            subject_sample_time_map[subject_id][visno].append([sample_id,study_day])
        
        
        
    return subject_sample_time_map 


###################################################################################
# reformat the subject_sample_time_map so that the samples that have been combined (because of visno replicates) are now in the dictionary instead of old samples that were merged into one.
##################################################################

def collapse_visno_reps_subject_sample_time_map(subject_sample_time_map):

    new_subject_sample_time_map={}
    for subject in subject_sample_time_map:
        new_subject_sample_time_map[subject]={}
        for visno in subject_sample_time_map[subject].keys(): # loop over samples
            new_subject_sample_time_map[subject][visno]=[]
            if len(subject_sample_time_map[subject][visno]) >1:
                sample_name=subject_sample_time_map[subject][visno][0][0] + 'c'
                day=subject_sample_time_map[subject][visno][0][1]
            else:
                sample_name=subject_sample_time_map[subject][visno][0][0]
                day =subject_sample_time_map[subject][visno][0][1]

            new_subject_sample_time_map[subject][visno].append([sample_name,day])

    return new_subject_sample_time_map
########################################################################################
#
# Prunes time data for HMP samples 
# If more than one sample is present for a visno, return the one with the most coverage. 
#
#######################################################################################
def prune_subject_sample_time_map(subject_sample_time_map_all_samples,sample_coverage_map):

    subject_sample_time_map={}

    for subject in subject_sample_time_map_all_samples.keys(): # loop over subjects (hosts)
        for visno in subject_sample_time_map_all_samples[subject].keys(): # loop over samples
            keep_coverage=0
            keep_sample=''
            keep_day=0
            if len(subject_sample_time_map_all_samples[subject][visno]) >1: # find the sample with highest cov
                
                for i in range(0,len(subject_sample_time_map_all_samples[subject][visno])): 
                    sample = subject_sample_time_map_all_samples[subject][visno][i][0]
                    day=subject_sample_time_map_all_samples[subject][visno][i][1]
                    if sample in sample_coverage_map.keys():
                        coverage=sample_coverage_map[sample]
                        if coverage>keep_coverage:
                            keep_coverage=coverage
                            keep_sample=sample
                            keep_day=day
            else:
                keep_sample = subject_sample_time_map_all_samples[subject][visno][0][0]
                keep_day=subject_sample_time_map_all_samples[subject][visno][0][1]
                
            if keep_sample !='':
                if subject not in subject_sample_time_map.keys():
                    subject_sample_time_map[subject]={}
                subject_sample_time_map[subject][visno]=[[keep_sample,keep_day]]
                        
    return subject_sample_time_map 


########################################################################################
#
# Returns index pairs for time points corresponding to the same subject_id.
# pairs corresponding to the first visno (stored in index1) and a subsequent visno (either visno 2 and/or 3) are returned. Visno pairs 2 and 3 are not returned. 
# Also returns the corresponding visnos and days. 
#
#######################################################################################

def calculate_time_pairs(subject_sample_time_map, samples):
    index1=[]
    index2=[]
    visno=[]
    day=[]

    for subject_id in subject_sample_time_map.keys():
        visnos=subject_sample_time_map[subject_id].keys() #visit numbers
        if (len(visnos) > 1) and (1 in visnos):           
            if (subject_sample_time_map[subject_id][1][0][0] in samples): #check if first visit in samples 
                #iterate through visit numbers. Append the index, day, and visnos to their lists
                for i in visnos:        
                    if (subject_sample_time_map[subject_id][i][0][0] in samples) and (i !=1):
                        index1.append(samples.tolist().index(subject_sample_time_map[subject_id][1][0][0]))
                        index2.append(samples.tolist().index(subject_sample_time_map[subject_id][i][0][0]))
                        visno.append(i)
                        day.append(subject_sample_time_map[subject_id][i][0][1])
        
    time_pair_idxs = (numpy.array(index1,dtype=numpy.int32), numpy.array(index2,dtype=numpy.int32))

    return time_pair_idxs, visno, day



########################################################################################
#
# Returns index pairs for time points corresponding to the same subject_id.
# pairs corresponding to the first visno (visno 1 or 2) (stored in index1) and a subsequent visno (either visno 2 and/or 3) are returned.  
# Also returns the corresponding visnos and days between the visnos. 
#
#######################################################################################

def calculate_all_time_pairs(subject_sample_time_map, samples):
    index1=[]
    index2=[]
    visno1=[]
    visno2=[]
    day=[]

    for subject_id in subject_sample_time_map.keys():
        visnos=subject_sample_time_map[subject_id].keys() #visit numbers
        if (len(visnos) > 1):            
            #iterate through visit numbers. Append the index, day, and visnos to their lists
            for i in range(0, len(visnos)):        
                for j in range(i+1, len(visnos)):
                    if (subject_sample_time_map[subject_id][visnos[i]][0][0] in samples and subject_sample_time_map[subject_id][visnos[j]][0][0] in samples):
                        if visnos[i] < visnos[j]:
                            first=visnos[i]
                            second=visnos[j]
                        else:
                            first=visnos[j]
                            second=visnos[i]
                        index1.append(samples.tolist().index(subject_sample_time_map[subject_id][first][0][0]))
                        index2.append(samples.tolist().index(subject_sample_time_map[subject_id][second][0][0]))
                        visno1.append(first)
                        visno2.append(second)
                        day.append(subject_sample_time_map[subject_id][second][0][1]-subject_sample_time_map[subject_id][first][0][1])
        
    time_pair_idxs = (numpy.array(index1,dtype=numpy.int32), numpy.array(index2,dtype=numpy.int32))

    return time_pair_idxs, visno1, visno2, day




########################################################################################
#
# Returns index pairs for time points corresponding to the same subject_id.
# Only one time pair per subject is returned (so, either 1 vs2, or 1 vs 3 is returned. If both exist, 1 vs 2 is returned).  Visno pairs 2 and 3 are not returned. 
# Also returns the corresponding visnos and days. 
#
#######################################################################################

def calculate_unique_time_pairs(subject_sample_time_map, samples):
    index1=[]
    index2=[]
    visno=[]
    day=[]

    for subject_id in subject_sample_time_map.keys():
        visnos=subject_sample_time_map[subject_id].keys() #visit numbers
        if (len(visnos) > 1) and (1 in visnos):           
            if (subject_sample_time_map[subject_id][1][0][0] in samples): #check if first visit in samples 
                #iterate through visit numbers. Append the index, day, and visnos to their lists
                unique_pair_found=False
                for i in [2,3]:
                    if i in subject_sample_time_map[subject_id].keys() and subject_sample_time_map[subject_id][i][0][0] in samples and unique_pair_found==False:
                        index1.append(samples.tolist().index(subject_sample_time_map[subject_id][1][0][0])) 
                        index2.append(samples.tolist().index(subject_sample_time_map[subject_id][i][0][0]))
                        visno.append(i)
                        day.append(subject_sample_time_map[subject_id][i][0][1])
                        unique_pair_found=True
                
    time_pair_idxs = (numpy.array(index1,dtype=numpy.int32), numpy.array(index2,dtype=numpy.int32))

    return time_pair_idxs, visno, day


########################################################################################  
# 
# Find the representative genome id in PATRIC for a given species 
#
#########################################################################################

def representative_genome_id(desired_species_name):
    species_info = open("%sspecies_info.txt" % midas_directory)
    species_info.readline() #header
    
    genome_id_to_return=''
    for line in species_info:
        items = line.split("\t")
        species_id = items[0].strip()
        rep_genome = items[1].strip()
        if desired_species_name == species_id:
            genome_id_to_return=rep_genome
            return genome_id_to_return


##########################################################
#
# Get a list of all genome_ids corresponding to reference genomes for a given species in the midas_db. 
#
#########################################################
def get_ref_genome_ids(desired_species_name):
    
    genome_ids=[]
    genome_info = open("%sgenome_info.txt" % midas_directory)
    genome_info.readline() #header
    for line in genome_info:
        items = line.split("\t")
        genome_id = items[0].strip()
        species_id=items[5].strip() 
        if desired_species_name == species_id:
            genome_ids.append(genome_id)
    return genome_ids


##########################################################
#
# Return a dictionary of all genome IDs and thier species id
#
#########################################################
def genome_ids_dictionary():
    
    genome_ids={}
    genome_info = open("%sgenome_info.txt" % midas_directory)
    genome_info.readline() #header
    for line in genome_info:
        items = line.split("\t")
        genome_id = items[0].strip()
        species_id=items[5].strip() 
        rep_genome=items[2]
        genome_ids[genome_id]=[species_id, rep_genome]
    return genome_ids



##########################################################
#
# parse the intermediate files of MIDAS to obtain CNV counts for 99% gene centroids
#
#########################################################
def parse_99_percent_genes(desired_species_name,samples, allowed_genes=[]):

    # dictionary to store all the data (key=sampleID, value={gene, numreads})
    data={}

    for sample in samples:
        if sample[-1]=='c':
            pollard_dir='/pollard/home/ngarud/BenNanditaProject/MIDAS_intermediate_files_hmp/MIDAS_1.2.2_samples_combined_output/'+sample+'/genes/output'
        else:
            pollard_dir='/pollard/home/ngarud/BenNanditaProject/MIDAS_intermediate_files_hmp/MIDAS_1.2.2_output/'+sample+'/genes/output'
        data[sample]={}

        file = gzip.open("%s/%s.genes.gz" % (pollard_dir, desired_species_name), 'r')
        file.readline() #header
        for line in file:
            items = line.split()
            gene = items[0]
            count_reads= items[1]            
            if gene in allowed_genes:
                data[sample][gene]=int(count_reads)
           
#    # create a numpy array with the data
    data_numpy_array_dict={}
    data_numpy_array=numpy.asarray([])
    for gene in allowed_genes:
        
        data_numpy_array_dict[gene]=[]
        #for sample in data.keys(): old version: may not preserve order of samples
        for sample in samples:
            if gene in data[sample]:
                data_numpy_array_dict[gene].append(data[sample][gene])
            else:
                data_numpy_array_dict[gene].append(0)
        data_numpy_array_dict[gene]=numpy.asarray(data_numpy_array_dict[gene])

    #consolidate into a single ref genome esp. if multiple genes are in a single ref. 
    ref_genome_dict={}
    for gene in allowed_genes:
        #modify the gene name to just have the refgenome
        ref_genome='.'.join(gene.split('.')[0:2])
        if ref_genome not in ref_genome_dict:
            ref_genome_dict[ref_genome]=data_numpy_array_dict[gene]
        else:
            ref_genome_dict[ref_genome]+=data_numpy_array_dict[gene]

    return ref_genome_dict
                
                

##########################################################
#
# parse the pre-merged species file to get a list of species at at least 3x coverage
#
#########################################################

def parse_intermediate_species_file(sample_id, inFN):
    inFile=open(inFN,'r')

    species_list=[]
    inFile.readline() #header
    for line in inFile:
        items=line.strip().split('\t')
        species_id=items[0]
        coverage=float(items[2])
        if coverage >=3.0: 
            species_list.append(species_id)

    return set(species_list)


#########################################
# 
# Read in centroids fasta sequences
#
#########################################
def load_centroid_fasta(species_name):
    centroid_fastas={} # key=gene, value=sequence
    centroid_file = gzip.open("%span_genomes/%s/centroids.ffn.gz" % (midas_directory, species_name), 'r')
    for line in centroid_file:
        line=line.strip()
        if line[0]=='>':
            gene_name=line[1:len(line)]
            centroid_fastas[gene_name]=''
        else:
            centroid_fastas[gene_name]+= line

    return centroid_fastas
    


#######################    

if __name__=='__main__':

    pass
    
