import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
from numpy.random import choice
import os
import stats_utils
import os.path
import config
################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size

#############################################
# Minimum median coverage of sample to look at
min_coverage = config.min_median_coverage

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
        
# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
# Only consider one sample per person
desired_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=desired_samples)
sys.stderr.write("Done!\n")

######################
# Load coverage data #
######################

# Load genomic coverage distributions
#sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
#median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
#sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
#median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))]) # I don't follow yet why this shows up 2x differently. 


###############################################################
# Compute Pi within patients to figure out which are haploid  #
###############################################################

# identify which samples have low piS -- use these to construct haplotypes
#pi_matrix_syn, avg_pi_matrix_syn, passed_sites = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')
#low_diversity_samples = (numpy.diag(avg_pi_matrix_syn/(passed_sites+(passed_sites==0)))<1e-03)# the diagonal of the above has the within-patient pi value. 
#unique_samples = parse_midas_data.calculate_unique_samples(subject_sample_map, samples)
#desired_samples = unique_samples*low_diversity_samples

for gene_name in allele_counts_map.keys():
#for gene_name in [allele_counts_map.keys()[0]]:
    locations_4D = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['4D']['locations']])*1.0
    locations_1D = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['1D']['locations']])*1.0

    #check if there are at least 100 bps in this gene
    gene_length=len(locations_4D) + len(locations_1D)
    if gene_length >=100:

        # create a dictionary keeping track of whether a location is 1D or 4D (key=location, value=[index,1D/4D]).
        location_dictionary={}
        for loc in range(0, len(locations_4D)):
            location_dictionary[locations_4D[loc]]=[loc,'4D']
        for loc in range(0, len(locations_1D)):
            location_dictionary[locations_1D[loc]]=[loc,'1D']

        allele_counts_4D = allele_counts_map[gene_name]['4D']['alleles']
        allele_counts_1D = allele_counts_map[gene_name]['1D']['alleles']

        #allele_counts_4D = allele_counts_4D[:,desired_samples,:]
        #allele_counts_1D = allele_counts_1D[:,desired_samples,:]

        # Create two files for a gene haplotype:
        # (1) A file with the haplotypes themselves
        # (1) annotation file with numbers 0,1,2,3,4 to indicate:
    
        #0: p=0, site is unmutated
        #1: p=1, site is a fixed syn difference (4D)
        #2: p=1, site is a fixed nonsyn difference (1D)
        #3: p<1, site is a polymorphic syn within patient (4D)
        #4: p>1, site is a polymorphic nonsyn within patient (1D)

        # generate_haplotype produces the following files:
        # (1) tmp_consensus.txt
        # (2) tmp_anno.txt
        diversity_utils.generate_haplotype(allele_counts_4D, allele_counts_1D, location_dictionary, species_name)

        # create consensus allele file and an annotation file. 

        # Cluster the haplotypes by identity
        num_samples=len(desired_samples)
        os.system('python ~/ben_nandita_hmp_scripts/H12_H2H1_MIDAS.py ~/tmp_intermediate_files/tmp_consensus_'+ species_name +'.txt ' + str(num_samples) + ' -o ~/tmp_intermediate_files/tmp_cluster_'+ species_name +'.txt -g ' + gene_name)  

        # Plot the haplotypes with R
        os.system('mkdir -p ~/ben_nandita_hmp_analysis/hap_plots/' + species_name)
        os.system('Rscript ~/ben_nandita_hmp_scripts/visualizeHaplotypesMicrobiome3.R ~/tmp_intermediate_files/tmp_cluster_'+ species_name+'.txt '+ str(num_samples) + ' ' +  str(gene_length) + ' ~/tmp_intermediate_files/tmp_consensus_'+ species_name+'.txt ~/tmp_intermediate_files/tmp_anno_' + species_name +'.txt ~/ben_nandita_hmp_analysis/hap_plots/' + species_name +'/hap_plot_' + gene_name +'.pdf ' + gene_name)
