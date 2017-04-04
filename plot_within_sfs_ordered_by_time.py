import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os


# plotting tools
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'




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

#############################################################################

# Minimum median coverage of sample to look at
min_coverage = 20

#################
# Load metadata #
#################

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load time metadata
subject_sample_time_map_all_samples = parse_midas_data.parse_subject_sample_time_map()

######################
# Load coverage data #
######################

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
# prune time meta data so that the highest coverage sample is retained for those subjects with >1 sample per time pt
subject_sample_time_map = parse_midas_data.prune_subject_sample_time_map(subject_sample_time_map_all_samples,sample_coverage_map)


###############################################################
# Compute Pi within patients to figure out which are haploid  #
###############################################################

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities
clipped_pis = (total_pis+1)/(total_pi_opportunities+1)

######################
# compute median cov #
######################

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

###############################################################
# Indexes for SNP samples that have high coverage #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
high_cov_samples = samples[(median_coverages>=min_coverage)]
#high_cov_pis     = pis[(median_coverages>=min_coverage)]

# get the time info for the snp_samples (this is only for visno 1 vs visno 2 or 3
#time_pair_idxs, visno, day = parse_midas_data.calculate_time_pairs(subject_sample_time_map, high_cov_samples)

# get the time info for the snp_samples -- this is for all visno combos
time_pair_idxs, visno1, visno2, day = parse_midas_data.calculate_all_time_pairs(subject_sample_time_map, high_cov_samples)


#### time pair idxs where patients can have exactly 1 time point (so that points plotted are iid)
#time_pair_idxs_unique, visno_snps_genes_unique, day_snps_genes_unique = parse_midas_data.calculate_unique_time_pairs(subject_sample_time_map, high_cov_samples)


### different patient idx: 
# to compare results to time_pair idxs, we want different patient pair idxs. This helps us to contextualize if we are seeing events within patients that resemble replacements or modifications. 

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
#snp_same_sample_idxs, snp_same_subject_idxs, snp_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, high_cov_samples)





##########################################################
# load SNP info
# compute folded SFS
##########################################################

sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug)
sys.stderr.write("Done!\n")

#folded SFS
sys.stderr.write("Calculate within person SFS...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs( allele_counts_map, passed_sites_map, variant_type='4D')

sfss = []
#pi_withins = []
bins = numpy.linspace(0.04,0.51,11)
bin_locations = bins[1:]-(bins[1]-bins[0])/2
    
for j in xrange(0,len(samples)):
    
    #pi_within = pis[j]
        
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
        
    sfss.append(sfs)
    #pi_withins.append(pi_within)
        

sys.stderr.write("Done!\n")

##########
# Plot:  #
##########

# Folded SFSs

# iterate through time pairs

color=['b','r']
for j in range(0, len(time_pair_idxs[0])):
    pylab.figure()
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,0.5])
    pylab.ylim([3e-06,3e-02])
    pylab.title(species_name)
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=samples[idx1]
    sample_name2=samples[idx2]
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    
    pylab.savefig('%s/%s_within_person_sfs_time_pair_folded_%s_%s.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')


#################
#polarized SFS  #
#################

sys.stderr.write("Calculate within person SFS...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs( allele_counts_map, passed_sites_map, variant_type='4D',fold=False)

sfss= []
bins = numpy.linspace(0.04,0.95,21)
bin_locations = bins[1:]-(bins[1]-bins[0])/2
    
for j in xrange(0,len(samples)):
    
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
        
    sfss.append(sfs)
        

sys.stderr.write("Done!\n")



# plot Polarized SFSs

# iterate through time pairs

color=['b','r']
for j in range(0, len(time_pair_idxs[0])):
    pylab.figure()
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,1])
    pylab.ylim([3e-06,3e-02])
    pylab.title(species_name)
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=samples[idx1]
    sample_name2=samples[idx2]
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    
    pylab.savefig('%s/%s_within_person_sfs_time_pair_polarized_%s_%s.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')






###################################################################
# Plot 2d SFS for time pairs (polarized)
###################################################################

# these sfs are polarized based on teh consensus allele. 
    
sample_freqs_2D, passed_sites_2D = diversity_utils.calculate_sample_freqs_2D( allele_counts_map, passed_sites_map, variant_type='4D', fold=False)


xbins = numpy.linspace(0,1,21) 
ybins = numpy.linspace(0,1,21)

for j in range(0, len(time_pair_idxs[0])):
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=samples[idx1]
    sample_name2=samples[idx2]
    counts, xbins, ybins = numpy.histogram2d(sample_freqs_2D[idx1], sample_freqs_2D[idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.figure()
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(species_name)
    im=pylab.imshow(sfs_2D,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1e-2), cmap='jet')
    pylab.colorbar(im)
    pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_polarized_%s_%s.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')




###################################################################
# Plot 2d SFS for time pairs (folded)
###################################################################

    
sample_freqs_2D, passed_sites_2D = diversity_utils.calculate_sample_freqs_2D( allele_counts_map, passed_sites_map, variant_type='4D')


xbins = numpy.linspace(0,1.1,21) 
ybins = numpy.linspace(0,1.1,21)

for j in range(0, len(time_pair_idxs[0])):
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=samples[idx1]
    sample_name2=samples[idx2]
    counts, xbins, ybins = numpy.histogram2d(sample_freqs_2D[idx1], sample_freqs_2D[idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.figure()
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(species_name)
    im=pylab.imshow(sfs_2D,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1e-2), cmap='jet')
    pylab.colorbar(im)
    pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_folded_%s_%s.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')

