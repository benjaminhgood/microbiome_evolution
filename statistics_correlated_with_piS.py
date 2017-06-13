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
mpl.rcParams['font.size'] = 13
mpl.rcParams['lines.linewidth'] = 4.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'large'




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

sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! %d core genes\n" % len(core_genes))

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
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_variant_types=set(['4D']), allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

######################
# compute median cov #
######################

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

##########################################################
# load SNP info
##########################################################

# note that this loads info for all samples. Later the desired samples are selected out. 
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug)
sys.stderr.write("Done!\n")


###################################################################
# compute 1D SFSs 
###################################################################                    

sys.stderr.write("Calculate within person SFS...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs(allele_counts_map, passed_sites_map, variant_type='4D',fold=True)

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
        
sfss=numpy.asarray(sfss)

sys.stderr.write("Done!\n")


###################################################################
# create subplots to plot 1D SFSs
###################################################################                    

#indexes of sorted pis array
idxs = numpy.argsort(pis)

#pylab.figure(figsize=(4, len(pis)*3))   
pylab.figure(figsize=(100,100))   
plot_no=1
#for j in range(0, len(pis)):     
for j in range(0,len(pis)):
    idx=idxs[j]
    # plot the 1D SFS
    pylab.subplot(18, 17, plot_no)
    #pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.title(pis[idx])
    pylab.xlim([0,0.5])
    pylab.ylim([3e-06,3e-02])
    if sfss[j].sum()!=0: 
        normalized_sfs = sfss[idx]
        pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=1)
    pylab.axhline(y=0.001, color='r', linestyle='-')
    plot_no+=1

pylab.savefig('%s/%s_1D_sfs_polarized.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')



###############################################################
# compute the fraction of SFS bins with >10^-3 SNPs
###############################################################

num_bins_surpassing_threshold=(sfss>0.001).sum(axis=1) 

pylab.figure(figsize=(8,8))
mpl.rcParams['font.size'] = 10
pylab.xlabel('pi_S')
pylab.ylabel('Number of SFS bins surpassing threshold')
pylab.title('Descriptive statistics of SFS vs piS')

pylab.plot(pis, num_bins_surpassing_threshold, "b.")
pylab.axvline(x=0.001, color='r', linestyle='-')
#pylab.plot([0,500],[0,500] , 'k-')
pylab.savefig('%s/descriptive_statistics_sfs_vs_piS.png' % (parse_midas_data.analysis_directory),bbox_inches='tight') 


##################################################################
# Compute Tajima's D
################################################################## 
denom=[1]
for i in range(2,100):
    denom.append(denom[i-2]+ 1/float(i))

TajD=[]
for j in xrange(0,len(samples)): 
    TajD.append((pis[j]- len(sample_freqs[j])/passed_sites[j]/denom[98]))

# plot
pylab.figure(figsize=(8,8))
mpl.rcParams['font.size'] = 10
pylab.xlabel('pi_S')
pylab.ylabel("Tajima's D")
pylab.title("Tajima's D vs piS")
pylab.plot(pis, TajD, "b.")
pylab.axvline(x=0.001, color='r', linestyle='-')
pylab.savefig('%s/TajD_vs_piS.png' % (parse_midas_data.analysis_directory),bbox_inches='tight') 
    
