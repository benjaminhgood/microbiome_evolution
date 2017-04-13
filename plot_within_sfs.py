import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
from numpy.random import normal
#from calculate_pi_matrix import calculate_self_pis
import diversity_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint


mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

vmin = log10(1e-04)
vmax = log10(1e-02)
cmap='jet'
jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
pi_scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

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
################################################################################


min_coverage = 20

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
       
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
samples = numpy.array(samples)

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold
snp_samples = samples[(median_coverages>=min_coverage)]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)

binss=None
bin_locationss=None
countss=None
pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])
pi_opportunities_syn = numpy.array([])
snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])
desired_median_coverages = numpy.array([])

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    desired_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    if binss==None:
        desired_median_coverages = numpy.array([sample_coverage_map[desired_samples[i]] for i in xrange(0,len(desired_samples))])

        binss = [numpy.arange(0,long(desired_median_coverages[i])+1)*1.0/long(desired_median_coverages[i]) for i in xrange(0,len(desired_samples))]
        for i in xrange(0,len(binss)):
            binss[i][0] = -0.01
            binss[i][-1] = 1.01
        
        bin_locationss = [binss[i][1:]-(binss[i][2]-binss[i][1])/2 for i in xrange(0,len(desired_samples))]
        
        countss = [numpy.zeros_like(bin_locationss[i]) for i in xrange(0,len(desired_samples))]

        pi_matrix_syn = numpy.zeros((len(desired_samples), len(desired_samples)))*1.0
        avg_pi_matrix_syn = numpy.zeros_like(pi_matrix_syn)
        pi_opportunities_syn = numpy.zeros_like(pi_matrix_syn)
        
    # Calculate full matrix of synonymous pairwise differences
    sys.stderr.write("Calculate synonymous pi matrix...\n")
    chunk_pi_matrix_syn, chunk_avg_pi_matrix_syn, chunk_passed_sites = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')

    pi_matrix_syn += chunk_pi_matrix_syn
    avg_pi_matrix_syn += chunk_avg_pi_matrix_syn
    pi_opportunities_syn += chunk_passed_sites
    

    sys.stderr.write("Calculate within person SFS...\n")
    chunk_sample_freqs, chunk_passed_sites = diversity_utils.calculate_sample_freqs( allele_counts_map, passed_sites_map, variant_type='4D')

    for i in xrange(0,len(desired_samples)):
        
        chunk_counts,dummy = numpy.histogram(chunk_sample_freqs[i],binss[i])
        chunk_counts[0] += (chunk_passed_sites[i]-len(chunk_sample_freqs[i]))
        
        countss[i] += chunk_counts
        
sys.stderr.write("Done!\n")   

avg_pi_matrix_syn = avg_pi_matrix_syn/(pi_opportunities_syn+(pi_opportunities_syn==0))

avg_pis = numpy.diag(avg_pi_matrix_syn)

# Now sort based on avg pi
sorted_idxs = numpy.arange(0,len(avg_pis))

avg_pis, sorted_idxs = (numpy.array(x) for x in zip(*sorted(zip(avg_pis,sorted_idxs),reverse=True)))

# This figure has everything on the same one
pylab.figure(1)
pylab.xlabel('Minor allele frequency')
pylab.ylabel('SFS')
pylab.xlim([0,0.5])
pylab.ylim([3e-06,3e-02])
pylab.title(species_name)
   
# This figure spreads them all out
figure_width=3.42
figure_height=2*len(avg_pis)
pylab.figure(2,figsize=(figure_width,figure_height))
fig = pylab.gcf()

# make three panels panels
outer_grid  = gridspec.GridSpec(len(avg_pis), 1, height_ratios=([1]*len(avg_pis)), hspace=0.1)

for i in xrange(0,len(avg_pis)):

    pi = avg_pis[i]
    D = desired_median_coverages[sorted_idxs[i]]
    counts = numpy.array(countss[sorted_idxs[i]])
    bin_locations = bin_locationss[sorted_idxs[i]]
    
    sfs = counts/counts.sum()
    
    # First do figure 1
    colorVal = pi_scalarMap.to_rgba(numpy.log10(pi))
    
    pylab.figure(1)
    pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, sfs,'.-',alpha=0.5, color=colorVal)
    
    
    # Then do figure 2
    pylab.figure(2)
    
    axis = plt.Subplot(fig, outer_grid[i])
    fig.add_subplot(axis)
    axis.set_xlim([0,0.5])
    axis.set_xticks(numpy.arange(0,11)*0.05)
    axis.set_xticklabels([])
    axis.set_ylim([0,counts[1:].max()+1])
    
    axis.plot(bin_locations[1:], counts[1:], label=('%d: pi=%g, D=%g' % (i, pi, D)))
    axis.legend(frameon=False,fontsize=7)
     
    
pylab.figure(1)   
m = pylab.scatter([-1],[1],c=[-1], vmin=vmin, vmax=vmax, cmap=cmap,    marker='^')
    
fig = pylab.gcf()
cax = fig.add_axes([0.95, 0.1, 0.02, 0.80])
cbar = fig.colorbar(m,cax=cax,orientation='vertical',ticks=[-4,-3,-2])
    
cbar.set_ticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$'])
#cl = pylab.getp(cbar.ax, 'ymajorticklabels')
#pylab.setp(cl, fontsize=9) 
fig.text(0.947,0.035,'$\\pi_s$',fontsize=12)
    
pylab.savefig('%s/%s_within_person_sfs.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_within_person_sfs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

pylab.figure(2)
pylab.savefig('%s/%s_within_person_sfs_separate.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_within_person_sfs_separate.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
  
