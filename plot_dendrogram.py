import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
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
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates
allowed_variant_types = set(['1D','2D','3D','4D'])

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sample_country_map = parse_midas_data.parse_sample_country_map()
sys.stderr.write("Done!\n")

# Load core gene set
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
 
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

print "All samples:", len(samples)

# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
# Restrict to single timepoint single timepoints per person

print "Haploid samples:", len(snp_samples)

unique_subject_idxs = parse_midas_data.calculate_unique_samples(subject_sample_map, snp_samples)
snp_samples = snp_samples[unique_subject_idxs]

print "Uniqued num samples:", len(snp_samples)

united_states_idxs = parse_midas_data.calculate_country_samples(sample_country_map, snp_samples, allowed_countries=set(['United States']))

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %d samples for %s...\n" % (len(snp_samples), species_name))

snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_genes=core_genes, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=allowed_variant_types, allowed_genes=core_genes, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
     
substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

substitution_rate = numpy.clip(substitution_rate,1e-09,10)

# calculate compressed distance matrix suitable for agglomerative clustering
Y = []
for i in xrange(0,substitution_rate.shape[0]):
    for j in xrange(i+1,substitution_rate.shape[1]):
        Y.append(substitution_rate[i,j]) 
Y = numpy.array(Y) 
    
Z = linkage(Y, method='average')        
    
c, coph_dists = cophenet(Z, Y)
sys.stderr.write("cophenetic correlation: %g\n" % c)

# Make a dendrogram
pylab.figure(1, figsize=(15, 5))
pylab.title('UPGMA dendrogram for %s' % species_name)
pylab.xlabel('sample index')
pylab.ylabel('distance')
ddata = dendrogram(Z, no_plot=True)

# calculate second minimum y value
ys = []
xs = []
for i, d in zip(ddata['icoord'], ddata['dcoord']):
    ys.extend(d)
    xs.extend(i)
    
xs = list(set(xs))
xs.sort()
xs = numpy.array(xs)

dx = xs[-1]-xs[0]
xmin = xs[0]-dx*0.05
xmax = xs[-1]+dx*0.05

ys = list(set(ys))
ys.sort()
ys = numpy.array(ys)

if ys[0]<1e-09:
    y_penultimin = ys[1]/2
else:
    y_penultimin = ys[0]/2

y_penultimax = ys[-1]

ymin=2e-10
ymax=1e-01

yplotmin = 1e-06
yplotmax = 1e-01


#print ymin

leaf_xs = []

for icoord, dcoord in zip(ddata['icoord'], ddata['dcoord']):
    for idx in xrange(0,len(icoord)-1):
        x0 = icoord[idx]
        y0 = dcoord[idx]
        if y0<1e-10:
            y0 = ymin
        x1 = icoord[idx+1]
        y1 = dcoord[idx+1]
        if y1<1e-10:
            y1 = ymin
        
        if (y0==ymin):
            leaf_xs.append(x0)
        
        if (y1==ymin):
            leaf_xs.append(x1)
        
        #print x0, '->', x1, '\t',y0, '->', y1       
        pylab.semilogy([x0,x1],[y0,y1],'b-')
        
        if (y0==y_penultimax) and (y1==y_penultimax):
            # it's the cross bar that bridges the two most-diverged clades
            # so plot a root branch to the top of the plot
            xavg = (x0+x1)*0.5
            
            pylab.semilogy([xavg,xavg],[y_penultimax, ymax],'b-')

leaf_xs = list(sorted(set(leaf_xs)))

print len(leaf_xs)
print len(ddata['ivl'])

print ddata['leaves']


for i in xrange(0,len(ddata['ivl'])):
    
    idx = long(ddata['ivl'][i])
    x = leaf_xs[i]
    y = yplotmin
    sample = snp_samples[idx]
    
    
    if united_states_idxs[idx]:
        color = 'b'
    else:
        color = 'r'
        
    pylab.plot([x],[y],'o',color=color,markeredgewidth=0)
    
pylab.xticks([])
pylab.xlim([xmin,xmax])
pylab.ylim([yplotmin,yplotmax])
pylab.xlabel('Samples')
pylab.ylabel('SNP divergence')

pylab.savefig('%s/%s_full_dendrogram.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_full_dendrogram.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

