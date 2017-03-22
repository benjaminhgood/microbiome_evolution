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

########################################################################################
#
# Standard header to read in argument information
#
########################################################################################
if len(sys.argv)>1:
    if len(sys.argv)>2:
        debug=True # debug does nothing in this script
        sample_name = sys.argv[2]
    else:
        debug=False
        sample_name = sys.argv[1]
else:
    sys.stderr.write("Usage: python plot_within_sfs_for_sample.py [debug] sample_name")
########################################################################################


species_names = parse_midas_data.parse_good_species_list()

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")    

pylab.figure()

# Set up figure
fig = plt.figure(figsize=(10, 4))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[6,3], wspace=0.1)

###################
#
# SFS panel
#
###################

sfs_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(sfs_axis)

sfs_axis.set_xlabel('Allele frequency of global minor variant')
sfs_axis.set_ylabel('Fraction of sites')
sfs_axis.set_xlim([-0.03,1.03])


legend_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  
    
min_coverage = 10
allowed_variant_types = set(['4D'])

for species_name in species_names:

    # Load genomic coverage distributions for this species
    sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
    median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])    
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

    desired_sample = parse_midas_data.sample_name_lookup(sample_name, samples)

    if desired_sample=="":
        sys.stderr.write("%s not present!\n" % (species_name))
        continue
    
    # Make sure there's enough coverage to plot an SFS
    median_coverage = sample_coverage_map[desired_sample]
    if median_coverage < min_coverage:
        sys.stderr.write("%s at too low coverage in %s (D=%g)!\n" % (species_name, desired_sample, median_coverage))
        continue


    # Load SNP information for species_name
    sys.stderr.write("Loading %s (D=%g)...\n" % (species_name, median_coverage))
    samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug)
    sys.stderr.write("Done!\n")
    
    # Make sure sample is actually in the SNP matrix!
    if desired_sample not in samples:
        continue
            
    desired_idx = list(samples).index(desired_sample)
    
    # Now calculate all allele frequencies for this sample
    sample_freqs = []
    passed_sites = 0
    for variant_type in allowed_variant_types:
        for gene_name in allele_counts_map.keys():
    
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

            if len(allele_counts)==0:
                continue
            
            depths = allele_counts[:,desired_idx,:].sum(axis=1)
            freqs = allele_counts[:,desired_idx,0]*1.0/(depths+(depths==0))
            # Make them minor allele freqs
            #freqs = numpy.fmin(freqs,1-freqs)
            sample_freqs.extend( freqs[freqs>0] )
            passed_sites += passed_sites_map[gene_name][variant_type]['sites'][desired_idx,desired_idx]
        
    bins = numpy.arange(0,long(median_coverage)+1)*1.0/long(median_coverage)
    bins[0] = -0.01
    bins[-1] = 1.01
    #bins = bins[bins<=0.5]
    bin_locations = bins[1:]-(bins[2]-bins[1])/2
      
    counts,dummy = numpy.histogram(sample_freqs,bins)
    counts[0] += (passed_sites-len(sample_freqs))
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/counts.sum()
    
    xs = bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1
        
    line, = sfs_axis.semilogy(xs[sfs>0], sfs[sfs>0],'.-',alpha=0.5,markersize=2)
    colorVal = pylab.getp(line,'color')
    legend_axis.plot([-2,-1],[-2,-1],'.-',color=colorVal,alpha=0.5, label=('%s (D=%g)' % (species_name,median_coverage)))
    
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   

fig.savefig('%s/%s_sample_sfs.pdf' % (parse_midas_data.analysis_directory,sample_name),bbox_inches='tight')
fig.savefig('%s/%s_sample_sfs.png' % (parse_midas_data.analysis_directory,sample_name),bbox_inches='tight')
     
