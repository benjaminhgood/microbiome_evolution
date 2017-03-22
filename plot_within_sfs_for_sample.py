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
pylab.xlabel('Minor allele frequency')
pylab.ylabel('SFS')
pylab.xlim([0,0.5])
#pylab.ylim([3e-06,3e-02])
    
min_coverage = 20
allowed_variant_types = set(['4D'])

for species_name in species_names:

    # Load genomic coverage distributions for this species
    sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
    median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])    
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

    desired_sample = parse_midas_data.sample_name_lookup(sample_name, samples)

    if desired_sample=="":
        sys.stderr.write("Sample %s not present for species %s!\n" % (sample_name, species_name))
        continue
    
    # Make sure there's enough coverage to plot an SFS
    median_coverage = sample_coverage_map[desired_sample]
    if median_coverage < min_coverage:
        continue

    # Load SNP information for species_name
    sys.stderr.write("Loading %s...\n" % species_name)
    samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug)
    sys.stderr.write("Done!\n")
    
    # Make sure sample is actually in the SNP matrix!
    if desired_sample not in samples:
        continue
            
    desired_idx = list(samples).index(desired_sample)
    
    # Now calculate all allele frequencies for this sample
    sample_freqs = []
    for variant_type in allowed_variant_types:
        for gene_name in allele_counts_map.keys():
    
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

            if len(allele_counts)==0:
                continue
            
            depths = allele_counts[:,desired_idx,:].sum(axis=1)
            freqs = allele_counts[:,desired_idx,0]*1.0/(depths+(depths==0))
            # Make them minor allele freqs
            freqs = numpy.fmin(freqs,1-freqs)
            sample_freqs.extend( freqs[freqs>0] )
        
    
    bins = numpy.arange(0,long(median_coverage)+1)*1.0/long(median_coverage)
    bins = bins[bins<=0.5]
    bin_locations = bins[1:]-(bins[1]-bins[0])/2
      
    counts,dummy = numpy.histogram(sample_freqs,bins)
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/counts.sum()
        
    pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, sfs,'.-',alpha=0.5,label=species_name)
    
pylab.legend(loc='center right',frameon=False,fontsize='8')
pylab.savefig('%s/%s_sample_sfs.pdf' % (parse_midas_data.analysis_directory,sample_name),bbox_inches='tight')
pylab.savefig('%s/%s_sample_sfs.png' % (parse_midas_data.analysis_directory,sample_name),bbox_inches='tight')
     
