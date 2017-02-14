import matplotlib  
matplotlib.use('Agg') 
import os
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
from numpy.random import choice
species_name=sys.argv[1]


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
sys.stderr.write("Done!\n")
    
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')
    
low_diversity_samples = (numpy.diag(avg_pi_matrix_syn)<1e-03)
    
unique_samples = parse_midas_data.calculate_unique_samples(subject_sample_map, samples)
    
desired_samples = unique_samples*low_diversity_samples

# initialize distance bins for LD computations
distance_bins = numpy.logspace(0,4,20) # bins start from 1 to 10^4 and there are 20 evenly spaced bins log(1)=0, log(10^4)-4
distance_bin_locations = numpy.array(distance_bins[:-1],copy=True) # shifted one to avoid edge effects for plotting.
distance_bins[0] = 0.5 # made smallest bin 0.5 to avoid edge effects
distance_bins[-1] = 1e09 # made largest bin very large to catch anything >10^4. 
    
low_freq=0.3
high_freq=0.5

sys.stderr.write("Calculating intra-gene LD...\n")

###############
# calculate LD
##############
binned_rsquareds_dict={}
control_rsquareds_dict={}
binned_rsquared_denominators_dict={}

for variant_type in ['1D','4D']:
    
    binned_rsquared_numerators = numpy.zeros_like(distance_bin_locations)
    binned_rsquared_denominators = numpy.zeros_like(distance_bin_locations)
    
    total_control_rsquared_numerators = 0
    total_control_rsquared_denominators = 0
    
    for gene_name in allele_counts_map.keys():
        locations = numpy.array([location for chromosome, location in allele_counts_map[gene_name][variant_type]['locations']])*1.0
        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']
        
        if len(allele_counts)==0:
            # no diversity to look at!
            continue
        
        # pick a random gene somewhere else as a control
        control_gene_name = gene_name
        control_allele_counts = []
        while gene_name==control_gene_name or len(control_allele_counts)==0:
            control_gene_name = choice(allele_counts_map.keys())
            control_allele_counts = allele_counts_map[control_gene_name][variant_type]['alleles']
        
        
        allele_counts = allele_counts[:,desired_samples,:]
        control_allele_counts = control_allele_counts[:,desired_samples,:]
        #compute the distances between all pairs of sites 
        # None in the two index positions results in a transpose of the vector relative to each other
        # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
        distances = numpy.fabs(locations[:,None]-locations[None,:])
  
        
        rsquared_numerators, rsquared_denominators = diversity_utils.calculate_rsquared_condition_freq(allele_counts, allele_counts, low_freq, high_freq)
        control_rsquared_numerators, control_rsquared_denominators = diversity_utils.calculate_rsquared_condition_freq(allele_counts, control_allele_counts, low_freq, high_freq)
        
        # get the indices of the upper diagonal of the distance matrix
        # numpy triu_indices returns upper diagnonal including diagonal
        # the 1 inside the function excludes diagonal. Diagnonal has distance of zero.
        desired_idxs = numpy.triu_indices(distances.shape[0],1)
        
        # fetch the distances and rsquared vals corresponding to the upper diagonal. 
        distances = distances[desired_idxs]
        rsquared_numerators = rsquared_numerators[desired_idxs]
        rsquared_denominators = rsquared_denominators[desired_idxs]
        
        # fetch entries where denominator != 0 (remember, denominator=pa*(1-pa)*pb*(1-pb). If zero, then at least one site is invariant)
        distances = distances[rsquared_denominators>0]
        rsquared_numerators = rsquared_numerators[rsquared_denominators>0] 
        rsquared_denominators = rsquared_denominators[rsquared_denominators>0]
        
        if len(distances) == 0:
            continue
            
        # numpy.digitize: For each distance value, return the bin index it belongs to in distances_bins. 
        bin_idxs = numpy.digitize(distances,bins=distance_bins)-1
            
        for i in xrange(0,len(bin_idxs)):
        
            binned_rsquared_numerators[bin_idxs[i]] += rsquared_numerators[i] 
            
            binned_rsquared_denominators[bin_idxs[i]] += rsquared_denominators[i]
        
        control_rsquared_numerators = control_rsquared_numerators[control_rsquared_denominators>0]
        control_rsquared_denominators = control_rsquared_denominators[control_rsquared_denominators>0]
        
        total_control_rsquared_numerators += (control_rsquared_numerators).sum()
        total_control_rsquared_denominators += (control_rsquared_denominators).sum()
            
    binned_rsquareds_dict[variant_type] = binned_rsquared_numerators/(binned_rsquared_denominators+(binned_rsquared_denominators==0))
    
    control_rsquareds_dict[variant_type] = total_control_rsquared_numerators/(total_control_rsquared_denominators+(total_control_rsquared_denominators==0))
    
    binned_rsquared_denominators_dict[variant_type]=binned_rsquared_denominators
    

# write to an intermediate file so that I can plot all species' LD decays on one plot (separate script: plot_gene_ld_multispecies.py)
#numpy.savez(os.path.expanduser('~/tmp_intermediate_files/LD_%s.npz' % species_name),binned_rsquareds=binned_rsquareds, binned_rsquared_denominators=binned_rsquared_denominators, distance_bin_locations=distance_bin_locations, control_rsquareds=control_rsquareds)

        
pylab.figure()
pylab.xlabel('Distance between SNPs')
pylab.ylabel("Ohta and Kimura's $\\sigma^2_d$")
pylab.suptitle("LD among common SNPs (MAF between 0.3 and 0.5)")
pylab.xlim([1,1e04])
pylab.ylim([1e-02,1])
variant_type='1D'
nonsyn=pylab.loglog(distance_bin_locations[binned_rsquared_denominators_dict[variant_type]>0], binned_rsquareds_dict[variant_type][binned_rsquared_denominators_dict[variant_type]>0],'k.-', label='nonsyn')
pylab.loglog(distance_bin_locations, numpy.ones_like(distance_bin_locations)*control_rsquareds_dict[variant_type],'k:')

variant_type='4D'
syn=pylab.loglog(distance_bin_locations[binned_rsquared_denominators_dict[variant_type]>0], binned_rsquareds_dict[variant_type][binned_rsquared_denominators_dict[variant_type]>0],'k.-', color='r', label='syn')
pylab.loglog(distance_bin_locations, numpy.ones_like(distance_bin_locations)*control_rsquareds_dict[variant_type],'k:', color='r')

handles=[nonsyn,syn]
labels=['nonsyn','syn']
pylab.legend(handles,labels,'lower left',prop={'size':6})  

pylab.savefig('%s/%s_intragene_ld_syn_nonsyn_common.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
pylab.savefig('%s/%s_intragene_ld_syn_nonsyn_common.png' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight', dpi=300)
    
#pylab.show()
    
