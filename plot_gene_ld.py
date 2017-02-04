import matplotlib  
matplotlib.use('Agg') 
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
    
distance_bins = numpy.logspace(0,4,20)
distance_bin_locations = numpy.array(distance_bins[:-1],copy=True)
distance_bins[0] = 0.5
distance_bins[-1] = 1e09
    
binned_rsquared_numerators = numpy.zeros_like(distance_bin_locations)
binned_rsquared_denominators = numpy.zeros_like(distance_bin_locations)
    
total_control_rsquared_numerators = 0
total_control_rsquared_denominators = 0
    
sys.stderr.write("Calculating intra-gene LD...\n")
# calculate SFS
for gene_name in allele_counts_map.keys():
        
    locations = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['4D']['locations']])*1.0
    allele_counts = allele_counts_map[gene_name]['4D']['alleles']
        
    if len(allele_counts)==0:
        # no diversity to look at!
        continue
        
    # pick a random gene somewhere else as a control
    control_gene_name = gene_name
    control_allele_counts = []
    while gene_name==control_gene_name or len(control_allele_counts)==0:
        control_gene_name = choice(allele_counts_map.keys())
        control_allele_counts = allele_counts_map[control_gene_name]['4D']['alleles']
        
        
    allele_counts = allele_counts[:,desired_samples,:]
    control_allele_counts = control_allele_counts[:,desired_samples,:]
        
    distances = numpy.fabs(locations[:,None]-locations[None,:])
    
    rsquared_numerators, rsquared_denominators = diversity_utils.calculate_rsquared(allele_counts, allele_counts)
    control_rsquared_numerators, control_rsquared_denominators = diversity_utils.calculate_rsquared(allele_counts, control_allele_counts)
        
    desired_idxs = numpy.triu_indices(distances.shape[0],1)
        
    distances = distances[desired_idxs]
    rsquared_numerators = rsquared_numerators[desired_idxs]
    rsquared_denominators = rsquared_denominators[desired_idxs]
        
    distances = distances[rsquared_denominators>0]
    rsquared_numerators = rsquared_numerators[rsquared_denominators>0]
    rsquared_denominators = rsquared_denominators[rsquared_denominators>0]
        
    if len(distances) == 0:
        continue
            
    bin_idxs = numpy.digitize(distances,bins=distance_bins)-1
            
    for i in xrange(0,len(bin_idxs)):
        
        binned_rsquared_numerators[bin_idxs[i]] += rsquared_numerators[i]
            
        binned_rsquared_denominators[bin_idxs[i]] += rsquared_denominators[i]
        
    control_rsquared_numerators = control_rsquared_numerators[control_rsquared_denominators>0]
    control_rsquared_denominators = control_rsquared_denominators[control_rsquared_denominators>0]
        
    total_control_rsquared_numerators += (control_rsquared_numerators).sum()
    total_control_rsquared_denominators += (control_rsquared_denominators).sum()
            
binned_rsquareds = binned_rsquared_numerators/(binned_rsquared_denominators+(binned_rsquared_denominators==0))
    
control_rsquareds = total_control_rsquared_numerators/(total_control_rsquared_denominators+(total_control_rsquared_denominators==0))
        
pylab.figure()
pylab.xlabel('Distance between SNPs')
pylab.ylabel("Ohta and Kimura's $\\sigma^2_d$")
pylab.xlim([1,1e04])
pylab.ylim([1e-02,1])
pylab.loglog(distance_bin_locations[binned_rsquared_denominators>0], binned_rsquareds[binned_rsquared_denominators>0],'k.-')
pylab.loglog(distance_bin_locations, numpy.ones_like(distance_bin_locations)*control_rsquareds,'k:')
pylab.savefig('%s/%s_intragene_ld.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
pylab.savefig('%s/%s_intragene_ld.png' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight', dpi=300)
    
#pylab.show()
    
