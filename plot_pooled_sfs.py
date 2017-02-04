import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
species=sys.argv[1]



# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species, combination_type="sample", debug=True)
sys.stderr.write("Done!\n")
    
sys.stderr.write("Calculating synonymous SFS...\n")
# calculate SFS
pooled_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, variant_type='4D')
pooled_freqs = numpy.fmin(pooled_freqs,1-pooled_freqs)
    
bins = numpy.linspace(0,0.5,51)
bins -= (bins[1]-bins[0])/2
xs = bins[1:]-(bins[1]-bins[0])/2
    
sfs_syn, dummy = numpy.histogram(pooled_freqs, bins=bins) 
    
sys.stderr.write("Calculating nonsynonymous SFS...\n")
pooled_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, variant_type='1D')
pooled_freqs = numpy.fmin(pooled_freqs,1-pooled_freqs)
    
sfs_non, dummy = numpy.histogram(pooled_freqs, bins=bins) 
    
pylab.figure()
pylab.xlabel('Minor allele frequency $f$')
pylab.ylabel('$f(1-f) p(f)$')
pylab.title(species)

pylab.plot(xs, sfs_syn*xs*(1-xs), 'b.-',label='synonymous (4D)')
pylab.plot(xs, sfs_non*xs*(1-xs), 'r.-',label='nonsynonymous (1D)')
pylab.legend(loc='upper right',frameon=False)
pylab.savefig('%s/%s_pooled_sfs.pdf' % (parse_midas_data.analysis_directory, species), bbox_inches='tight')
pylab.savefig('%s/%s_pooled_sfs.png' % (parse_midas_data.analysis_directory, species), bbox_inches='tight', dpi=300)
    
#pylab.show()
    
