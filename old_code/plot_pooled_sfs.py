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
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species, debug=False)
sys.stderr.write("Done!\n")
    
sys.stderr.write("Calculating synonymous SFS...\n")
# calculate SFS
pooled_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_variant_types=['4D'])
pooled_freqs = numpy.fmin(pooled_freqs,1-pooled_freqs)
    
bins = numpy.linspace(0,0.5,51)
bins -= (bins[1]-bins[0])/2
xs = bins[1:]-(bins[1]-bins[0])/2
    
sfs_syn, dummy = numpy.histogram(pooled_freqs, bins=bins) 
    
sys.stderr.write("Calculating nonsynonymous SFS...\n")
pooled_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_variant_types=['1D'])
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
    
# plot the observed vs expected (1/x) SFS

pylab.figure()
pylab.xlabel('Minor allele frequency $f$')
pylab.ylabel('fraction of snps')
pylab.title(species)

observed_syn=sfs_syn[1:]/float(sum(sfs_syn[1:]))
observed_non=sfs_non[1:]/float(sum(sfs_syn[1:]))
expected=1/(xs[1:]*(1-xs[1:]))/(float(sum(1/(xs[1:]*(1-xs[1:])))))

pylab.semilogy(xs[1:], observed_syn, label='observed synonymous (4D)', color='b')
pylab.semilogy(xs[1:], observed_non, label='observed non synonymous (1D)', color='r')
pylab.semilogy(xs[1:], expected, label='expected (1/(f(1-f)))', color='k')
pylab.legend(loc='upper right',frameon=False)
pylab.savefig('%s/%s_pooled_sfs_obs_exp.pdf' % (parse_midas_data.analysis_directory, species), bbox_inches='tight')
