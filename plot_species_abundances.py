import numpy
import matplotlib  
matplotlib.use('Agg') 
import pylab
import parse_midas_data

# Load marker gene coverages
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
total_coverages = species_coverage_matrix.sum(axis=0)

for j in xrange(0,len(samples)):
    if samples[j].startswith('SRR'):
        print "Kuleshov coverages:"
        for i in xrange(0,len(species)):
            species_coverage_list = species_coverage_matrix[i,:]
            print species[i], species_coverage_list[j], numpy.median(species_coverage_list[species_coverage_list>=1])
        
# Only keep those samples with nonzero coverage!
nonzero_sample_idxs = (total_coverages>0.5)
species_coverage_matrix = species_coverage_matrix[:,nonzero_sample_idxs]
total_coverages = total_coverages[nonzero_sample_idxs]

print "Species to keep!"
good_species=0
for i in xrange(0,len(species)):
    coverages = species_coverage_matrix[i,:]
    if len(coverages[coverages>=5]) >= 10:
        good_species+=1
        print species[i], len(coverages[coverages>=1]), numpy.median(coverages[coverages>=1])
    
print good_species, "good species"

species_abundance_matrix = species_coverage_matrix / total_coverages[None,:]

pylab.figure(figsize=(14,2))
for i in xrange(0,10):
    pylab.semilogy(species_abundance_matrix[i,:],'-')
    
pylab.savefig('%sspecies_abundances.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')