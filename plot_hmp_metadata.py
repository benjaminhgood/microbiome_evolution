import matplotlib  
matplotlib.use('Agg') 
import pylab
import numpy
import parse_midas_data
import stats_utils
import sys

combination_type="sample"

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

# Calculate num timepoints per sample
num_timepoints_per_subject = []
for subject in subject_sample_map.keys():
    num_timepoints_per_subject.append( len(subject_sample_map[subject].keys()) )

num_timepoints_per_subject.sort()
num_timepoints_per_subject = numpy.array(num_timepoints_per_subject)

num_timepoints, num_subjects = stats_utils.calculate_unnormalized_survival_from_vector(num_timepoints_per_subject)

pylab.figure(1,figsize=(5,3))
pylab.step(num_timepoints+0.25, num_subjects,where='pre')
pylab.semilogy([0],[1])
pylab.xlim([0.5,9])
pylab.ylim([0.3,300])
pylab.xlabel('Num timepoints, $T$')
pylab.ylabel('Num subjects with $\geq T$')
print len(num_timepoints_per_subject), max(num_timepoints_per_subject)
pylab.savefig('%s/num_timepoints_per_subject.pdf' % parse_midas_data.analysis_directory,bbox_inches='tight')


# Load marker gene coverages
species_coverage_matrix, sample_list, species_list = parse_midas_data.parse_global_marker_gene_coverages(combination_type=combination_type)
total_coverage_vector = species_coverage_matrix.sum(axis=0)
# Calculate frequencies
species_frequencies = species_coverage_matrix / (total_coverage_vector+(total_coverage_vector<0.5))[None,:]

pylab.figure(2,figsize=(5,3))
pylab.xlabel('Coverage per sample, $D$')
pylab.ylabel('Num samples >= $D$')
pylab.loglog([1],[1])
pylab.xlim([5,1000])
pylab.ylim([10,500])

pylab.figure(3,figsize=(5,3))
pylab.xlabel('Relative abundance, $f$')
pylab.ylabel('Num samples >= $f$')
pylab.loglog([1],[1])
pylab.xlim([1e-05,1])
pylab.ylim([10,500])

pylab.figure(4,figsize=(24,3))
pylab.xlabel('Sample')
pylab.ylabel('Relative abundance, $f$')
pylab.semilogy([-1],[1])
pylab.xlim([-0.5,len(species_list)-0.5])
pylab.ylim([1e-05, 1])
pylab.gca().set_xticklabels([])

allowed_species = set(["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318", "Bacteroides_ovatus_58035", "Prevotella_copri_61740"])


for species_idx in xrange(0,len(species_list)):
    species_name = species_list[species_idx]
    sample_coverages = species_coverage_matrix[species_idx,:]
    sample_frequencies = species_frequencies[species_idx,:]
    
    num_above_threshold = (sample_coverages>=10).sum()
    
    if num_above_threshold < 10:
        continue 
    
    print species_name, numpy.median(sample_coverages[sample_coverages>=10]), num_above_threshold
    
    
    
    #if species_name == "Bacteroides_uniformis_57318":
    if species_name == "Butyrivibrio_crossotus_61674":
    #if species_name.startswith('Prevotella'):
    #if species_name.startswith('Ruminococcus_torques'):
        print species_name, sample_coverages[sample_coverages>5], num_above_threshold
        linewidth=3
    else:
        linewidth=1

    coverages, coverage_survival_function = stats_utils.calculate_unnormalized_survival_from_vector(sample_coverages)
    frequencies, frequency_survival_function = stats_utils.calculate_unnormalized_survival_from_vector(sample_frequencies)
    
    frequencies[0] = 1e-06
    
    pylab.figure(2)
    line, = pylab.step(coverages+0.25, coverage_survival_function, where='pre',alpha=0.5,linewidth=linewidth, zorder=-species_idx)
    color = pylab.getp(line,'color')
    pylab.figure(3)
    pylab.step(frequencies, frequency_survival_function, where='pre',alpha=0.5,linewidth=linewidth,color=color, zorder=-species_idx)
    
    pylab.figure(4)
    if species_name in allowed_species:
        pylab.plot(numpy.arange(0,len(sample_frequencies)), sample_frequencies, '.-', alpha=0.5, linewidth=linewidth, color=color, zorder=-species_idx)
    
pylab.figure(2)
pylab.savefig('%s/species_marker_coverage_distribution.pdf' % parse_midas_data.analysis_directory,bbox_inches='tight')
pylab.figure(3)
pylab.savefig('%s/species_abundance_distribution.pdf' % parse_midas_data.analysis_directory,bbox_inches='tight')
pylab.figure(4)
pylab.savefig('%s/species_abundance_vector.pdf' % parse_midas_data.analysis_directory,bbox_inches='tight')
