import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils

species_names = ["Bacteroides_uniformis_57318"]

for species_name in species_names:

    # Load subject and sample metadata
    sys.stderr.write("Loading HMP metadata...\n")
    subject_sample_map = parse_midas_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    # Load SNP information for species_name
    sys.stderr.write("Loading %s...\n" % species_name)
    samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
    sys.stderr.write("Done!\n")
    
    # Calculate full matrix of synonymous pairwise differences
    sys.stderr.write("Calculate synonymous pi matrix...\n")
    pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')
    
    within_pis = numpy.diag(pi_matrix_syn)
    idxs = numpy.arange(0,len(within_pis))
    
    within_pis, idxs = zip(*sorted(zip(within_pis, idxs)))
    
    for idx in idxs:
        print "\t".join(["%g" % D for D in pi_matrix_syn[idx,idxs]])
    