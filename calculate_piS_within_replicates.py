import parse_midas_data
import pylab
import sys
import numpy
from calculate_pi_matrix import calculate_self_pis
import os

species_names = ["Bacteroides_uniformis_57318"]
default_directory_prefix =  os.path.expanduser("~/Documents/files_too_big_for_Dropbox/midas_output_120816/")

for species_name in species_names:

    print species_name

    sys.stderr.write("Loading %s...\n" % species_name)
    
    samples, allele_counts_syn, locations_syn, genes_syn, passed_sites_syn, allele_counts_non, locations_non, genes_non, passed_sites_non = parse_midas_data.parse_snps(species_name, site_depth_threshold=15, directory_prefix=default_directory_prefix)
    
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Calculating pis...\n")
    piS = calculate_self_pis(allele_counts_syn)
    piS /= (passed_sites_syn+(passed_sites_syn==0))
    sys.stderr.write("Done!\n")
    
    for sample,pi in zip(samples,piS):
        print sample,pi