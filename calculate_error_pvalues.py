import parse_midas_data
import sys
import os

#default_directory_prefix =  os.path.expanduser("~/Documents/files_too_big_for_Dropbox/midas_output_121816/")
default_directory_prefix = parse_midas_data.default_directory_prefix
combination_type = "sample"
snp_prefix = parse_midas_data.get_snp_prefix_from_combination_type(combination_type)

#species_names = ["Bacteroides_uniformis_57318"]
#species_names = ["Roseburia_intestinalis_56239", "Eubacterium_eligens_61678"]
species_names = ["Prevotella_copri_61740"]


for species in species_names:

    sys.stderr.write("Processing %s...\n" % species)
 
    output_filename = "%ssnps/%s/%sannotated_snps.txt.bz2" % (default_directory_prefix, species, snp_prefix)
 
    os.system('python pipe_midas_data.py %s | ./annotate_pvalue | bzip2 -c > %s' % (species,output_filename) )  
 