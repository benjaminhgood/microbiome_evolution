#!/usr/bin/env python 
### This script runs the necessary post-processing of the MIDAS output 
### across all species in serial. 

import os
import sys
import config

if len(sys.argv) > 1:
    argument=sys.argv[1]
else:
    argument = 'all'

# First calculate core genes for each species
#os.system('python core_gene_utils.py')

# Call postprocess_midas_data.py for each species
os.system('python loop_over_species_wrapper.py %s python postprocess_midas_data.py' % argument)

# Calculate substitution rates for the most prevalent species
#os.system('python calculate_substitution_rates.py')

# Calculate temporal changes for the most prevalent species
#os.system('python calculate_private_snvs.py')

# Calculate temporal changes for the most prevalent species
#os.system('python calculate_temporal_changes.py')

# Calculate linkage disequilibria for the most prevalent species
#os.system('python calculate_singletons.py')

# This one has to go here because we have to estimate clades first...
# Calculate linkage disequilibria for the most prevalent species
#os.system('python calculate_linkage_disequilibria.py')
