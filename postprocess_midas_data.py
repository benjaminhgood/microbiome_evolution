### This script runs the necessary post-processing of the MIDAS output so that we can start analyzing
import os
import sys

#species_names=['Prevotella_copri_61740','Butyrivibrio_crossotus_61674','Bacteroides_uniformis_57318','Roseburia_intestinalis_56239','Eubacterium_eligens_61678']

species_names=['Butyrivibrio_crossotus_61674']

for species in species_names:
    sys.stderr.write('species: ' + species + '\n')

    sys.stderr.write('Calculating species-specific marker gene coverages...\n')
    os.system('python ~/projectBenNandita/calculate_marker_gene_coverage.py ' + species)   
    sys.stderr.write('Done!\n')

    sys.stderr.write('Merging technical replicates from same timepoint...\n')
    os.system('python ~/projectBenNandita/combine_replicates.py' + species)   
    sys.stderr.write('Done!\n')

    sys.stderr.write('Calculating coverage distributions...\n')
    os.system('python ~/projectBenNandita/calculate_coverage_distribution.py' + species)
    sys.stderr.write('Done!\n')

    # Calculate core genome set 
    # (currently doing nothing here)

    # Calculate error pvalues
    sys.stderr.write('Calculating error pvalues...\n')
    os.system('python ~/projectBenNandita/calculate_error_pvalues.py' + species)
    sys.stderr.write('Done!\n')

    sys.stderr.write('Finished postprocessing MIDAS data!\n')
