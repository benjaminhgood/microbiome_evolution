#!/usr/bin/env python 
### This script runs the necessary post-processing of the MIDAS output so that we can start analyzing
import os
import sys
import parse_midas_data

########################################################################################
#
# Standard header to read in argument information
#
########################################################################################
if len(sys.argv)>1:
    if len(sys.argv)>2:
        debug=True # debug does nothing in this script
        species_name=sys.argv[2]
    else:
        debug=False
        species_name=sys.argv[1]
else:
    sys.stderr.write("Usage: python postprocess_midas_data.py [debug] species_name")
########################################################################################


sys.stderr.write('Postprocessing species: %s\n' % species_name)

sys.stderr.write('Calculating species-specific marker gene coverages...\n')
os.system('python %scalculate_marker_gene_coverage.py %s' % (parse_midas_data.scripts_directory, species_name))   
sys.stderr.write('Done calculating species-specific marker gene coverages!\n')

sys.stderr.write('Calculating coverage distributions...\n')
os.system('python %scalculate_coverage_distribution.py %s' % (parse_midas_data.scripts_directory, species_name))
sys.stderr.write('Done calculating coverage distribution!\n')

# Calculate core genome set 
# (currently doing nothing here)

# Calculate error pvalues
sys.stderr.write('Calculating error pvalues...\n')
os.system('python %scalculate_error_pvalues.py %s' % (parse_midas_data.scripts_directory, species_name))
sys.stderr.write('Done calculating error pvalues!\n')
    
sys.stderr.write("Done postprocessing %s!\n\n" % species_name)
