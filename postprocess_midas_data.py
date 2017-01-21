### This script runs the necessary post-processing of the MIDAS output so that we can start analyzing
import os
import sys


sys.stderr.write('Calculating species-specific marker gene coverages...\n')
os.system('python calculate_marker_gene_coverage.py')   
sys.stderr.write('Done!\n')

sys.stderr.write('Merging technical replicates from same timepoint...\n')
os.system('python combine_replicates.py')   
sys.stderr.write('Done!\n')

sys.stderr.write('Calculating coverage distributions...\n')
os.system('python calculate_coverage_distribution.py')
sys.stderr.write('Done!\n')

# Calculate core genome set 
# (currently doing nothing here)

# Calculate error pvalues
sys.stderr.write('Calculating error pvalues...\n')
os.system('python calculate_error_pvalues.py')
sys.stderr.write('Done!\n')

sys.stderr.write('Finished postprocessing MIDAS data!\n')