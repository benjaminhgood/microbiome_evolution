### This script runs the necessary post-processing of the MIDAS output so that we can start analyzing
import os
import sys
import parse_midas_data

if len(sys.argv) > 1:
    if sys.argv[1]=='debug':
        species_names = [parse_midas_data.debug_species_name] 
    else:
        species_names = sys.argv[1:]
else:
    species_names = parse_midas_data.parse_good_species_list()

sys.stderr.write("Running postprocessing routines on %d species...\n\n" % len(species_names))

for species_name in species_names:

    # Note: when we want to parallelize this portion of the pipeline
    # We will put the following steps in a wrapper script that can either be
    # sent to QSUB or run on the command line (depending on whether we run 
    # the script in cluster mode or not)
    
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

sys.stderr.write('Finished postprocessing MIDAS data!\n')
