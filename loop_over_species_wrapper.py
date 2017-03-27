#!/usr/bin/env python
### This script is used to run a subscript over a list of species. 
### Meant as a general purpose wrapper for species looping functionality.
### Maybe we'll add a cluster specific version later. 
import os
import sys
import parse_midas_data

if len(sys.argv) < 3:
    sys.stderr.write("Usage: python loop_over_species_wrapper.py all|debug|species command...\n")
    sys.exit(1)

# First argument is either 'all', 'debug', or a species name

debug_flag = ""
if sys.argv[1]=='debug':
    species_names = [parse_midas_data.debug_species_name] 
    debug_flag = "--debug"
elif sys.argv[1]=='all':
    species_names = parse_midas_data.parse_good_species_list()
else:
    good_species_names = parse_midas_data.parse_good_species_list()
    species_names = []
    pattern = sys.argv[1]
    for species_name in good_species_names:
        if species_name.startswith(pattern):
            species_names.append(species_name)

# Remaining arguments are command to run, with species name appended as last argument
command = " ".join(sys.argv[2:])

sys.stderr.write("Running command: %s\n" % command)
sys.stderr.write("for %d species...\n\n" % len(species_names))

for species_name in species_names:

    sys.stderr.write('Processing species: %s\n' % species_name)
    os.system('%s %s %s' % (command, debug_flag, species_name))
    sys.stderr.write("Finished processing %s!\n\n" % species_name)

sys.stderr.write('Finished looping over species!\n')
