### This script is used to run a subscript over a list of species. 
### Meant as a general purpose wrapper for species looping functionality.
### Maybe we'll add a cluster specific version later. 
import os
import sys
import parse_midas_data

if len(sys.argv) < 2:
    sys.stderr.write("Usage: python print_good_species_list.py [num|names]\n")
    sys.exit(1)
    
species_names = parse_midas_data.parse_good_species_list()    
if sys.argv[1]=='num':
    print len(species_names)
else:
    for species_name in species_names:
        print species_name
