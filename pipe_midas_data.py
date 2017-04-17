import sys
import os
import bz2
import parse_midas_data

if len(sys.argv) > 1:
    species_name=sys.argv[1]
else:
    species_name=parse_midas_data.debug_species_name

parse_midas_data.pipe_snps(species_name)
    