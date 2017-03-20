import parse_midas_data
import sys
import numpy
import os

species_name = sys.argv[1]

parse_midas_data.pipe_snps(species_name,combination_type=None,avg_depth_threshold=10)
    
