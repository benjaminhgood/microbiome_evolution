import matplotlib 
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import os
species_name=sys.argv[1]   


os.system('python ~/projectBenNandita/print_distance_matrix.py ' + species_name)    

os.system('Rscript ~/projectBenNandita/plot_tree.R ' + species_name)
