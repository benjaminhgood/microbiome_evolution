import parse_midas_data
import sys

file = open(parse_midas_data.scripts_directory+"manual_clade_definitions.txt","r")
    
line = file.readline().strip()
    
# Just put some default values there in case of issue
current_species = ""
current_clade = 1

output_strs = []
# Header
output_strs.append("\t".join(["Species", "Sample", "Clade"]))
    
while line!="":
    
    items = line.split()
        
    if items[0].isdigit():
        # is a sample entry
        sample_name = items[1].strip()
        output_strs.append("\t".join([current_species, sample_name, str(current_clade)]))
        
    elif items[0].startswith('-'):
        # delimits a clade
        current_clade+=1
    else:
        # new species
        current_species = items[0].strip()
        current_clade = 1
    
    line = file.readline().strip()
    
file.close()

output_file = open('%s/clade_assignments.txt' % (parse_midas_data.analysis_directory),"w")
output_file.write("\n".join(output_strs))
output_file.close()