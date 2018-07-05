import os.path
import bz2

# read in the PATRIC genomes
species_list=['B_vulgatus', 'B_fragilis', 'P_distasonis']

PATRIC_db={} # key=genome_id, value=[name,abbreviation]

for species in species_list:
    PATRIC_inFN='/pollard/home/ngarud/BenNanditaProject/fake_data/genome_lists/PATRIC_genome_%s.txt' %species
    PATRIC_infile=open(PATRIC_inFN, 'r')
    for line in PATRIC_infile:
        items=line.strip('\n').split('\t')
        genome_id = items[0][1:]
        genome_id = genome_id.replace('"', '')
        genome_name= items[1]
        genome_name = genome_name.replace('"', '')
        PATRIC_db[genome_id]=[genome_name,species]


# read in the species abundance files
infile=bz2.BZ2File(os.path.expanduser('~/ben_nandita_hmp_data/simulations/species/relative_abundance.txt.bz2'),"r")

header=infile.readline().strip('\n').split('\t')

species_nonzero=['Bacteroides_fragilis_54507', 'Bacteroides_vulgatus_57955', 'Parabacteroides_distasonis_56985', 'Bacteroides_fragilis_56548', 'Bacteroides_ovatus_58035', 'Bacteroides_xylanisolvens_57185', 'Citrobacter_freundii_56776', 'Citrobacter_freundii_56148'] # this is the list of species we want to know the relative abundance of in our data

# store the abundances in this dictionary:
abundances={}
for sample in header[1:]:
    abundances[sample]={}
    for species in species_nonzero:
        abundances[sample][species]=0
    # also parse the sample:
    if '_' in sample: # mixture
        substrs=sample.split('_')
        abundances[sample]['genome1']=PATRIC_db[substrs[0]][0]
        abundances[sample]['genome2']=PATRIC_db[substrs[1]][0]
        abundances[sample]['abbreviation1']=PATRIC_db[substrs[0]][1]
        abundances[sample]['abbreviation2']=PATRIC_db[substrs[1]][1]
        abundances[sample]['iteration']=substrs[2]
        abundances[sample]['coverage']=substrs[3]
    else: # isolate
        genome_name=PATRIC_db[sample][0]
        abbreviation=PATRIC_db[sample][1]
        abundances[sample]['genome1']=genome_name
        abundances[sample]['genome2']='-'
        abundances[sample]['abbreviation1']=abbreviation
        abundances[sample]['abbreviation2']='-'
        abundances[sample]['iteration']='1'
        abundances[sample]['coverage']='100'


for line in infile:
    items=line.strip('\n').split('\t')
    species=items[0]
    if species in species_nonzero:
        for i in range(1, len(items)):
            sample=header[i]
            abundance=items[i]
            abundances[sample][species]=abundance

outFile=open(os.path.expanduser('~/ben_nandita_hmp_analysis/simulations_species_abundances.txt'),'w')


# output a table:
header='sample\tgenome1\tgenome2\tabbreviation1\tabbreviation2\titeration\tcoverage\t'

for species in species_nonzero:
    header += species +'\t'

outFile.write(header +'\n')

for sample in abundances:
        genome1=abundances[sample]['genome1']
        genome2=abundances[sample]['genome2']
        abbreviation1=abundances[sample]['abbreviation1']
        abbreviation2=abundances[sample]['abbreviation2']
        iteration=abundances[sample]['iteration']
        coverage=abundances[sample]['coverage']
        s= sample + '\t' + genome1 + '\t' + genome2 + '\t' + abbreviation1 + '\t' + abbreviation2 + '\t' + iteration + '\t' + coverage +'\t' 
        for species in species_nonzero:
            s+= abundances[sample][species] +'\t'
        outFile.write(s + '\n')


outFile.close()
