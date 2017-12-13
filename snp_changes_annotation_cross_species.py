# Within-host gene change annotation

import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import os
import parse_HMP_data


import pylab
import sys
import numpy
import random

import diversity_utils
import gene_diversity_utils
import sfs_utils
import calculate_substitution_rates
import calculate_temporal_changes
import parse_patric
import species_phylogeny_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice
import operator
import pickle
import seaborn as sns; sns.set()

mpl.rcParams['font.size'] = 5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'



#######################################
# code for loading cross-species data #
#######################################
num_trials=100

good_species_list = parse_midas_data.parse_good_species_list() 


all_data={}
for species_name in good_species_list:
    print species_name
    if os.path.isfile("/pollard/home/ngarud/tmp_intermediate_files/%s_gene_changes.p" %species_name):
        all_data_species = pickle.load( open( "/pollard/home/ngarud/tmp_intermediate_files/%s_gene_changes.p" %species_name, "rb" ))
        if (len(all_data_species[species_name].keys()) >0): # check if there were any gene changes to be outputted. 
            all_data[species_name]=all_data_species[species_name]

#for gene in all_data[species_name]['gene_changes'].keys():
#    print gene + str(all_data[species_name]['gene_changes'][gene]['snps']) + '\t' + str(all_data[species_name]['gene_changes'][gene]['all']) + '\t' + str(all_data[species_name]['gene_changes'][gene]['gains']) + '\t' + str(all_data[species_name]['gene_changes'][gene]['losses'])

#  sum  all gene changes  across species
all_data['all_species']={}
#for gene_type in ['gene_changes','gene_changes_category']:
for gene_type in ['gene_changes']:
    all_data['all_species'][gene_type]={}
    #all_data['all_species']['gene_changes_category']={}

for species_name in all_data.keys():
    if species_name != 'all_species':
#        for gene_type in ['gene_changes','gene_changes_category']:
        for gene_type in ['gene_changes']:
            for gene in all_data[species_name][gene_type].keys():
                if gene not in all_data['all_species'][gene_type].keys():
                    all_data['all_species'][gene_type][gene] = {'snps':[0,0,0,0]} #[obs, exp_between, exp_present, exp_pangeome]
                for change_type in ['snps']:  
                    all_data['all_species'][gene_type][gene][change_type][0] +=all_data[species_name][gene_type][gene][change_type]
                                                                 

# compute a null for all species by combining all trials across species:
all_data['all_species']['null']={}
#for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes', 'between_host_category','present_category','pangenome_category']:

for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
    all_data['all_species']['null'][null_type]={}

for species_name in all_data.keys(): 
    print species_name
    if species_name != 'all_species' and species_name !='Escherichia_coli_58110':
        for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
            for gene in all_data[species_name]['null'][null_type].keys():
                if gene not in  all_data['all_species']['null'][null_type].keys():
                    all_data['all_species']['null'][null_type][gene]={'snps':[0]*num_trials}
                for change_type in ['snps']:
                    for i in range(0, num_trials):
                        all_data['all_species']['null'][null_type][gene][change_type][i]+=all_data[species_name]['null'][null_type][gene][change_type][i]

                

# compute the expected number of changes under different nulls by averaging
for gene in all_data['all_species']['gene_changes'].keys():
    for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
        #for change_type in ['gains','losses','all']:
        for change_type in ['snps']: 
            if gene in all_data['all_species']['null'][null_type].keys():
                print all_data['all_species']['null'][null_type][gene].keys()
                expectation = numpy.array(all_data['all_species']['null'][null_type][gene][change_type]).mean()
                #prob=sum(all_data['all_species']['gene_changes'][gene][change_type] <= tmp_null)/float(num_trials)
            else:
                expectation=0
            if null_type=='between_host_genes':
                all_data['all_species']['gene_changes'][gene][change_type][1]=expectation
            elif null_type=='present_genes':
                all_data['all_species']['gene_changes'][gene][change_type][2]=expectation
            else:
                all_data['all_species']['gene_changes'][gene][change_type][3]=expectation   

# print the observed vs expected values of the genes in the all species gene changes:
outFile_gene_change=open('%ssnp_changes_accross_species.txt' %  parse_midas_data.analysis_directory,'w')

# store data for plotting cdf: (this says p-val, but actually expectations are being stored)
p_val_arrays={}
for change_type in ['snps']: 
    p_val_arrays[change_type]={'between':[],'present':[],'pangenome':[]}

gene_order_sorted = sorted(all_data['all_species']['gene_changes'].items(), key=operator.itemgetter(1))
for i in range(0, len(gene_order_sorted)):
    gene=gene_order_sorted[i][0]
    string=gene 
    for change_type in ['snps']:
        counts=all_data['all_species']['gene_changes'][gene][change_type][0]
        p_val_between=all_data['all_species']['gene_changes'][gene][change_type][1]
        p_val_arrays[change_type]['between'].append(p_val_between)
        p_val_present=all_data['all_species']['gene_changes'][gene][change_type][2]
        p_val_arrays[change_type]['present'].append(p_val_present)
        p_val_pangenome=all_data['all_species']['gene_changes'][gene][change_type][3]
        p_val_arrays[change_type]['pangenome'].append(p_val_pangenome)
        string += '\t' + str(counts) +'\t' +str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) 
    print string
    outFile_gene_change.write(string + '\n' )


outFile_gene_change.close()


##################################################
# list of common genes that show up              #
##################################################

keywords={}
keywords['ABC transporter']=['ABC']
keywords['phage']=['hage']
keywords['transposon']=['onjugati','anspos']
keywords['mobilization']=['mob','obilization','obile']
keywords['integrase']=['ntegrase']
keywords['plasmid']=['plasmid']
keywords['recombinase']=['ecombinase']
keywords['tRNA']=['tRNA']
keywords['ATP']=['ATP']
keywords['excisionase']=['xcisionase']
keywords['transmembrane']=['embrane']
keywords['replication']=['eplication']
keywords['regulator']=['egulator']
keywords['transcription']=['anscription']
keywords['toxin']=['toxin']
keywords['restriction']=['estriction']
keywords['replication']=['eplication']
keywords['transferase']=['ansferase']
keywords['reductase']=['eductase']
keywords['phosphatase']=['phosphatase']
keywords['helicase']=['elicase']
keywords['kinase']=['kinase']
keywords['dehydrogenase']=['dehydrogenase']
keywords['drug']=['drug']
keywords['cell wall']=['ell wall']
keywords['primase']=['imase']
keywords['resistance']=['resistance']
keywords['hydrolase']=['ydrolase']
keywords['topoisomerase']=['opoisomerase']
keywords['hypothetical'] = ['ypothetical']

# since this is a greedy algorithm, order the more important keywords first
keyword_order=['ABC transporter','phage','transposon','mobilization','integrase', 'plasmid','recombinase','tRNA','ATP','excisionase','transmembrane','replication','regulator','transcription','toxin','restriction','replication','transferase','reductase','phosphatase','helicase','kinase','dehydrogenase','drug','cell wall','primase','topoisomerase','hypothetical']

common_genes={}

# key=keyword
# value={}
# num={}
# genes={}
for keyword in keywords.keys():
    #common_genes[keyword]={'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}
    common_genes[keyword]={'snps':[0,0,0,0], 'genes':[]}


for gene in all_data['all_species']['gene_changes']:
  if all_data['all_species']['gene_changes'][gene]['snps'][0] > 0: 
    keyword_found=False
    for keyword in keyword_order:
        for regexp in keywords[keyword]:
            if regexp in gene and keyword_found==False:
                for change_type in ['snps']:
                    for i in range(0,4):
                        common_genes[keyword][change_type][i]+=all_data['all_species']['gene_changes'][gene][change_type][i]
                common_genes[keyword]['genes'].append(gene)
                keyword_found=True
    if keyword_found==False and gene !='':
        if gene not in common_genes.keys():
            common_genes[gene]={'snps':[0,0,0,0], 'genes':[]}
        for change_type in ['snps']:
            for i in range(0,4):
                common_genes[gene][change_type][i]+=all_data['all_species']['gene_changes'][gene][change_type][i]
        common_genes[gene]['genes'].append(gene)



#    if keyword_found==False:
#        for change_type in ['snps']:
#            for i in range(0,4):
#                common_genes['other'][change_type][i]+=all_data['all_species']['gene_changes'][gene][change_type][i]
#        common_genes['other']['genes'].append(gene)

#sorted list by total num
genes_sorted={}
for gene in common_genes.keys():
    genes_sorted[gene]=common_genes[gene]['snps'][0]

sorted_genes = sorted(genes_sorted.items(), key=operator.itemgetter(1), reverse=True)


outFile_keywords=open('%ssnp_changes_accross_species_keywords.txt' %  parse_midas_data.analysis_directory,'w')

#outFile_keywords.write('keyword\tnum_snps\texp_snps_between\texp_snps_present\texp_snps_pangenome\tgene_names\n')
outFile_keywords.write('keyword\tnum_snps\texp_snps_between\texp_snps_present\tgene_names\n')

print common_genes
for i in range (0, len(sorted_genes)):
    gene=sorted_genes[i][0]    
    string=gene
    if common_genes[gene][change_type][0] >0:
        for change_type in ['snps']: 
            #for i in range(0,4):
            for i in range(0,3):
                string += '\t' + str(common_genes[gene][change_type][i])
        string += '\t' + ';'.join(common_genes[gene]['genes'])
        print string
        outFile_keywords.write(string+'\n')

outFile_keywords.close()