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


#  sum  all gene changes  across species
all_data['all_species']={}
#for gene_type in ['gene_changes','gene_changes_category']:
for gene_type in ['gene_changes']:
    all_data['all_species'][gene_type]={}
    #all_data['all_species']['gene_changes_category']={}

for species_name in all_data.keys():
    print species_name
    if species_name != 'all_species':
#        for gene_type in ['gene_changes','gene_changes_category']:
        for gene_type in ['gene_changes']:
            for gene in all_data[species_name][gene_type].keys():
                if gene not in all_data['all_species'][gene_type].keys():
                    all_data['all_species'][gene_type][gene] = {'all':[0,0,0,0],'losses':[0,0,0,0],'gains':[0,0,0,0]} #[obs, exp_between, exp_present, exp_pangeome]
                for change_type in ['gains','losses','all']:  
                    all_data['all_species'][gene_type][gene][change_type][0] +=all_data[species_name][gene_type][gene][change_type]
                                                                 

# compute a null for all species by combining all trials across species:
all_data['all_species']['null']={}
#for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes', 'between_host_category','present_category','pangenome_category']:

for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
    all_data['all_species']['null'][null_type]={}

for species_name in all_data.keys(): 
    print species_name
    if species_name != 'all_species' : #CHANGE THIS
        for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
            for gene in all_data[species_name]['null'][null_type].keys():
                if gene not in  all_data['all_species']['null'][null_type].keys():
                    all_data['all_species']['null'][null_type][gene]={'all':[0]*num_trials,'gains':[0]*num_trials,'losses':[0]*num_trials}
                for change_type in ['gains','losses','all']:
                    for i in range(0, num_trials):
                        #print change_type
                        #print all_data['all_species']['null'][null_type][gene][change_type]
                        #print all_data[species_name]['null'][null_type][gene][change_type]
                        all_data['all_species']['null'][null_type][gene][change_type][i]+=all_data[species_name]['null'][null_type][gene][change_type][i]

                

# compute the expected number of changes under different nulls by averaging
for gene in all_data['all_species']['gene_changes'].keys():
    for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
        for change_type in ['gains','losses','all']:
            if gene in all_data['all_species']['null'][null_type].keys():
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

'''
# repeat for gene_category

for gene in all_data['all_species']['gene_changes_category'].keys():
    for null_type in ['between_host_category','present_category','pangenome_category']:
        if gene in all_data['all_species']['null'][null_type].keys():
            tmp_null = numpy.array(all_data['all_species']['null'][null_type][gene])
            prob=sum(all_data['all_species']['gene_changes_category'][gene][0] <= tmp_null)/float(num_trials)
        else:
            prob=0
        if null_type=='between_host_category':
            all_data['all_species']['gene_changes_category'][gene][1]=prob
        elif null_type=='present_category':
            all_data['all_species']['gene_changes_category'][gene][2]=prob
        else:
            all_data['all_species']['gene_changes_category'][gene][3]=prob
        
'''

# print the observed vs expected values of the genes in the all species gene changes:
outFile_gene_change=open('%sgene_changes_accross_species.txt' %  parse_midas_data.analysis_directory,'w')

# store data for plotting cdf: (this says p-val, but actually expectations are being stored)
p_val_arrays={}
for change_type in ['all','gains','losses']: 
    p_val_arrays[change_type]={'between':[],'present':[],'pangenome':[]}

gene_order_sorted = sorted(all_data['all_species']['gene_changes'].items(), key=operator.itemgetter(1))
for i in range(0, len(gene_order_sorted)):
    gene=gene_order_sorted[i][0]
    string=gene 
    for change_type in ['all','gains','losses']:
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


'''
p_val_between_array=numpy.asarray(sorted(p_val_between_array))
p_val_present_array=numpy.asarray(sorted(p_val_present_array))
p_val_pangenome_array=numpy.asarray(sorted(p_val_pangenome_array))

# print the p-values of the genes in the all species category:
outFile_category=open('%sgene_changes_accross_species_category.txt' %  parse_midas_data.analysis_directory,'w')

gene_order_sorted = sorted(all_data['all_species']['gene_changes_category'].items(), key=operator.itemgetter(1))
for i in range(0, len(gene_order_sorted)):
    gene=gene_order_sorted[i][0]
    counts=all_data['all_species']['gene_changes_category'][gene][0]
    p_val_between=all_data['all_species']['gene_changes_category'][gene][1]
    p_val_present=all_data['all_species']['gene_changes_category'][gene][2]
    p_val_pangenome=all_data['all_species']['gene_changes_category'][gene][3]
    print gene + '\t' + str(counts) + '\t' +str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) 
    outFile_category.write(gene + '\t' +str(counts) +'\t'+ str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) + '\n' )

outFile_category.close()
'''

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
#keywords['other']=['other']

# since this is a greedy algorithm, order the more important keywords first
keyword_order=['ABC transporter','phage','transposon','mobilization','integrase', 'plasmid','recombinase','tRNA','ATP','excisionase','transmembrane','replication','regulator','transcription','toxin','restriction','replication','transferase','reductase','phosphatase','helicase','kinase','dehydrogenase','drug','cell wall','primase','resistance','hydrolase','topoisomerase','hypothetical']

common_genes={}

# key=keyword
# value={}
# num={}
# genes={}
for keyword in keywords.keys():
    common_genes[keyword]={'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}


for gene in all_data['all_species']['gene_changes']:
  if all_data['all_species']['gene_changes'][gene]['all'][0] > 0: 
    keyword_found=False
    for keyword in keyword_order:
        for regexp in keywords[keyword]:
            if regexp in gene and keyword_found==False:
                for change_type in ['all','gains','losses']:
                    for i in range(0,4):
                        common_genes[keyword][change_type][i]+=all_data['all_species']['gene_changes'][gene][change_type][i]
                common_genes[keyword]['genes'].append(gene)
                keyword_found=True
    if keyword_found==False and gene !='':
        if gene not in common_genes.keys():
            common_genes[gene]={'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}
        for change_type in ['all','gains','losses']:
            for i in range(0,4):
                common_genes[gene][change_type][i]+=all_data['all_species']['gene_changes'][gene][change_type][i]
        common_genes[gene]['genes'].append(gene)
#    if keyword_found==False:
#        for change_type in ['all','gains','losses']:
#            for i in range(0,4):
#                common_genes['other'][change_type][i]+=all_data['all_species']['gene_changes'][gene][change_type][i]
#        common_genes['other']['genes'].append(gene)

#sorted list by total num

genes_sorted={}
for gene in common_genes.keys():
    genes_sorted[gene]=common_genes[gene]['all'][0]

sorted_genes = sorted(genes_sorted.items(), key=operator.itemgetter(1), reverse=True)


outFile_keywords=open('%sgene_changes_accross_species_keywords.txt' %  parse_midas_data.analysis_directory,'w')

outFile_keywords.write('keyword\tnum_all\texp_all_between\texp_all_present\texp_all_pangenome\tnum_gains\texp_gains_between\texp_gains_present\texp_gains_pangenome\tnum_loss\texp_loss_between\texp_loss_present\texp_loss_pangenome\tgene_names\n')

for i in range (0, len(sorted_genes)):
    gene=sorted_genes[i][0]    
    string=gene
    for change_type in ['all','gains','losses']: 
        for i in range(0,4):
            string += '\t' + str(common_genes[gene][change_type][i])
    string += '\t' + ';'.join(common_genes[gene]['genes'])
    print string
    outFile_keywords.write(string+'\n')

outFile_keywords.close()


###################################################
# plot CDF of p-values for the different nulls 
###################################################

pylab.figure(figsize=(6,6))
prevalence_axis = pylab.subplot(111)

prevalence_axis.set_ylabel('Fraction genes $\leq p$',labelpad=2)
prevalence_axis.set_xlabel('Expected number, $p$',labelpad=2)
prevalence_axis.set_xlim([0,20])
#prevalence_axis.set_ylim([0,1.1])

prevalence_axis.spines['top'].set_visible(False)
prevalence_axis.spines['right'].set_visible(False)
prevalence_axis.get_xaxis().tick_bottom()
prevalence_axis.get_yaxis().tick_left()

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.asarray(p_val_arrays['losses']['present']))
prevalence_axis.step(xs,1-ns*1.0/ns[0],'b-',label='Within-host present genes',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.asarray(p_val_arrays['losses']['pangenome']))
prevalence_axis.step(xs,1-ns*1.0/ns[0],'r-',label='Between-host gene changes',zorder=1)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.asarray(p_val_arrays['losses']['between']))
prevalence_axis.step(xs,1-ns*1.0/ns[0],'k-',label='Pangenome',zorder=0)

#prevalence_axis.set_ylim([0,1.1])
#prevalence_axis.set_xlim([0,1.05])

prevalence_axis.legend(loc='upper right',frameon=False,fontsize=4)

pylab.savefig('%s/expected_gene_change_annotation_losses.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)




'''


##############################################################
# 2D matrix of the p-values (genes x species) for a heatmap:
##############################################################
num_species=len(all_data.keys())
num_genes=len(all_data['all_species']['gene_changes'].keys())-1 # remove hypothetical protein
gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# sort genes based on how abundant they are accross species:
gene_order_sorted = sorted(all_data['all_species']['gene_changes'].items(), key=operator.itemgetter(1))
gene_order=[]
for i in range(0, len(gene_order_sorted)):
    if gene_order_sorted[i][0] !='hypothetical protein':
        gene_order.append(gene_order_sorted[i][0])

species_order=species_phylogeny_utils.sort_phylogenetically(all_data.keys())
species_order.append('all_species')



# fill matrix with probability under the present gene null 
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes']:
                if all_data[species_name]['gene_changes'][gene][2] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx]=2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx]=1
            gene_idx +=1
    species_idx +=1



pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_present_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)



gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# fill matrix with probability under the between-host gene changes null 
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes']:
                if all_data[species_name]['gene_changes'][gene][1] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx]=2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx]=1
            gene_idx +=1
    species_idx +=1



pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_between_host_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)


#####
# repeat for gene_changes_category

# put everything in a 2D matrix so that I can plot a heatmap:
num_species=len(all_data.keys())
num_genes=len(all_data['all_species']['gene_changes_category'].keys())-1 # remove hypothetical protein
gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# sort genes based on how abundant they are accross species:
gene_order_sorted = sorted(all_data['all_species']['gene_changes_category'].items(), key=operator.itemgetter(1))
gene_order=[]
for i in range(0, len(gene_order_sorted)):
    if gene_order_sorted[i][0] !='hypothetical protein':
        gene_order.append(gene_order_sorted[i][0])

species_order=species_phylogeny_utils.sort_phylogenetically(all_data.keys())
species_order.append('all_species')



# fill based on the probability under the present gene model
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes_category']:
                if all_data[species_name]['gene_changes_category'][gene][2] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx] = 2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx] = 1
            gene_idx +=1
    species_idx +=1



pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_category_present_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)


# fill based on the probability under the between-host model
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes_category']:
                if all_data[species_name]['gene_changes_category'][gene][1] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx] = 2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx] = 1
            gene_idx +=1
    species_idx +=1



pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_category_between_host_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)


'''
