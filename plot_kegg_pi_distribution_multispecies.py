import matplotlib  
matplotlib.use('Agg') 
import os
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
from numpy.random import choice
import matplotlib.cm as cmx
import matplotlib.colors as colors
import pickle
import pandas
import seaborn as sns
import pandas
import random

# plotting tools
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

#############

#species_names=['Alistipes_putredinis_61533','Bacteroides_uniformis_57318']
species_names = parse_midas_data.parse_good_species_list()

files={}

for species_name in  species_names:
    if (os.path.exists(os.path.expanduser('~/tmp_intermediate_files/kegg_pi_%s.dat' % species_name))):
        files[species_name]=pickle.load(open(os.path.expanduser('~/tmp_intermediate_files/kegg_pi_%s.dat' % species_name),'rb'))
    else:
        print species_name

colors=['#a1d99b','#c994c7']  

####################################################################
# Order the species according to their phylogenetic relationships. #
#################################################################### 

# to make this plot I will have to feed the list of species to R. When I am in R, I will drop the species that are not in the list. Then, R will output a tree. 

# output a file to the temp dir with the species number only
outFN=os.path.expanduser("~/tmp_intermediate_files/species_names.txt")
out_file=open(outFN, 'w')
for species in species_names:
    print species
    id_no=species.strip().split('_')[2]
    out_file.write(species + '\t' + id_no +'\n')
    
# run the R script to prune the species tree
os.system('Rscript ~/ben_nandita_hmp_scripts/prune_species_tree.R ' +  outFN)

# read in the order of the species that are on the tree:
file=open(os.path.expanduser("~/tmp_intermediate_files/species_order.txt"))
species_order=[]
for line in file:
    species_order.append(line.strip())


######################################################################
# plot a comparison of pi for core vs variable genes accross species
#####################################################################
# try ordering species based on median piS in core
species_list=[]
mean_avg_pi=[]
sample_size=[]
for species_name in species_names:
    if species_name in files.keys():
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        avg_pi_per_pathway_core=files[species_name]['avg_pi_per_pathway_core']
        avg_pi_per_pathway_variable=files[species_name]['avg_pi_per_pathway_variable']
        passed_sites_per_pathway_core=files[species_name]['passed_sites_per_pathway_core']
        passed_sites_per_pathway_variable=files[species_name]['passed_sites_per_pathway_variable'] #need to fix this
        num_genes_per_pathway_core=files[species_name]['num_genes_per_pathway_core']
        num_genes_per_pathway_variable=files[species_name]['num_genes_per_pathway_variable']
        num_people_with_data_pathway_core=files[species_name]['num_people_with_data_pathway_core']
        num_people_with_data_pathway_variable=files[species_name]['num_people_with_data_pathway_variable']
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core'] 
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable'] 
        fraction_nonsynonymous_per_pathway_core=files[species_name]['fraction_nonsynonymous_per_pathway_core']
        fraction_nonsynonymous_per_pathway_variable=files[species_name]['fraction_nonsynonymous_per_pathway_variable']
        fixation_opportunities_per_pathway_syn_non_core=files[species_name]['fixation_opportunities_per_pathway_syn_non_core']
        fixation_opportunities_per_pathway_syn_non_variable=files[species_name]['fixation_opportunities_per_pathway_syn_non_variable']
        num_genes_per_pathway_syn_non_core=files[species_name]['num_genes_per_pathway_syn_non_core']
        num_genes_per_pathway_syn_non_variable=files[species_name]['num_genes_per_pathway_syn_non_variable']
        num_people_with_data_per_pathway_fixations_core=files[species_name]['num_people_with_data_per_pathway_fixations_core']
        num_people_with_data_per_pathway_fixations_variable=files[species_name]['num_people_with_data_per_pathway_fixations_variable']
        dtot_core=files[species_name]['dtot_core'] 
        dtot_variable=files[species_name]['dtot_variable']
        dtot_per_pathway_core=files[species_name]['dtot_per_pathway_core']
        dtot_per_pathway_variable=files[species_name]['dtot_per_pathway_variable']
        fixation_opportunities_per_pathway_all_core=files[species_name]['fixation_opportunities_per_pathway_all_core']
        fixation_opportunities_per_pathway_all_variable=files[species_name]['fixation_opportunities_per_pathway_all_variable']
        num_genes_per_pathway_tot_core=files[species_name]['num_genes_per_pathway_tot_core']
        num_genes_per_pathway_tot_variable=files[species_name]['num_genes_per_pathway_tot_variable']
        num_people_with_data_per_pathway_tot_core=files[species_name]['num_people_with_data_per_pathway_tot_core']
        num_people_with_data_per_pathway_tot_variable=files[species_name]['num_people_with_data_per_pathway_tot_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        avg_pi=numpy.median(avg_pi_matrix_core[same_sample_idxs])
        species_list.append(species_name)
        mean_avg_pi.append(avg_pi)
        sample_size.append(len(same_sample_idxs[0]))

df=pandas.DataFrame({'species':species_list,'sample_size':sample_size})
df=df.sort('sample_size')
sorted_species=list(df['species']) # this has a list of species sorted by sample size 

######################################################################
# plot a comparison of pi for core vs variable genes accross species
#####################################################################

pylab.figure(figsize=(6,20))
pylab.xlabel('Pi/bp')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1e-1)

data=[]    
labels=[]   
for species_name in sorted_species:
    if species_name in files.keys():
#for species_name in sorted_species:
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        data.append(avg_pi_matrix_core[same_sample_idxs]) 
        labels.append(species_name + '_core, m=' + str(len(same_sample_idxs[0]))) 
        data.append(avg_pi_matrix_variable[same_sample_idxs]) 
        labels.append(species_name + '_variable, m=' + str(len(same_sample_idxs[0]))) 

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
k=0  
for patch in bp['boxes']: 
    patch.set_facecolor(colors[k%2])
    k+=1

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   
pylab.axvline(x=0.001, ymin=0, ymax=1, hold=None)

pylab.savefig('%s/core_vs_variable_genes_pi_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)



######################################################################
# plot a comparison of fraction nonsyn for core vs variable genes accross species
#####################################################################

pylab.figure(figsize=(6,20))
pylab.xlabel('Fraction nonsynonymous fixations')
pylab.ylabel("Species/gene type")
pylab.xlim(0,1)

data=[]    
labels=[]   
for species_name in sorted_species:
    if species_name in files.keys():
#for species_name in sorted_species:
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        same_sample_idxs=files[species_name]['same_sample_idxs']
        data.append(fraction_nonsynonymous_core[diff_subject_idxs]) 
        labels.append(species_name + '_core, n=' + str(len(same_sample_idxs[0]))) 
        data.append(fraction_nonsynonymous_variable[diff_subject_idxs]) 
        labels.append(species_name + '_variable, n=' + str(len(same_sample_idxs[0]))) 

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
k=0  
for patch in bp['boxes']: 
    patch.set_facecolor(colors[k%2])
    k+=1

locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/core_vs_variable_genes_fraction_nonsyn_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)




######################################################################
# plot a comparison of fraction nonsyn for core vs variable genes for just A. putredinis and B. vulgatus
#####################################################################

fig=pylab.figure(figsize=(4,6))

pylab.ylabel('Fraction nonsynonymous fixations')
pylab.xlabel("Species/gene type")
pylab.ylim(0,0.4)

data=[]    
labels=[]   
colors=['#de2d26','#3182bd']

ax = fig.add_subplot(111)
ax.patch.set_facecolor('white')
ax.patch.set_alpha(0)
 
for species_name in ['Bacteroides_vulgatus_57955', 'Alistipes_putredinis_61533']:
    if species_name in files.keys():
#for species_name in sorted_species:
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        same_sample_idxs=files[species_name]['same_sample_idxs']
        data.append(fraction_nonsynonymous_core[diff_subject_idxs]) 
        labels.append(species_name + '_core, n=' + str(len(same_sample_idxs[0]))) 
        data.append(fraction_nonsynonymous_variable[diff_subject_idxs]) 
        labels.append(species_name + '_variable, n=' + str(len(same_sample_idxs[0]))) 

bp=pylab.boxplot(data,patch_artist=True)
k=0  
for patch in bp['boxes']: 
    patch.set_facecolor(colors[k%2])
    k+=1

locs, dummy_labels = pylab.xticks()  
pylab.xticks(locs, labels, fontsize=9)   

#fig.set_facecolor("white")

fig.savefig('%s/core_vs_variable_genes_fraction_nonsyn_B_vul_A_put.png' % (parse_midas_data.analysis_directory),facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=300)




######################################################################
# plot a comparison of total fixations for core vs variable genes accross species
#####################################################################


pylab.figure(figsize=(6,20))
pylab.xlabel('Total fixations')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1e-1)

data=[]    
labels=[]   
for species_name in sorted_species:
    if species_name in files.keys():
#for species_name in sorted_species:
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        data.append(dtot_core[diff_subject_idxs]) 
        labels.append(species_name + '_core, m=' + str(len(same_sample_idxs[0]))) 
        data.append(dtot_variable[diff_subject_idxs]) 
        labels.append(species_name + '_variable, m=' + str(len(same_sample_idxs[0]))) 

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
k=0  
for patch in bp['boxes']: 
    patch.set_facecolor(colors[k%2])
    k+=1

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/core_vs_variable_genes_fixations_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


################################################################
# Plot fraction dN vs dtot for all species core vs var 
################################################################  


pylab.figure(figsize=(8,len(files.keys()*3)))

plot_no=1
for species_name in species_names:
    if species_name in files.keys():
#for species_name in sorted_species:
        pylab.subplot(len(files.keys()), 1, plot_no)
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        pylab.scatter(dtot_core[diff_subject_idxs],fraction_nonsynonymous_core[diff_subject_idxs], color='red',marker='.', label=species_name)
        pylab.scatter(dtot_variable[diff_subject_idxs],fraction_nonsynonymous_variable[diff_subject_idxs], color='blue',marker='.')
        pylab.title(species_name)
        pylab.ylabel("Fraction nonsynonymous")
        pylab.xlim(1e-7,1e-1)
        pylab.ylim(0,1)
        pylab.xscale('log')
        plot_no+=1

pylab.xlabel('Total fixations')
pylab.legend(loc='lower right',frameon=False)

pylab.savefig('%s/core_vs_variable_genes_fraction_nonsyn_vs_fixations_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)



###################################################
# Bar plot of number of species with each pathway #
###################################################


# iterate through the species
# add the kegg pathway to a dictionary
# How many species share each pathway?
# make a bar plot showing this info

pathways_variable={}
pathways_core={}

for species_name in species_names:
    if species_name in files.keys():
        avg_pi_per_pathway_core=files[species_name]['avg_pi_per_pathway_core']
        avg_pi_per_pathway_variable=files[species_name]['avg_pi_per_pathway_variable']
        passed_sites_per_pathway_core=files[species_name]['passed_sites_per_pathway_core']
        num_genes_per_pathway_core=files[species_name]['num_genes_per_pathway_core']
        num_people_with_data_pathway_core=files[species_name]['num_people_with_data_pathway_core']
        #iterate through the pathways
        for pathway in avg_pi_per_pathway_variable.keys():
            if pathway not in pathways_variable:
                pathways_variable[pathway]=[]
            pathways_variable[pathway].append(species_name)
        for pathway in avg_pi_per_pathway_core.keys():
            if pathway not in pathways_core:
                pathways_core[pathway]=[]
            pathways_core[pathway].append(species_name)

#plot for core pathways
species_counts=[]
pathway_labels=[]
for pathway in pathways_core:
    species_counts.append(len(pathways_core[pathway]))
    pathway_labels.append(pathway)

table=[pathway_labels,species_counts]
df=pandas.DataFrame({'num_species':species_counts,'pathways':pathway_labels})
df_sorted=df.sort(['num_species'], ascending=0)

ypos=numpy.arange(len(df_sorted['num_species']))
pylab.figure(figsize=(12,15)) 
pylab.barh(ypos, df_sorted['num_species']) 
pylab.title("Number of species sharing each pathway (core genes)", fontsize=8)
pylab.yticks(ypos, df_sorted['pathways'], fontsize=9)
pylab.xlabel('Number of species')
pylab.savefig('%s/num_species_per_pathway_core.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# plot for variable pathways
species_counts=[]
pathway_labels=[]
for pathway in pathways_variable:
    species_counts.append(len(pathways_variable[pathway]))
    pathway_labels.append(pathway)

table=[pathway_labels,species_counts]
df=pandas.DataFrame({'num_species':species_counts,'pathways':pathway_labels})
df_sorted=df.sort(['num_species'], ascending=0)

ypos=numpy.arange(len(df_sorted['num_species']))
pylab.figure(figsize=(12,15))
pylab.barh(ypos, df_sorted['num_species']) 
pylab.title("Number of species sharing each pathway (variable genes)", fontsize=8)
pylab.yticks(ypos, df_sorted['pathways'], fontsize=9)
pylab.xlabel('Number of species')
pylab.savefig('%s/num_species_per_pathway_variable.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)



#########################################################################
# plot pi within pathways with the most data, with species side-by-side #
#########################################################################

pylab.figure(figsize=(6,100))
pylab.xlabel('pi/bp')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1e-1)

pathways_core={}
num_genes={}
num_people={}
min_passed_sites_per_person=100
pathways_core['All core genes']=[]
pathways_core['All variable genes']=[]
min_mean_number_genes=10 # this is the min mean number of genes that need to be in a pathway in order for it to be plotted. 

for species_name in species_names:
    if species_name in files.keys():
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        avg_pi_per_pathway_core=files[species_name]['avg_pi_per_pathway_core']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        passed_sites_per_pathway_core=files[species_name]['passed_sites_per_pathway_core']
        num_genes_per_pathway_core=files[species_name]['num_genes_per_pathway_core']
        num_people_with_data_pathway_core=files[species_name]['num_people_with_data_pathway_core']
        # add the full genome data:
        pathways_core['All core genes'].append(avg_pi_matrix_core[same_sample_idxs])
        pathways_core['All variable genes'].append(avg_pi_matrix_variable[same_sample_idxs])
        for pathway in avg_pi_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(passed_sites_per_pathway_core[pathway][same_sample_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if pathway not in pathways_core:
                    pathways_core[pathway]=[]
                    num_genes[pathway]=[]
                    num_people[pathway]=[]
                pathways_core[pathway].append(avg_pi_per_pathway_core[pathway][same_sample_idxs][high_num_sites_idx]) 
                num_genes[pathway].append(num_genes_per_pathway_core[pathway])
                num_people[pathway].append(num_people_with_data_pathway_core[pathway])

# iterate through all pathways and concatenate into a single list
data=[]    
labels=[]   
color_list=[]
k=0
for pathway in ['All core genes','All variable genes','','Annotated pathways']:
    for i in range(0, len(pathways_core[pathway])):
        data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
        color_list.append(colors[k%2])
        if i ==0:
            if pathway =='':
                labels.append('Unannotated pathways')
            else:
                labels.append(pathway) 
        else:
            labels.append('')
    k+=1

k=0
for pathway in pathways_core:
    if pathway !='' and pathway !='Annotated pathways' and pathway !='All core genes' and pathway != 'All variable genes':
        # check if the mean_num_genes is >=min_mean_number_genes:
        mean_num_genes=numpy.mean(numpy.asarray(num_genes[pathway]))
        mean_num_people=numpy.mean(numpy.asarray(num_people[pathway]))
        if mean_num_genes>=min_mean_number_genes:
            for i in range(0, len(pathways_core[pathway])):
                data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
                color_list.append(colors[k%2])
                if i ==0:
                    labels.append(pathway+ ', n=' + str(int(mean_num_genes))+ ', m='+str(int(mean_num_people))) 
                else:
                    labels.append('')
            k+=1

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
i=0
for patch in bp['boxes']: 
    patch.set_facecolor(color_list[i])
    i+=1

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/pi_per_pathway_core_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


#########################################################################
# plot pi within pathways with the most data, with species side-by-side #
# repeat for variable genes
#########################################################################

pylab.figure(figsize=(6,25))
pylab.xlabel('pi/bp')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1e-1)

pathways_core={}
num_genes={}
num_people={}
min_passed_sites_per_person=100
pathways_core['All core genes']=[]
pathways_core['All variable genes']=[]
min_mean_number_genes=5 # this is the min mean number of genes that need to be in a pathway in order for it to be plotted. 

for species_name in species_names:
    if species_name in files.keys():
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        avg_pi_per_pathway_core=files[species_name]['avg_pi_per_pathway_variable']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        passed_sites_per_pathway_core=files[species_name]['passed_sites_per_pathway_variable']
        num_genes_per_pathway_core=files[species_name]['num_genes_per_pathway_variable']
        num_people_with_data_pathway_core=files[species_name]['num_people_with_data_pathway_variable']
        # add the full genome data:
        pathways_core['All core genes'].append(avg_pi_matrix_core[same_sample_idxs])
        pathways_core['All variable genes'].append(avg_pi_matrix_variable[same_sample_idxs])
        for pathway in avg_pi_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(passed_sites_per_pathway_core[pathway][same_sample_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if pathway not in pathways_core:
                    pathways_core[pathway]=[]
                    num_genes[pathway]=[]
                    num_people[pathway]=[]
                pathways_core[pathway].append(avg_pi_per_pathway_core[pathway][same_sample_idxs][high_num_sites_idx]) 
                num_genes[pathway].append(num_genes_per_pathway_core[pathway])
                num_people[pathway].append(num_people_with_data_pathway_core[pathway])

# iterate through all pathways and concatenate into a single list
data=[]    
labels=[]   
color_list=[]
k=0
for pathway in ['All core genes','All variable genes','','Annotated pathways']:
    for i in range(0, len(pathways_core[pathway])):
        data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
        color_list.append(colors[k%2])
        if i ==0:
            if pathway =='':
                labels.append('Unannotated pathways')
            else:
                labels.append(pathway) 
        else:
            labels.append('')
    k+=1

k=0
for pathway in pathways_core:
    if pathway !='' and pathway !='Annotated pathways' and pathway !='All core genes' and pathway != 'All variable genes':
        # check if the mean_num_genes is >=min_mean_number_genes:
        mean_num_genes=numpy.mean(numpy.asarray(num_genes[pathway]))
        mean_num_people=numpy.mean(numpy.asarray(num_people[pathway]))
        if mean_num_genes>=min_mean_number_genes:
            for i in range(0, len(pathways_core[pathway])):
                data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
                color_list.append(colors[k%2])
                if i ==0:
                    labels.append(pathway+ ', n=' + str(int(mean_num_genes))+ ', m='+str(int(mean_num_people))) 
                else:
                    labels.append('')
            k+=1

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
i=0
for patch in bp['boxes']: 
    patch.set_facecolor(color_list[i])
    i+=1

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/pi_per_pathway_variable_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)




################################################################################################
# plot fraction nonsyn fixations within pathways with the most data, with species side-by-side #
################################################################################################

pylab.figure(figsize=(6,100))
pylab.xlabel('Fraction of fixations that are nonsynonymous')
pylab.ylabel("Species/gene type")
pylab.xlim(0,1)

pathways_core={}
num_genes={}
num_people={}
min_passed_sites_per_person=100
pathways_core['All core genes']=[]
pathways_core['All variable genes']=[]
min_mean_number_genes=10 # this is the min mean number of genes that need to be in a pathway in order for it to be plotted. 

for species_name in species_names:
    if species_name in files.keys():
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        fraction_nonsynonymous_per_pathway_core=files[species_name]['fraction_nonsynonymous_per_pathway_core']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_syn_non_core=files[species_name]['fixation_opportunities_per_pathway_syn_non_core']
        num_genes_per_pathway_syn_non_core=files[species_name]['num_genes_per_pathway_syn_non_core']
        num_people_with_data_per_pathway_fixations_core=files[species_name]['num_people_with_data_per_pathway_fixations_core']
        # add the full genome data:
        pathways_core['All core genes'].append(fraction_nonsynonymous_core[diff_subject_idxs])
        pathways_core['All variable genes'].append(fraction_nonsynonymous_variable[diff_subject_idxs])
        for pathway in fraction_nonsynonymous_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if pathway not in pathways_core:
                    pathways_core[pathway]=[]
                    num_genes[pathway]=[]
                    num_people[pathway]=[]
                pathways_core[pathway].append(fraction_nonsynonymous_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                num_genes[pathway].append(num_genes_per_pathway_syn_non_core[pathway])
                num_people[pathway].append(num_people_with_data_per_pathway_fixations_core[pathway])

# iterate through all pathways and concatenate into a single list
data=[]    
labels=[]   
color_list=[]
k=0
for pathway in ['All core genes','All variable genes','','Annotated pathways']:
    for i in range(0, len(pathways_core[pathway])):
        data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
        color_list.append(colors[k%2])
        if i ==0:
            if pathway =='':
                labels.append('Unannotated pathways')
            else:
                labels.append(pathway) 
        else:
            labels.append('')
    k+=1

k=0
for pathway in pathways_core:
    if pathway !='' and pathway !='Annotated pathways' and pathway !='All core genes' and pathway != 'All variable genes':
        # check if the mean_num_genes is >=min_mean_number_genes:
        mean_num_genes=numpy.mean(numpy.asarray(num_genes[pathway]))
        mean_num_people=numpy.mean(numpy.asarray(num_people[pathway]))
        if mean_num_genes>=min_mean_number_genes:
            for i in range(0, len(pathways_core[pathway])):
                data.append(pathways_core[pathway][i])
                color_list.append(colors[k%2])
                if i ==0:
                    labels.append(pathway+ ', n=' + str(int(mean_num_genes))+ ', m='+str(int(mean_num_people))) 
                else:
                    labels.append('')
            k+=1

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
i=0
for patch in bp['boxes']: 
    patch.set_facecolor(color_list[i])
    i+=1

locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/fraction_nonsynonymous_per_pathway_core_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)



################################################################################################
# plot fraction nonsyn fixations within pathways with the most data, with species side-by-side #
# repeat for variable genes
################################################################################################

pylab.figure(figsize=(6,25))
pylab.xlabel('Fraction of fixations that are nonsynonymous')
pylab.ylabel("Species/gene type")
pylab.xlim(0,1)

pathways_core={}
num_genes={}
num_people={}
min_passed_sites_per_person=100
pathways_core['All core genes']=[]
pathways_core['All variable genes']=[]
min_mean_number_genes=5 # this is the min mean number of genes that need to be in a pathway in order for it to be plotted. 

for species_name in species_names:
    if species_name in files.keys():
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        fraction_nonsynonymous_per_pathway_core=files[species_name]['fraction_nonsynonymous_per_pathway_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_syn_non_core=files[species_name]['fixation_opportunities_per_pathway_syn_non_variable']
        num_genes_per_pathway_syn_non_core=files[species_name]['num_genes_per_pathway_syn_non_variable']
        num_people_with_data_per_pathway_fixations_core=files[species_name]['num_people_with_data_per_pathway_fixations_variable']
        # add the full genome data:
        pathways_core['All core genes'].append(fraction_nonsynonymous_core[diff_subject_idxs])
        pathways_core['All variable genes'].append(fraction_nonsynonymous_variable[diff_subject_idxs])
        for pathway in fraction_nonsynonymous_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if pathway not in pathways_core:
                    pathways_core[pathway]=[]
                    num_genes[pathway]=[]
                    num_people[pathway]=[]
                pathways_core[pathway].append(fraction_nonsynonymous_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                num_genes[pathway].append(num_genes_per_pathway_syn_non_core[pathway])
                num_people[pathway].append(num_people_with_data_per_pathway_fixations_core[pathway])

# iterate through all pathways and concatenate into a single list
data=[]    
labels=[]   
color_list=[]
k=0
for pathway in ['All core genes','All variable genes','','Annotated pathways']:
    for i in range(0, len(pathways_core[pathway])):
        data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
        color_list.append(colors[k%2])
        if i ==0:
            if pathway =='':
                labels.append('Unannotated pathways')
            else:
                labels.append(pathway) 
        else:
            labels.append('')
    k+=1

k=0
for pathway in pathways_core:
    if pathway !='' and pathway !='Annotated pathways' and pathway !='All core genes' and pathway != 'All variable genes':
        # check if the mean_num_genes is >=min_mean_number_genes:
        mean_num_genes=numpy.mean(numpy.asarray(num_genes[pathway]))
        mean_num_people=numpy.mean(numpy.asarray(num_people[pathway]))
        if mean_num_genes>=min_mean_number_genes:
            for i in range(0, len(pathways_core[pathway])):
                data.append(pathways_core[pathway][i])
                color_list.append(colors[k%2])
                if i ==0:
                    labels.append(pathway+ ', n=' + str(int(mean_num_genes))+ ', m='+str(int(mean_num_people))) 
                else:
                    labels.append('')
            k+=1

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
i=0
for patch in bp['boxes']: 
    patch.set_facecolor(color_list[i])
    i+=1

locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/fraction_nonsynonymous_per_pathway_variable_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)






################################################################################################
# plot fixations within pathways with the most data, with species side-by-side #
################################################################################################

pylab.figure(figsize=(6,100))
pylab.xlabel('Total divergence')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1)

pathways_core={}
num_genes={}
num_people={}
min_passed_sites_per_person=100
pathways_core['All core genes']=[]
pathways_core['All variable genes']=[]
min_mean_number_genes=10 # this is the min mean number of genes that need to be in a pathway in order for it to be plotted. 

for species_name in species_names:
    if species_name in files.keys():
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        dtot_per_pathway_core=files[species_name]['dtot_per_pathway_core']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_all_core=files[species_name]['fixation_opportunities_per_pathway_all_core']
        num_genes_per_pathway_tot_core=files[species_name]['num_genes_per_pathway_tot_core']
        num_people_with_data_per_pathway_tot_core=files[species_name]['num_people_with_data_per_pathway_tot_core']
        # add the full genome data:
        pathways_core['All core genes'].append(dtot_core[diff_subject_idxs])
        pathways_core['All variable genes'].append(dtot_variable[diff_subject_idxs])
        for pathway in dtot_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_all_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if pathway not in pathways_core:
                    pathways_core[pathway]=[]
                    num_genes[pathway]=[]
                    num_people[pathway]=[]
                pathways_core[pathway].append(dtot_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                num_genes[pathway].append(num_genes_per_pathway_tot_core[pathway])
                num_people[pathway].append(num_people_with_data_per_pathway_tot_core[pathway])

# iterate through all pathways and concatenate into a single list
data=[]    
labels=[]   
color_list=[]
k=0
for pathway in ['All core genes','All variable genes','','Annotated pathways']:
    for i in range(0, len(pathways_core[pathway])):
        data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
        color_list.append(colors[k%2])
        if i ==0:
            if pathway =='':
                labels.append('Unannotated pathways')
            else:
                labels.append(pathway) 
        else:
            labels.append('')
    k+=1

k=0
for pathway in pathways_core:
    if pathway !='' and pathway !='Annotated pathways' and pathway !='All core genes' and pathway != 'All variable genes':
        # check if the mean_num_genes is >=min_mean_number_genes:
        mean_num_genes=numpy.mean(numpy.asarray(num_genes[pathway]))
        mean_num_people=numpy.mean(numpy.asarray(num_people[pathway]))
        if mean_num_genes>=min_mean_number_genes:
            for i in range(0, len(pathways_core[pathway])):
                data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
                color_list.append(colors[k%2])
                if i ==0:
                    labels.append(pathway+ ', n=' + str(int(mean_num_genes))+ ', m='+str(int(mean_num_people))) 
                else:
                    labels.append('')
            k+=1

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
i=0
for patch in bp['boxes']: 
    patch.set_facecolor(color_list[i])
    i+=1


pylab.xscale('log')
locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/fixations_per_pathway_core_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)

################################################################################################
# plot fixations within pathways with the most data, with species side-by-side #
# repeat for variable genes
################################################################################################

pylab.figure(figsize=(6,25))
pylab.xlabel('Total divergence')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1)

pathways_core={}
num_genes={}
num_people={}
min_passed_sites_per_person=100
pathways_core['All core genes']=[]
pathways_core['All variable genes']=[]
min_mean_number_genes=5 # this is the min mean number of genes that need to be in a pathway in order for it to be plotted. 

for species_name in species_names:
    if species_name in files.keys():
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        dtot_per_pathway_core=files[species_name]['dtot_per_pathway_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_all_core=files[species_name]['fixation_opportunities_per_pathway_all_variable']
        num_genes_per_pathway_tot_core=files[species_name]['num_genes_per_pathway_tot_variable']
        num_people_with_data_per_pathway_tot_core=files[species_name]['num_people_with_data_per_pathway_tot_variable']
        # add the full genome data:
        pathways_core['All core genes'].append(dtot_core[diff_subject_idxs])
        pathways_core['All variable genes'].append(dtot_variable[diff_subject_idxs])
        for pathway in dtot_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_all_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if pathway not in pathways_core:
                    pathways_core[pathway]=[]
                    num_genes[pathway]=[]
                    num_people[pathway]=[]
                pathways_core[pathway].append(dtot_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                num_genes[pathway].append(num_genes_per_pathway_tot_core[pathway])
                num_people[pathway].append(num_people_with_data_per_pathway_tot_core[pathway])

# iterate through all pathways and concatenate into a single list
data=[]    
labels=[]   
color_list=[]
k=0
for pathway in ['All core genes','All variable genes','','Annotated pathways']:
    for i in range(0, len(pathways_core[pathway])):
        data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
        color_list.append(colors[k%2])
        if i ==0:
            if pathway =='':
                labels.append('Unannotated pathways')
            else:
                labels.append(pathway) 
        else:
            labels.append('')
    k+=1

k=0
for pathway in pathways_core:
    if pathway !='' and pathway !='Annotated pathways' and pathway !='All core genes' and pathway != 'All variable genes':
        # check if the mean_num_genes is >=min_mean_number_genes:
        mean_num_genes=numpy.mean(numpy.asarray(num_genes[pathway]))
        mean_num_people=numpy.mean(numpy.asarray(num_people[pathway]))
        if mean_num_genes>=min_mean_number_genes:
            for i in range(0, len(pathways_core[pathway])):
                data.append(numpy.clip(pathways_core[pathway][i],5*1e-7,1))
                color_list.append(colors[k%2])
                if i ==0:
                    labels.append(pathway+ ', n=' + str(int(mean_num_genes))+ ', m='+str(int(mean_num_people))) 
                else:
                    labels.append('')
            k+=1

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
i=0
for patch in bp['boxes']: 
    patch.set_facecolor(color_list[i])
    i+=1


pylab.xscale('log')
locs, dummy_labels = pylab.yticks()  
pylab.yticks(locs, labels, fontsize=9)   

pylab.savefig('%s/fixations_per_pathway_variable_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)




########################
# HEATMAPS             #
########################

###################
# pi per pathway  #
###################

min_number_genes=10
min_passed_sites_per_person=100

data_dict={}
# load the data into data_dict, where keys are species_name, values are mean_pi
for species_name in species_names:
    if species_name in files.keys():
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        avg_pi_per_pathway_core=files[species_name]['avg_pi_per_pathway_core']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        passed_sites_per_pathway_core=files[species_name]['passed_sites_per_pathway_core']
        num_genes_per_pathway_core=files[species_name]['num_genes_per_pathway_core']
        num_people_with_data_pathway_core=files[species_name]['num_people_with_data_pathway_core']
        mean_pi_dict={}
        mean_pi_dict['All core genes']=numpy.mean(avg_pi_matrix_core[same_sample_idxs])
        mean_pi_dict['All variable genes']=numpy.mean(avg_pi_matrix_variable[same_sample_idxs])
        for pathway in avg_pi_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(passed_sites_per_pathway_core[pathway][same_sample_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if num_genes_per_pathway_core[pathway]>=min_number_genes:
                    mean_pi=numpy.clip(numpy.mean(avg_pi_per_pathway_core[pathway][same_sample_idxs][high_num_sites_idx]),1e-7,1) 
                    if pathway!='':
                        mean_pi_dict[pathway]=mean_pi
                    else:
                        mean_pi_dict['Unannotated pathways']=mean_pi  
        data_dict[species_name]=mean_pi_dict

# convert data_dict into a pandas dataframe
df=pandas.DataFrame(data_dict)
pathway_order=[]
pathway_order.append('All core genes')
pathway_order.append('All variable genes')
pathway_order.append('Annotated pathways')
pathway_order.append('Unannotated pathways')
pathway_names=list(df.index)
for pathway in pathway_names:
    if pathway not in pathway_order:
        pathway_order.append(pathway)

df=df.reindex(pathway_order)

# use seaborn to plot a heat map of the pi values
pylab.figure(figsize=(24,8))
pylab.subplot(1,2,1)
sns.heatmap(df, cmap='RdYlGn_r',norm=LogNorm(vmin=1e-7, vmax=1e-2))
#pylab.savefig('%s/pi_per_pathway_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# Iterate through each species again and compute the rank order of each column. 
for species_name in species_names:   
    if species_name in files.keys():
        df[species_name]=df[species_name].rank(ascending=True)

# use seaborn to plot a heat map of the rank orders
#pylab.figure(figsize=(8,8))
pylab.subplot(1,2,2)
sns.heatmap(df, cmap='RdYlGn_r',yticklabels=False)
pylab.savefig('%s/pi_per_pathway_ranking_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)



################################
# fraction nonsyn per pathway  #
# core genes                   #
################################

min_number_genes=10
min_passed_sites_per_person=100

data_dict={}
# load the data into data_dict, where keys are species_name, values are mean_pi
for species_name in species_names:
    if species_name in files.keys():
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        fraction_nonsynonymous_per_pathway_core=files[species_name]['fraction_nonsynonymous_per_pathway_core']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_syn_non_core=files[species_name]['fixation_opportunities_per_pathway_syn_non_core']
        num_genes_per_pathway_syn_non_core=files[species_name]['num_genes_per_pathway_syn_non_core']
        num_people_with_data_per_pathway_fixations_core=files[species_name]['num_people_with_data_per_pathway_fixations_core']
        mean_pi_dict={}
        mean_pi_dict['All core genes']=numpy.mean(fraction_nonsynonymous_core[diff_subject_idxs])
        mean_pi_dict['All variable genes']=numpy.mean(fraction_nonsynonymous_variable[diff_subject_idxs])
        for pathway in fraction_nonsynonymous_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if num_genes_per_pathway_syn_non_core[pathway]>=min_number_genes:
                    mean_pi=numpy.mean(fraction_nonsynonymous_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                    if pathway!='':
                        mean_pi_dict[pathway]=mean_pi
                    else:
                        mean_pi_dict['Unannotated pathways']=mean_pi  
        data_dict[species_name]=mean_pi_dict

# convert data_dict into a pandas dataframe
df=pandas.DataFrame(data_dict)
pathway_order=[]
pathway_order.append('All core genes')
pathway_order.append('All variable genes')
pathway_order.append('Annotated pathways')
pathway_order.append('Unannotated pathways')
pathway_names=list(df.index)
for pathway in pathway_names:
    if pathway not in pathway_order:
        pathway_order.append(pathway)

df=df.reindex(pathway_order)

# use seaborn to plot a heat map of the pi values
pylab.figure(figsize=(24,8))
pylab.subplot(1,2,1)
sns.heatmap(df, cmap='RdYlGn_r',vmin=0.0, vmax=1)
#pylab.savefig('%s/pi_per_pathway_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# Iterate through each species again and compute the rank order of each column. 
for species_name in species_names:   
    if species_name in files.keys():
        df[species_name]=df[species_name].rank(ascending=True)

# use seaborn to plot a heat map of the rank orders
#pylab.figure(figsize=(8,8))
pylab.subplot(1,2,2)
sns.heatmap(df, cmap='RdYlGn_r',yticklabels=False)
pylab.savefig('%s/fraction_nonsyn_per_pathway_core_ranking_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


################################
# fraction nonsyn per pathway  #
# variable genes               #
################################

min_number_genes=5
min_passed_sites_per_person=100

data_dict={}
# load the data into data_dict, where keys are species_name, values are mean_pi
for species_name in species_names:
    if species_name in files.keys():
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        fraction_nonsynonymous_per_pathway_core=files[species_name]['fraction_nonsynonymous_per_pathway_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_syn_non_core=files[species_name]['fixation_opportunities_per_pathway_syn_non_variable']
        num_genes_per_pathway_syn_non_core=files[species_name]['num_genes_per_pathway_syn_non_variable']
        num_people_with_data_per_pathway_fixations_core=files[species_name]['num_people_with_data_per_pathway_fixations_variable']
        mean_pi_dict={}
        mean_pi_dict['All core genes']=numpy.mean(fraction_nonsynonymous_core[diff_subject_idxs])
        mean_pi_dict['All variable genes']=numpy.mean(fraction_nonsynonymous_variable[diff_subject_idxs])
        for pathway in fraction_nonsynonymous_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if num_genes_per_pathway_syn_non_core[pathway]>=min_number_genes:
                    mean_pi=numpy.mean(fraction_nonsynonymous_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                    if pathway!='':
                        mean_pi_dict[pathway]=mean_pi
                    else:
                        mean_pi_dict['Unannotated pathways']=mean_pi  
        data_dict[species_name]=mean_pi_dict

# convert data_dict into a pandas dataframe
df=pandas.DataFrame(data_dict)
pathway_order=[]
pathway_order.append('All core genes')
pathway_order.append('All variable genes')
pathway_order.append('Annotated pathways')
pathway_order.append('Unannotated pathways')
pathway_names=list(df.index)
for pathway in pathway_names:
    if pathway not in pathway_order:
        pathway_order.append(pathway)

df=df.reindex(pathway_order)

# use seaborn to plot a heat map of the pi values
pylab.figure(figsize=(24,8))
pylab.subplot(1,2,1)
sns.heatmap(df, cmap='RdYlGn_r')
#pylab.savefig('%s/pi_per_pathway_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# Iterate through each species again and compute the rank order of each column. 
for species_name in species_names:   
    if species_name in files.keys():
        df[species_name]=df[species_name].rank(ascending=True)

# use seaborn to plot a heat map of the rank orders
#pylab.figure(figsize=(8,8))
pylab.subplot(1,2,2)
sns.heatmap(df, cmap='RdYlGn_r',yticklabels=False)
pylab.savefig('%s/fraction_nonsyn_per_pathway_variable_ranking_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)



################################
# fraction nonsyn per pathway  #
# core and variable genes      #
################################

min_number_genes=5
min_passed_sites_per_person=100

data_dict={}
# load the data into data_dict, where keys are species_name, values are mean_pi
for species_name in species_names:
    if species_name in files.keys():
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        fraction_nonsynonymous_per_pathway_core=files[species_name]['fraction_nonsynonymous_per_pathway_core_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        fixation_opportunities_per_pathway_syn_non_core=files[species_name]['fixation_opportunities_per_pathway_syn_non_core_variable']
        num_genes_per_pathway_syn_non_core=files[species_name]['num_genes_per_pathway_syn_non_core_variable']
        num_people_with_data_per_pathway_fixations_core=files[species_name]['num_people_with_data_per_pathway_fixations_core_variable']
        mean_pi_dict={}
        mean_pi_dict['All core genes']=numpy.mean(fraction_nonsynonymous_core[diff_subject_idxs])
        mean_pi_dict['All variable genes']=numpy.mean(fraction_nonsynonymous_variable[diff_subject_idxs])
        for pathway in fraction_nonsynonymous_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if num_genes_per_pathway_syn_non_core[pathway]>=min_number_genes:
                    mean_pi=numpy.mean(fraction_nonsynonymous_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                    if pathway!='':
                        mean_pi_dict[pathway]=mean_pi
                    else:
                        mean_pi_dict['Unannotated pathways']=mean_pi  
        data_dict[species_name]=mean_pi_dict

# convert data_dict into a pandas dataframe
df=pandas.DataFrame(data_dict)
df.index.name='pathways'
# get the mean value across species:
df['mean']=df.mean(axis=1)
df=df.sort('mean')
pathway_names=list(df.index)
pathway_order=[]
pathway_order.append('All core genes')
pathway_order.append('All variable genes')
pathway_order.append('Annotated pathways')
pathway_order.append('Unannotated pathways')
for pathway in pathway_names:
    if pathway not in pathway_order:
        pathway_order.append(pathway)

df=df.reindex(pathway_order)
del df['mean']

# drop Acidaminococcus_intestini_54097 because it messes up coloring
#df=df.drop('Acidaminococcus_intestini_54097',axis=1)

# make a boxplot where each species is a data point. Order the pathways based on mean fraction nonsyn
pylab.figure(figsize=(8,8))
df.T.boxplot(vert=False)
pylab.savefig('%s/fraction_nonsyn_per_pathway_boxplot_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)

#sns.pairplot(df)

# plot a scatter plot of every pair of points
# randomly sample n species to plot
n=8
random_species_names=random.sample(list(df.columns.values),n)

pylab.figure(figsize=(30,30))
plot_no=1
for i in range(0, n):
    for j in range(0,n):        
        pylab.subplot(n,n,plot_no) 
        x=list(df[random_species_names[i]])
        y=list(df[random_species_names[j]])
        pylab.scatter(x,y, color='red',marker='.')
        pylab.xlabel(random_species_names[i])
        pylab.ylabel(random_species_names[i])
        pylab.xlim(0,0.5)
        pylab.ylim(0,0.5)
        #ax=df.plot(x=random_species_names[i], y=random_species_names[j], style='o', xlim=[0,0.8], ylim=[0,0.8],legend=False,subplots=True)
        #ax.set_xlabel(random_species_names[i])
        #ax.set_ylabel(random_species_names[j])
        plot_no+=1

pylab.savefig('%s/fraction_nonsyn_per_pathway_core_variable_pairplot_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)

# plot the spearman correlation coefficients between every pair of species
df = df[species_order]
spearman_corr=df.corr(method='spearman')

# plot a heatmap of the spearman correlation coefficients
pylab.figure(figsize=(15,15))
sns.heatmap(spearman_corr, cmap='RdYlBu_r')
pylab.savefig('%s/fraction_nonsyn_per_pathway_core_variable_spearman_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# use seaborn to plot a heat map of the pi values
pylab.figure(figsize=(24,8))
pylab.subplot(1,2,1)
sns.heatmap(df, cmap='RdYlBu_r')

# Iterate through each species again and compute the rank order of each column. 
for species_name in species_names:   
    if species_name in df.columns:
        df[species_name]=df[species_name].rank(ascending=True)

# use seaborn to plot a heat map of the rank orders
pylab.subplot(1,2,2)
sns.heatmap(df, cmap='RdYlBu_r',yticklabels=False)
pylab.savefig('%s/fraction_nonsyn_per_pathway_core_variable_ranking_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


################################
# Total fixations per pathway  #
# core genes                   #
################################

min_number_genes=10
min_passed_sites_per_person=100

data_dict={}
# load the data into data_dict, where keys are species_name, values are mean_pi
for species_name in species_names:
    if species_name in files.keys():
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        fixation_opportunities_per_pathway_all_core=files[species_name]['fixation_opportunities_per_pathway_all_core']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        dtot_per_pathway_core=files[species_name]['dtot_per_pathway_core']
        num_genes_per_pathway_tot_core=files[species_name]['num_genes_per_pathway_tot_core']
        num_people_with_data_per_pathway_tot_core=files[species_name]['num_people_with_data_per_pathway_tot_core']
        mean_pi_dict={}
        mean_pi_dict['All core genes']=numpy.mean(dtot_core[diff_subject_idxs])
        mean_pi_dict['All variable genes']=numpy.mean(dtot_variable[diff_subject_idxs])
        for pathway in dtot_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_all_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if num_genes_per_pathway_tot_core[pathway]>=min_number_genes:
                    mean_pi=numpy.mean(dtot_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                    if pathway!='':
                        mean_pi_dict[pathway]=mean_pi
                    else:
                        mean_pi_dict['Unannotated pathways']=mean_pi  
        data_dict[species_name]=mean_pi_dict

# convert data_dict into a pandas dataframe
df=pandas.DataFrame(data_dict)
pathway_order=[]
pathway_order.append('All core genes')
pathway_order.append('All variable genes')
pathway_order.append('Annotated pathways')
pathway_order.append('Unannotated pathways')
pathway_names=list(df.index)
for pathway in pathway_names:
    if pathway not in pathway_order:
        pathway_order.append(pathway)

df=df.reindex(pathway_order)

# use seaborn to plot a heat map of the pi values
pylab.figure(figsize=(24,8))
pylab.subplot(1,2,1)
sns.heatmap(df, cmap='RdYlGn_r')
#pylab.savefig('%s/pi_per_pathway_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# Iterate through each species again and compute the rank order of each column. 
for species_name in species_names:   
    if species_name in files.keys():
        df[species_name]=df[species_name].rank(ascending=True)

# use seaborn to plot a heat map of the rank orders
#pylab.figure(figsize=(8,8))
pylab.subplot(1,2,2)
sns.heatmap(df, cmap='RdYlGn_r',yticklabels=False)
pylab.savefig('%s/fixations_per_pathway_core_ranking_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)




################################
# Total fixations per pathway  #
# variable genes               #
################################

min_number_genes=5
min_passed_sites_per_person=100

data_dict={}
# load the data into data_dict, where keys are species_name, values are mean_pi
for species_name in species_names:
    if species_name in files.keys():
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        fixation_opportunities_per_pathway_all_core=files[species_name]['fixation_opportunities_per_pathway_all_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        dtot_per_pathway_core=files[species_name]['dtot_per_pathway_variable']
        num_genes_per_pathway_tot_core=files[species_name]['num_genes_per_pathway_tot_variable']
        num_people_with_data_per_pathway_tot_core=files[species_name]['num_people_with_data_per_pathway_tot_variable']
        mean_pi_dict={}
        mean_pi_dict['All core genes']=numpy.mean(dtot_core[diff_subject_idxs])
        mean_pi_dict['All variable genes']=numpy.mean(dtot_variable[diff_subject_idxs])
        for pathway in dtot_per_pathway_core:
            #check which people have enough data. Otherwise don't add this.
            high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_all_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)
            if len(high_num_sites_idx)>0:
                if num_genes_per_pathway_tot_core[pathway]>=min_number_genes:
                    mean_pi=numpy.mean(dtot_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx]) 
                    if pathway!='':
                        mean_pi_dict[pathway]=mean_pi
                    else:
                        mean_pi_dict['Unannotated pathways']=mean_pi  
        data_dict[species_name]=mean_pi_dict

# convert data_dict into a pandas dataframe
df=pandas.DataFrame(data_dict)
pathway_order=[]
pathway_order.append('All core genes')
pathway_order.append('All variable genes')
pathway_order.append('Annotated pathways')
pathway_order.append('Unannotated pathways')
pathway_names=list(df.index)
for pathway in pathway_names:
    if pathway not in pathway_order:
        pathway_order.append(pathway)

df=df.reindex(pathway_order)

# use seaborn to plot a heat map of the pi values
pylab.figure(figsize=(24,8))
pylab.subplot(1,2,1)
sns.heatmap(df, cmap='RdYlGn_r')
#pylab.savefig('%s/pi_per_pathway_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)


# Iterate through each species again and compute the rank order of each column. 
for species_name in species_names:   
    if species_name in files.keys():
        df[species_name]=df[species_name].rank(ascending=True)

# use seaborn to plot a heat map of the rank orders
#pylab.figure(figsize=(8,8))
pylab.subplot(1,2,2)
sns.heatmap(df, cmap='RdYlGn_r',yticklabels=False)
pylab.savefig('%s/fixations_per_pathway_variable_ranking_heatmap_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)

