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

############
# color map
def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color


#############

#species_names=['Alistipes_putredinis_61533','Bacteroides_uniformis_57318']
species_names = parse_midas_data.parse_good_species_list()

files={}

for i in range(0, 30):
    species_name=species_names[i]
    if (os.path.exists(os.path.expanduser('~/tmp_intermediate_files/kegg_pi_%s.npz' % species_name))):
        files[species_name]=numpy.load(os.path.expanduser('~/tmp_intermediate_files/kegg_pi_%s.npz' % species_name))

colors=['#a1d99b','#c994c7']  



######################################################################
# plot a comparison of pi for core vs variable genes accross species
#####################################################################
# try ordering species based on median piS in core
species_list=[]
mean_avg_pi=[]
for species_name in species_names:
    if species_name in files.keys():
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        tmp1=numpy.array(same_sample_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(same_sample_idxs[1], dtype=numpy.int32)
        same_sample_idxs=(tmp1,tmp2)
        avg_pi=numpy.median(avg_pi_matrix_core[same_sample_idxs])
        species_list.append(species_name)
        mean_avg_pi.append(avg_pi)

sorted_species=[x for (y,x) in sorted(zip(mean_avg_pi,species_list))]


pylab.figure(figsize=(6,12))
pylab.xlabel('Pi/bp')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1e-1)

data=[]    
labels=[]   
for species_name in species_names:
    if species_name in files.keys():
#for species_name in sorted_species:
        avg_pi_matrix_core=files[species_name]['avg_pi_matrix_core']
        avg_pi_matrix_variable=files[species_name]['avg_pi_matrix_variable']
        same_sample_idxs=files[species_name]['same_sample_idxs']
        tmp1=numpy.array(same_sample_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(same_sample_idxs[1], dtype=numpy.int32)
        same_sample_idxs=(tmp1,tmp2)
        data.append(avg_pi_matrix_core[same_sample_idxs]) 
        labels.append(species_name + '_core, n=' + str(len(same_sample_idxs[0]))) 
        data.append(avg_pi_matrix_variable[same_sample_idxs]) 
        labels.append(species_name + '_variable, n=' + str(len(same_sample_idxs[0]))) 

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

# try ordering species based on median fraction nonsyn in core
species_list=[]
mean_avg_pi=[]
for species_name in species_names:
    if species_name in files.keys():
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        avg_pi=numpy.median(fraction_nonsynonymous_core[diff_subject_idxs])
        species_list.append(species_name)
        mean_avg_pi.append(avg_pi)

sorted_species=[x for (y,x) in sorted(zip(mean_avg_pi,species_list))]


pylab.figure(figsize=(6,12))
pylab.xlabel('Fraction nonsynonymous fixations')
pylab.ylabel("Species/gene type")
pylab.xlim(0,1)

data=[]    
labels=[]   
for species_name in species_names:
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
# plot a comparison of total fixations for core vs variable genes accross species
#####################################################################

# try ordering species based on median dtot_core
species_list=[]
mean_avg_pi=[]
for species_name in species_names:
    if species_name in files.keys():
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        avg_pi=numpy.median(dtot_core[diff_subject_idxs])
        species_list.append(species_name)
        mean_avg_pi.append(avg_pi)

sorted_species=[x for (y,x) in sorted(zip(mean_avg_pi,species_list))]


pylab.figure(figsize=(6,12))
pylab.xlabel('Total fixations')
pylab.ylabel("Species/gene type")
pylab.xlim(1e-7,1e-1)

data=[]    
labels=[]   
for species_name in species_names:
    if species_name in files.keys():
#for species_name in sorted_species:
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        same_sample_idxs=files[species_name]['same_sample_idxs']
        data.append(dtot_core[diff_subject_idxs]) 
        labels.append(species_name + '_core, n=' + str(len(same_sample_idxs[0]))) 
        data.append(dtot_variable[diff_subject_idxs]) 
        labels.append(species_name + '_variable, n=' + str(len(same_sample_idxs[0]))) 

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

pylab.figure(figsize=(8,8))
pylab.xlabel('Total fixations')
pylab.ylabel("Fraction nonsynonymous")
pylab.xlim(1e-7,1e-1)
pylab.ylim(0,1)

data=[]    
labels=[]   
for species_name in species_names:
    if species_name in files.keys():
#for species_name in sorted_species:
        fraction_nonsynonymous_core=files[species_name]['fraction_nonsynonymous_core']
        fraction_nonsynonymous_variable=files[species_name]['fraction_nonsynonymous_variable']
        dtot_core=files[species_name]['dtot_core']
        dtot_variable=files[species_name]['dtot_variable']
        diff_subject_idxs=files[species_name]['diff_subject_idxs']
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        same_sample_idxs=files[species_name]['same_sample_idxs']
        tmp1=numpy.array(same_sample_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(same_sample_idxs[1], dtype=numpy.int32)
        same_sample_idxs=(tmp1,tmp2)
        pylab.semilogx(dtot_core[diff_subject_idxs],fraction_nonsynonymous_core[diff_subject_idxs], 'r.', label='core')
        pylab.semilogx(dtot_variable[diff_subject_idxs],fraction_nonsynonymous_variable[diff_subject_idxs], 'b.', label='variable')

#pylab.legend(loc='lower right',frameon=False)
pylab.savefig('%s/core_vs_variable_genes_fraction_nonsyn_vs_fixations_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)




################################################
# repeat except color each species differently #
################################################

species_names_Bacteroides=['Bacteroides_uniformis_57318', 'Bacteroides_vulgatus_57955','Bacteroides_fragilis_54507','Bacteroides_thetaiotaomicron_56941','Bacteroides_ovatus_58035', 'Bacteroides_massiliensis_44749','Bacteroides_cellulosilyticus_58046','Bacteroides_finegoldii_57739','Bacteroides_stercoris_56735','Bacteroides_caccae_53434','Bacteroides_xylanisolvens_57185']

species_names_Alistipes=['Alistipes_finegoldii_56071','Alistipes_onderdonkii_55464','Alistipes_putredinis_61533', 'Alistipes_shahii_62199']




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
        tmp1=numpy.array(diff_subject_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(diff_subject_idxs[1], dtype=numpy.int32)
        diff_subject_idxs=(tmp1,tmp2)
        same_sample_idxs=files[species_name]['same_sample_idxs']
        tmp1=numpy.array(same_sample_idxs[0], dtype=numpy.int32)
        tmp2=numpy.array(same_sample_idxs[1], dtype=numpy.int32)
        same_sample_idxs=(tmp1,tmp2)
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
