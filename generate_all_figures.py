import os
##################################
#
# Plots the figures for our paper
#
##################################

##################################
#
# Figure 1: within-host diversity
# (figure_1.pdf)
#
# Supplemental Figure: average fixed differences vs polymorphism
# (supplemental_avg_distance.pdf)
#
# Supplemental Figure: correlation between diversity and fraction CP
# (supplemental_haploid_correlation.pdf)
#
# Supplemental Figure: distribution of CP samples per host. 
# (supplemental_haploid_distribution_fig)
#
# Supplemental Figure: temporal CP distribution
# (supplemental_temporal_haploid.pdf)
#
# Supplemental Figure: distribution of gene copynum for 4 example samples
# (supplemental_copynum_distributions.pdf)
#
##################################
os.system('python plot_figure_1.py')

##################################
#
# Figure 2: between-host divergence in 
# B. vulgatus and B. stercoris
# (figure_2.pdf)
#
##################################
os.system('python plot_figure_2.py')

##################################
#
# Figure 3: between-host divergence across prevalent species
# (figure_3.pdf)
#
##################################
os.system('python plot_figure_3.py')

##################################
#
# Figure 4: decay of linkage disequilibrium
# (figure_4.pdf)
# Supplemental Figure: estimate of r/mu
# (supplemental_rbymu.pdf)
##################################
os.system('python plot_figure_4.py')

##################################
#
# Figure 5: within-host changes in B. vulgatus
# (figure_5.pdf)
#
##################################
os.system('python plot_figure_5.py')

##################################
#
# Figure 6: within-host changes across prevalent species
# (figure_6.pdf)
#
##################################
os.system('python plot_figure_6.py')

############
#
# Figure 7 is hand-drawn schematic
#
############

##################################
#
# Supplemental Figure: clade / geographic structure 
# (supplemental_clade_city_correlation.pdf)
#
##################################
os.system('python plot_supplemental_clade_city_correlation.py') 

##################################
#
# Supplemental figure: singletons as function of min divergence
# (supplemental_singleton_distribution.pdf)
#
##################################
os.system('python plot_singleton_distribution.py') 


##################################
#
# Supplemental Figure: triplet timecourses
# (supplemental_triplet_trajectories.png)
#
##################################
# commented out to save time. uncomment when need to regenerate
#os.system('python plot_supplemental_triplets.py') 


##################################
#
# Supplemental Figure: Examples of clonal vs local sweep in B. vulgatus
#
##################################
# commented out to save time. uncomment when need to regenerate
#os.system('python plot_supplemental_clonal_local.py')
