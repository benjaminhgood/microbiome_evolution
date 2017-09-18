import os
##################################
#
# Plots the figures for our paper
#
##################################

##################################
#
# Figure 1: within-host diversity
# Supplemental Figure: temporal CPS distribution
#
##################################
os.system('python plot_figure_1.py')

##################################
#
# Figure 2: between-host divergence in B. vulgatus
#
##################################
os.system('python plot_figure_2.py')

##################################
#
# Figure 3: between-host divergence across prevalent species
#
##################################
os.system('python plot_figure_3.py')

##################################
#
# Figure 4: decay of linkage disequilibrium
# Supplemental Figure: estimate of r/mu
#
##################################
os.system('python plot_figure_4.py')

##################################
#
# Figure 5: within-host changes in B. vulgatus
#
##################################
os.system('python plot_figure_5.py')

##################################
#
# Figure 6: within-host changes across prevalent species
#
##################################
os.system('python plot_figure_7.py')

##################################
#
# Figure 7: Examples of clonal vs local sweep in B. vulgatus
#
##################################
os.system('python plot_supplemental_joint_sfs.py')




