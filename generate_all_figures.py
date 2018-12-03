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
os.system('python plot_supplemental_sfss.py')

##################################
#
# Figure 2: between-host divergence across prevalent species
# (figure_2.pdf)
#
# Supplemental Figure: Raw # of SNV and gene differences 
#                      between closely related strains
# (supplemental_closely_related_differences.pdf)
#
##################################
os.system('python plot_figure_2.py')

##################################
#
# Figure 3: dN/dS vs dS
# (figure_3.pdf)
#
# Supplemental Figure: Singleton dN/dS vs dS
# (supplemental_singleton_dNdS.pdf)
os.system('python plot_figure_3.py')

###################################
#
# Figure 4: Signatures of recombination across hosts
# (figure_4.pdf)
#
# Supplemental Figure: decay of LD in 3 example species
#
# Supplemental Figure: estimates of r/mu across species
# (supplemental_rbymu.pdf)
##################################
os.system('python plot_figure_4.py')

##################################
#
# Supplemental Figure: clade / geographic structure 
# (supplemental_clade_city_correlation.pdf)
#
##################################
os.system('python plot_supplemental_manual_clade_fst.py') 

##################################
#
# Supplemental Figure: QP status over time
# (supplemental_temporal_haploid.pdf)
#
##################################
os.system('plot_supplemental_qp_over_time.py')

##################################
#
# Supplemental Figure: Marker SNV sharing
# (supplemental_marker_sharing.pdf)
#
##################################
os.system('plot_supplemental_marker_snv_sharing.py')


##################################
#
# Figure 5: within-host changes across prevalent species
# (figure_5.pdf)
#
##################################
os.system('python plot_figure_5.py')



