args = commandArgs(TRUE)
species_list_file=args[1]
#species_list_file='~/tmp_intermediate_files/species_names.txt'

# read in the data
d=read.table(species_list_file, sep='\t')

# create a dictionary for labels
species_list=as.vector(d[,2])
species_full_name=as.vector(d[,1])
names(species_full_name)=species_list

# set the working directory
setwd('~/ben_nandita_hmp_analysis')

# install ape:
# install.packages("ape")

# load ape:
library(ape)

tree=read.tree(file='~/ben_nandita_hmp_data/midas_db/species_tree.newick')
tree_tips=tree$tip

#iterate through tree_tips to drop the tips that are not in the species_list
for (i in c(1:length(tree_tips))){
    species_id=tree_tips[i]
    if ((species_id %in% species_list) == FALSE){
	    tree=drop.tip(tree,c(species_id))
	}
    else {
	    print(species_id)
	}
}

# iterate through and rename the tips
# also output to a file the order of the tip labels.
tree_tips=tree$tip
for (i in c(1:length(tree_tips))){
    tree$tip.label[i]=species_full_name[tree$tip.label[i]]
}

# write the tip order to a file so that I can read it back into python. 
write(tree$tip, file = "~/tmp_intermediate_files/species_order.txt", append = FALSE, sep = "\n")

# print out the tree pdf
out_file=paste("~/ben_nandita_hmp_analysis/pruned_species_tree.pdf",sep='')
pdf(out_file,width=6,height=10, title = "Figure1",paper="special")   
plot(tree,show.tip.label = TRUE, main='Pruned species tree', cex=0.5)
add.scale.bar(cex = 0.7, font = 2)
dev.off()

